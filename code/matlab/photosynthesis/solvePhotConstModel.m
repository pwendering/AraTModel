function [solution,gurobiSolution] = solvePhotConstModel(qcp,params,nRxns)
%% [mu,A,solution] = solvePhotConstModel(qcp,params)
% Function to solve a FBA or an enzyme-constaint optimization problem with
% additional constraints on CO2 uptake and photosynthesis (see
% 'addCO2UptakeConstraints'). To this end, a quadratically-constraint
% optimization problem is solved using the Gurobi solver.
%
% Input:
%   struct qcp:             quadratically-constraint optimization problem
%                           formatted for the Gurobi solver (tested with
%                           v9.1.1) as returned from the 'addCO2UptakeConstraints'
%                           function
%   struct params:          (optional) parameters to be passed to the
%                           Gurobi solver
%                           By default, the following settings are applied:
%                           params.NonConvex = 2;
%                           params.OutputFlag = 0;
%                           params.NumericFocus = 3;
%                           params.MIPFocus = 3;
%                           params.PoolSearchMode = 2;
%                           params.PoolSolutions = 20;
%   double nRxns:           number of reaction flux variables in the
%                           problem; it is assumed that the reaction
%                           indices go from 1:nRxns
% Output:
%   struct solution:        contains fields
%                           mu:     average of predicted growth rates from
%                                   pool solutions [h^-1]
%                           A:      average of predicted net CO2 assimilation rate
%                                   from pool solutions [umol gDW^-1 h^-1]
%                           mu_min: minimum growth rate over all pool solutions
%                           mu_max: maximum growth rate over all pool solutions
%                           A_min:  minium A over all pool solutions
%                           A_max:  maximum A over all pool solutions
%                           x_step1:solution vector associated with the
%                                   optimal solution to the first step
%                           x_step2:solution vector associated with the
%                                   optimal solution to the second step
%                           objVal: average objective value; if the
%                                   biomass reaction is the objective,
%                                   this field contains the same value
%                                   as the "mu" field and the unit will 
%                                   be h^-1, otherwise  [mmol gDW^-1 h^-1]
%                           objVal_min: minimum objective value
%                                   [h^-1]/[mmol gDW^-1 h^-1]
%                           objVal_max: maximum objective value
%                                   [h^-1]/[mmol gDW^-1 h^-1]
%   struct gurobiSolution:  solution structure returned by the Gurobi
%                           solver

% set parameters to default settings if not specified
if nargin < 2 || isempty(params) || numel(fieldnames(params)) == 0
    params = config('solver_params');
end

% number of reactions in the model
A_IDX = find(ismember(qcp.varnames,'A'));
% photon uptake reaction index
ABS_IDX = find(ismember(qcp.varnames,'PHOTON_UPTAKE'));
% biomass reaction index
BIO_IDX = find(strcmp(qcp.varnames,'GROWTH'));
% QCP objective indices
OBJ_IDX = find(qcp.obj);
BIO_OBJ_FLAG = isequal(OBJ_IDX,BIO_IDX);

% function to validate the solution returned by the Gurobi solver
validateSolution = @(sol)isequal(sol.status, 'OPTIMAL');

% initialize solution structure
solution = struct('mu',[],'mu_min',[],'mu_max',[],...
    'A',[],'A_min',[],'A_max',[],'x_step1',[],'x_step2', [], 'status', '');

%% optimize for maximum growth
qcp.obj = zeros(size(qcp.A,2),1);
qcp.obj(BIO_IDX) = 1;

gurobiSolution = gurobi(qcp,params);
solution.status = gurobiSolution.status;

if validateSolution(gurobiSolution)
    
    numSol = length(gurobiSolution.pool);
    tmpMu = nan(numSol,1);
    for j=1:numSol
        tmpMu(j) = gurobiSolution.pool(j).xn(BIO_IDX);
    end
    
    solution.x_step1 = gurobiSolution.x(1:nRxns+2);
    solution.mu = mean(tmpMu);
    solution.mu_min = min(tmpMu);
    solution.mu_max = max(tmpMu);
    
    fprintf('%s - Growth optimization successful!\n', gurobiSolution.status)
else
    warning(['First problem has no valid solution, status: ' gurobiSolution.status])
    return
end

%% minimize the sum of fluxes at optimal growth
qcp.lb(BIO_IDX) = 0.99*gurobiSolution.x(BIO_IDX);
minObjIdx = 1:nRxns;
qcp = addPfbaConstraints(qcp, minObjIdx);

% solve problem
gurobiSolution = gurobi(qcp,params);
solution.status = gurobiSolution.status;

if validateSolution(gurobiSolution)
    numSol = length(gurobiSolution.pool);
    tmpA = nan(numSol,1);
    for j=1:numSol
        tmpA(j) = gurobiSolution.pool(j).xn(A_IDX);
    end
    solution.A = mean(tmpA);
    solution.A_min = min(tmpA);
    solution.A_max = max(tmpA);
    
    solution.x_step2 = gurobiSolution.x(1:nRxns+2); % reaction and variables A and Z
    
    fprintf('%s - pFBA successful!\n', gurobiSolution.status)
    
else
    warning(['Second problem has no valid solution, status: ' gurobiSolution.status])
    return
end

if ~BIO_OBJ_FLAG
    
    % fix minimum photon uptake at optimal growth
    qcp.ub(ABS_IDX) = 1.01*gurobiSolution.x(ABS_IDX);
    qcp.obj = zeros(size(qcp.A,2),1);
    qcp.obj(OBJ_IDX) = 1;
    qcp.modelsense = 'max';
    
    % solve problem
    gurobiSolution = gurobi(qcp,params);
    solution.status = gurobiSolution.status;
    
    if validateSolution(gurobiSolution)
        numSol = length(gurobiSolution.pool);
        tmpObj = nan(numSol,1);
        for j=1:numSol
            tmpObj(j) = sum(gurobiSolution.pool(j).xn(OBJ_IDX));
        end
        solution.objVal = mean(tmpObj);
        solution.objVal_min = min(tmpObj);
        solution.objVal_max = max(tmpObj);
        
        fprintf('%s - Third optimization successful!\n', gurobiSolution.status)
    else
        warning(['Third problem has no valid solution, status: ' gurobiSolution.status])
        return
    end
else
    solution.objVal = solution.mu;
    solution.objVal_min = solution.mu_min;
    solution.objVal_max = solution.mu_max;
end