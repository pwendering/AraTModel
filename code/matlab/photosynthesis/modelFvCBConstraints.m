function [solution,photParams,gurobiSolution,gurobiProblem] = ...
    modelFvCBConstraints(model,T,I,solverParams,nRxns,saBool,...
    saPercentage,saParameter)
%% modelFvCBConstraints(model,T,I,solverParams,nRxns,saBool,saPercentage,saParameter, gurobiProblem)
% Constraints derived from the FvCB model (Farquhar et al. 1980) and
% extensions are added to a metabolic model. The resulting quadratically-
% constraint optimization problem is then solved using the Gurobi solver.
% Finally, predicted growth and net CO2 assimilation rate are reported as
% well as temperature-adjusted parameters for the FvCB model.
% Input:
%   struct model:               metabolic model with additional fields
%                               added using the 'createTGEM' function
%   double scalar T:            temperature in Kelvin
%   double I:                   irradiance (umol/m2/s)
%   struct solverParams:        parameters for the Gurobi solver (tested
%                               with v9.1.1); if empty, default parameters
%                               are used as specified in 'solvePhotConstModel'
%   double nRxns:               number of reaction flux variables in the
%                               problem; it is assumed that the reaction
%                               indices go from 1:nRxns
%   logical saBool:             (optional) specifies whether or not a
%                               sensitivity analysis is performed (defaults
%                               to false)
%   double saPercentage:        (optional) if saBool is true, specifies the
%                               percentage for increase of the parameter
%                               given in saParameter at the given
%                               temperature
%   char saParameter:           (optional) name of the parameter that
%                               should be increased by saPercentage
% Output:
%   struct solution:            contains fields
%                       mu:     average of predicted growth rates from
%                               pool solutions [h^-1]
%                       A:      average of predicted net CO2 assimilation rate
%                               from pool solutions [umol gDW^-1 h^-1]
%                       mu_min: minimum growth rate over all pool solutions
%                       mu_max: maximum growth rate over all pool solutions
%                       A_min:  minium A over all pool solutions
%                       A_max:  maximum A over all pool solutions
%                       A_net:  net CO2 assimilation rate
%                               calculated by
%                               sum(fluxes cons. CO2) - sum(fluxes prod. CO2)
%                       x_step1:solution vector associated with the
%                               optimal solution to the first step
%                       x_step2:solution vector associated with the
%                               optimal solution to the second step
%                       objVal: average objective value; if the
%                               biomass reaction is the objective,
%                               this field contains the same value
%                               as the "mu" field and the unit will 
%                               be h^-1, otherwise  [mmol gDW^-1 h^-1]
%                       objVal_min: minimum objective value
%                                   [h^-1]/[mmol gDW^-1 h^-1]
%                       objVal_max: maximum objective value
%                                   [h^-1]/[mmol gDW^-1 h^-1]
%   struct photParams:          contains temperature-adjusted parameters of
%                               for the C3 photosynthesis model
%   struct gurobiSolution:      solution structure returned by the Gurobi
%                               solver (v9.1.1) after growth simulation
%   struct gurobiProblem:       optimization problem that was solved to
%                               obtain the solution in gurobiSolution

% define scaling factor for unit conversion between 
fvcb_scaling = config('fvcb_scaling');

if nargin < 5 || isempty(nRxns)
    fprintf('Number of reaction variables not given, using column dimension of S.\n')
    nRxns = size(model.S, 2);
end

if nargin < 6
    saBool = false;
elseif nargin >= 6 && nargin < 8 && saBool
    warning('Sensitivity analysis requested but additional values for percentage and parameter were not given, skipping...')
elseif nargin == 8 && saBool
    
    if saPercentage > 1 || saPercentage < -1
        saPercentage = saPercentage / 100;
    end
    
    if saPercentage > 1 || saPercentage < -1
        warning(['Invalid percentage for sensitivity analysis given: ' num2str(saPercentage) ', skipping...'])
        saBool=false;
    end
    
    if iscell(saParameter)
        saParameter = char(saParameter);
    end
end

% add CO2 uptake constraints derived from the C3 photosynthesis model
% (Farquhar et al. 1980)
if ~saBool
    [gurobiProblem, photParams] = addCO2UptakeConstraints(model, T, I);
else
    % if requested, perform a sensitivity analysis
    [gurobiProblem, photParams] = addCO2UptakeConstraints(model, T, I,...
        'saBool', saBool,...
        'saPercentage', saPercentage,...
        'saParameter', saParameter);
end

[solution,gurobiSolution] = solvePhotConstModel(gurobiProblem, solverParams, nRxns);

% change units of A from mmol gDW^-1 h^-1 to umol m^-2 s^-1
LMA = model.lma_model(T);
fvcb2fba = fvcb_scaling * 3600 / 1000 / LMA;
solution.A = solution.A / fvcb2fba;
solution.A_min = solution.A_min / fvcb2fba;
solution.A_max = solution.A_max / fvcb2fba;

%% calculate more exact net CO2 assimilation
% find all reactions that consume CO2
if isequal(solution.status ,'OPTIMAL')
    CO2_IDX = ismember(model.metFormulas,'CO2');
    exch_rxns = sum(model.S~=0)==1;
    co2_tr_rxns = sum(model.S(CO2_IDX,:)~=0)>1;
    co2_consuming = any(model.S(CO2_IDX,:)<0) & ~exch_rxns & ~co2_tr_rxns;
    co2_producing = any(model.S(CO2_IDX,:)>0) & ~exch_rxns & ~co2_tr_rxns;
    
    co2_net_assimilation = sum(solution.x_step2(co2_consuming)) - sum(solution.x_step2(co2_producing));
    % change units of A_net from mmol gDW^-1 h^-1 to umol m^-2 s^-1
    solution.A_net = co2_net_assimilation / fvcb2fba;
else
    solution.A_net = NaN;
end


end