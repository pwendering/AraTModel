function pfba_problem = addPfbaConstraints(problem, minVarIdx)
%% pfba_problem = addPfbaConstraints(problem, minVarIdx)
% Add additional variables and constraints to minimize the sum of absolute
% fluxe values:
%       min sum(|v|) = sum(delta(+) + delta(-))
% s.t.
%       *previous constraints*
%       v = delta(+) - delta(-)
% Input
%   struct problem:                 Gurobi-style optimization problem
%   integer minVarIdx:              indices of variables whose sum should
%                                   be minimized (i.e. reactions)
% Output
%   struct pfba_problem:            updated optimization problem

nMinVar = numel(minVarIdx);

% create additional linear constraint matrix for minimization of absolute
% values
minConstMat = sparse(nMinVar, size(problem.A,2)+2*nMinVar);
minConstMat(:, minVarIdx) = eye(nMinVar);
deltaPosIdx = size(problem.A,2)+1:size(problem.A,2)+nMinVar;
minConstMat(:, deltaPosIdx) = -eye(nMinVar);
deltaNegIdx = size(problem.A,2)+nMinVar+1:size(problem.A,2)+2*nMinVar;
minConstMat(:, deltaNegIdx) = eye(nMinVar);

% update linear problem
pfba_problem = problem;
% ~~~~~~~~~~~~~~~~~~~~~~~ linear constraint matrix ~~~~~~~~~~~~~~~~~~~~~~ %
pfba_problem.A = [...
    problem.A sparse(size(problem.A,1), 2*nMinVar);...
    minConstMat];

% ~~~~~~~~~~~~~~~~~~~~~~~ linear objective vector ~~~~~~~~~~~~~~~~~~~~~~~ %
pfba_problem.obj = [zeros(size(problem.obj)); ones(2*nMinVar,1)];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ objective sense ~~~~~~~~~~~~~~~~~~~~~~~~~~ %
pfba_problem.modelsense = 'min';

% ~~~~~~~~~~~~~~~~~~~~ sense of linear constraints ~~~~~~~~~~~~~~~~~~~~~~ %
pfba_problem.sense = [problem.sense; repelem('=', nMinVar, 1)];

% ~~~~~~~~~~~~~~~~ right-hand side of linear constraints ~~~~~~~~~~~~~~~~ %
pfba_problem.rhs = [problem.rhs; zeros(nMinVar,1)];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~ constraint names ~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
pfba_problem.constrnames = [problem.constrnames;...
    cellfun(@(x)strcat('pfba_constr_', x), problem.varnames(minVarIdx), 'un', 0)];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ lower bounds ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
pfba_problem.lb = [problem.lb; zeros(2*nMinVar,1)];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ upper bounds ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
pfba_problem.ub = [problem.ub; repelem(1000,2*nMinVar,1)];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ variable names ~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
pfba_problem.varnames = [problem.varnames;...
    cellfun(@(x)strcat('delta_pos_', x), problem.varnames(minVarIdx), 'un', 0);...
    cellfun(@(x)strcat('delta_neg_', x), problem.varnames(minVarIdx), 'un', 0)];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~ variable types ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
pfba_problem.vtype = [problem.vtype; repelem('C', 2*nMinVar, 1)];

% ~~~~~~~~~~~~~~~~~~~~~~~~~ quadratic constraints ~~~~~~~~~~~~~~~~~~~~~~~ %

% quadratic part
pfba_problem.quadcon.Qc = [problem.quadcon.Qc sparse(size(problem.quadcon.Qc,1),2*nMinVar);...
    sparse(2*nMinVar, size(problem.quadcon.Qc,2)+2*nMinVar)];

% linear part
pfba_problem.quadcon.q = [problem.quadcon.q; zeros(2*nMinVar,1)];

end