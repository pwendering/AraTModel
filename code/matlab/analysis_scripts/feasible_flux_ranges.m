% Determine feasible flux ranges upon temperature changes

clear; clc

%% 1) Load model and specify conditions

% read model
model_file = config('tgemFile');
tmp = load(model_file);
model = tmp.TGEM;
clear tmp

% add RubisCO activase module
model = addRCA(model);

% correct typo in subsystem name
model.subSystems = cellfun(@(x)strrep(x, 'actyl-coenzyme A synthesis',...
    'acetyl-coenzyme A synthesis'),...
    model.subSystems,'un',0);

% block uptake of H2S and chloroplastic alanine transaminase
model.ub(findRxnIDs(model, {'Im_H2S', 'Im_H2S_REV'})) = 0;
model.ub(findRxnIDs(model, {'AlaTA_h', 'AlaTA_h_REV'})) = 0;
model.ub(findRxnIDs(model, {'Bio_AA', 'Bio_CLim', 'Bio_NLim'})) = 0;

% correct THF transport reaction
model.S(findMetIDs(model, 'H_h'), findRxnIDs(model, {'Tr_THF_h', 'Tr_THF_h_REV'})) = [-1 1];

% light intensitiy
irradiance = config('I'); % [umol m^-2 s^-1]

% temperature range
temperatures = celsius2kelvin(linspace(10, 40, 10));

% whether to adjust the total protein content to temperature
adj_ptot = true;

% reaction IDs
reactions = model.rxns(~startsWith(model.rxns, {'draw_', 'prot_'}));
reactions(ismember(reactions, 'Bio_opt')) = {'GROWTH'};

% solver parameters
solver_params = config('solver_params');

% name of workspace to store sampling results
c = clock;
datestring = sprintf('%d%02.0f%02.0f', c(1:3));
result_ws = ['flux_ranges_I_' num2str(irradiance) '_' datestring];
    
%% 2) Determine flux ranges at each temperature

% initialize matrix of mininum and maximum fluxes and ranges
fs_ranges = nan(numel(reactions), numel(temperatures));
rgr = nan(1, numel(temperatures));

for i = 1:numel(temperatures)
    % get temperature-constrained optimization problem
    [tmp_sol, ~, ~, tmp_qcp] = simulateTempEffects(model, irradiance,...
        'tempRange', temperatures(i),...
        'adjPtotFlag', adj_ptot);
    rgr(i) = tmp_sol.mu_max;

    nrxns = numel(reactions);
    for j = 1:numel(reactions)
        
        if j > 100 && mod(j, 100) == 1
            fprintf('Done with %d reactions (%d%%)\n', j-1, 100*round((j-1)/nrxns, 2))
        end
        
        % reset objective
        tmp_qcp.obj(:) = 0;
        
        % minimization
        tmp_qcp.obj(ismember(tmp_qcp.varnames, reactions(j))) = 1;
        tmp_qcp.modelsense = 'min';
        sol = gurobi(tmp_qcp, solver_params);
        min_flux = sol.objval;
        
        % maximization
        tmp_qcp.modelsense = 'max';
        sol = gurobi(tmp_qcp, solver_params);
        max_flux = sol.objval;
        
        fs_ranges(j, i) = max_flux - min_flux;
    end
    
end

save(result_ws)