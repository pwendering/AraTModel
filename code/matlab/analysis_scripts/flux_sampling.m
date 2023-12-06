% Sample flux distributions at different temperatures
clear;clc

%% 1) Input variables
rng(42)
% number of samples
n_samples = 3e4;
% maximum number of tries
max_tries = 5*n_samples;
% percent of reactions to select for minimization
p_select = 0.01;
% minimum of optimal objective value to guarantee for all samplings
min_obj_pct = 0.9;

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

rxns = {'ATPase_h';'cplx5_m';'Tr_HPR1';...
    'GCEADH_c';'MalDH2_m';'MalDH4_h';'NDC1_1_h';'NDC1_2_h';'NDH1_h';...
    'NDH2_h';'Cytb6f1_h';'Cytb6f2_h';'AOX1A_m';'AOX4_h';'Tr_Asp_OAA_h';...
    'Tr_DIT1_OAA_h';'Tr_H_Na_h';'Tr_Pi_h';'Tr_Pyr_H_h';'Tr_Pyr_Mal_h';...
    'Tr_Pyr_Na_h';'Tr_THF_h';'Tr_Asp_Glu_m';'Tr_Asp_m';'Tr_Cit_iCit_m';...
    'Tr_DIC4_m';'Tr_KG_Mal_m';'Tr_PIC_m';'Tr_DTC5_m';'Tr_UCP_m';'Tr_PXN3_p';...
    'Tr_SCA_p';'Tr_iCit_p';'Si_H_c';'Si_H_h';'Si_H_m';'Si_H_p';'AspAT2_c';...
    'AspAT2_h';'AspAT2_m';'AspAT2_p';'CHRM_h';'KARI12_h';'KARI34_h';...
    'MetAdT_c';'NDA2_1_m';'NDA2_2_m';'PGR5PGRL11_h';'PGR5PGRL12_h';'iCitL_p'};

subsystems = {'', 'oxidative phosphorylation', 'transport', 'photorespiration, serine synthesis',...
    'tricarboxylic acid cycle, glyoxylate cycle, gluconeogenesis',...
    'gluconeogenesis', 'phylloquinol biosynthesis', 'phylloquinol biosynthesis',...
    'cyclic electron flow', 'cyclic electron flow', 'light reactions, cyclic electron flow',...
    'light reactions, cyclic electron flow', 'alternative respiration',...
    'alternative respiration', 'transport', 'transport', 'transport', 'transport',...
    'transport', 'transport', 'transport', 'transport', 'transport', 'transport',...
    'transport', 'transport', 'transport', 'transport', 'transport', 'transport',...
    'transport', 'transport', 'transport', 'proton sink/source', 'proton sink/source',...
    'proton sink/source', 'proton sink/source',...
    'aspartate degradation, aspartate synthesis, glutamate degradation',...
    'aspartate degradation, aspartate synthesis, glutamate degradation',...
    'aspartate degradation, aspartate synthesis, glutamate degradation',...
    'aspartate degradation, aspartate synthesis, glutamate degradation',...
    'phenylalanine synthesis, tyrosine synthesis', 'valine synthesis',...
    'valine synthesis', 'S-adenosyl-L-methionine salvage, methionine degradation',...
    '', '', 'cyclic electron flow, photoprotection',...
    'cyclic electron flow, photoprotection', 'glyoxylate cycle'};

for i = 1:numel(rxns)
    model.subSystems(startsWith(model.rxns, rxns(i))) = subsystems(i);
end

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

% reaction IDs
reactions = model.rxns(~startsWith(model.rxns, {'draw_', 'prot_'}));
reactions(ismember(reactions, 'Bio_opt')) = {'GROWTH'};
reactions(ismember(reactions, 'Im_hnu')) = {'PHOTON_UPTAKE'};
reactions(ismember(reactions, 'Im_H2O_REV')) = {'TRANSPIRATION'};
n_rxns = numel(reactions);

solver_params = config('solver_params');
solver_params = rmfield(solver_params, {'PoolSearchMode', 'PoolSolutions'});

% name of workspace to store sampling results
c = clock;
datestr = sprintf('%d%02.0f%02.0f', c(1:3));
result_ws = ['sampling_workspace_I_' num2str(irradiance) '_' datestr];

save(result_ws, 'model', 'irradiance', 'min_obj_pct', 'temperatures');

%% 2) Obtain operational ranges at 90% of the optimal objective value
min_fluxes = nan(numel(reactions), numel(temperatures));
max_fluxes = nan(numel(reactions), numel(temperatures));

disp('Operational ranges')
for i = 1:numel(temperatures)
    
    [tmp_sol, ~, ~, tmp_qcp] = simulateTempEffects(model, irradiance,...
        'tempRange', temperatures(i));
    rxn_idx = ismember(tmp_qcp.varnames, reactions);
    
    % get minimum sum of fluxes from pFBA solution
    sum_flux = sum(tmp_sol.x_step2(1:numel(model.rxns)));
    
    % add constraint to add an upper bound on the sum of fluxes
    tmp_qcp.A = [tmp_qcp.A; [ones(1, numel(model.rxns)) 0 0]];
    tmp_qcp.rhs = [tmp_qcp.rhs; sum_flux];
    tmp_qcp.sense = [tmp_qcp.sense; '<'];
    tmp_qcp.constrnames = [tmp_qcp.constrnames; 'UB_sum_flux'];
    
    % set lower bound of objective value to 90% of the optimum
    fva_tmp_problem = tmp_qcp;
    fva_tmp_problem.lb(ismember(tmp_qcp.varnames, 'GROWTH')) = min_obj_pct*tmp_sol.mu_max;
    
    for j = 1:numel(reactions)
        
        if j > 100 && mod(j, 100) == 1
            fprintf('Done with %d reactions (%d%%)\n', j-1, 100*round((j-1)/n_rxns, 2))
        end
        
        % minimization
        fva_tmp_problem.obj(:) = 0;
        fva_tmp_problem.obj(ismember(fva_tmp_problem.varnames, reactions(j))) = 1;
        fva_tmp_problem.modelsense = 'min';
        sol = gurobi(fva_tmp_problem, solver_params);
        min_flux = sol.objval;
        
        % maximization
        fva_tmp_problem.modelsense = 'max';
        sol = gurobi(fva_tmp_problem, solver_params);
        max_flux = sol.objval;
        
        min_fluxes(j,i) = min_flux;
        max_fluxes(j,i) = max_flux;
    end
end

save(result_ws, 'min_fluxes', 'max_fluxes', '-append');
clear min_flux max_flux sol fva_tmp_problem tmp_sol tmp_qcp

%% 3) Perform sampling at each temperature
% create a random vector of fluxes within the operational ranges and
% project it onto the feasible space

samples = cell(1, numel(temperatures));
for i = 1:numel(temperatures)
    % rows = number of model variables + 2 (A and Z)
    samples{i} = sparse(size(model.S, 2)+2, n_samples);
end

n_failed = zeros(numel(temperatures), 1);
rgr = nan(1, numel(temperatures));
coverages = nan(1, numel(temperatures));
tic;
for i = 1:numel(temperatures)
    t1 = toc;
    % get QCP and optimal RGR value
    [tmp_sol, ~, ~, tmp_qcp] = simulateTempEffects(model, irradiance,...
        'tempRange', temperatures(i));
    rxn_idx = ismember(tmp_qcp.varnames, reactions);
    rgr(i) = tmp_sol.mu_max;
    
    % get minimum sum of fluxes from pFBA solution
    sum_flux = sum(tmp_sol.x_step2(1:numel(model.rxns)));
    
    % current operational range
    curr_min_fluxes = min_fluxes(:, i);
    curr_max_fluxes = max_fluxes(:, i);
    op_ranges = curr_max_fluxes - curr_min_fluxes;
    
    % max_delta = op_ranges + 1;
    max_delta = 1000*ones(n_rxns, 1);
    
    % First norm minimization
    sample_qcp = tmp_qcp;
    
    % add constraint to add an upper bound on the sum of fluxes
    tmp_qcp.A = [tmp_qcp.A; [ones(1, numel(model.rxns)) 0 0]];
    tmp_qcp.rhs = [tmp_qcp.rhs; sum_flux];
    tmp_qcp.sense = [tmp_qcp.sense; '<'];
    tmp_qcp.constrnames = [tmp_qcp.constrnames; 'UB_sum_flux'];
    
    % add matrix to QCP for minimization of absolute distances to random
    % flux vector
    delta_mat = zeros(n_rxns, size(tmp_qcp.A, 2) + 2*n_rxns);
    delta_mat(:, cellfun(@(x)find(ismember(tmp_qcp.varnames, x)), reactions)) = -speye(n_rxns);
    delta_mat(:, size(tmp_qcp.A, 2)+1:end) = [speye(n_rxns) -speye(n_rxns)];
    
    sample_qcp.A = [tmp_qcp.A sparse(size(tmp_qcp.A, 1), 2*n_rxns); sparse(delta_mat)];
    sample_qcp.lb = [tmp_qcp.lb; zeros(2*n_rxns, 1)];
    sample_qcp.ub = [tmp_qcp.ub; max_delta; max_delta];
    sample_qcp.vtype = [tmp_qcp.vtype; repmat('C', 2*n_rxns, 1)];
    sample_qcp.varnames = [tmp_qcp.varnames; strcat('delta_plus_', reactions);...
        strcat('delta_minus_', reactions)];
    sample_qcp.rhs = [tmp_qcp.rhs; zeros(n_rxns, 1)];
    sample_qcp.sense = [tmp_qcp.sense; repmat('=', n_rxns, 1)];
    sample_qcp.constrnames = [tmp_qcp.constrnames; strcat('delta_constr_', reactions)];
    
    % initialize weights for first norm minimization
    weights = zeros(n_rxns, 1);
    sample_qcp.obj = [zeros(size(tmp_qcp.A, 2), 1); weights; weights];
    
    sample_qcp.modelsense = 'min';
    
    % continue updating fields for quadratic constraints
    n_quad = numel(sample_qcp.quadcon);
    for j = 1:n_quad
        sample_qcp.quadcon(j).Qc = [tmp_qcp.quadcon(j).Qc zeros(size(tmp_qcp.quadcon(j).Qc, 1), 2*n_rxns);
            zeros(2*n_rxns, size(tmp_qcp.quadcon(j).Qc, 2)+2*n_rxns)];
        sample_qcp.quadcon(j).q = [tmp_qcp.quadcon(j).q; zeros(2*n_rxns, 1)];
    end
    
    % set minimum flux through biomass reaction
    sample_qcp.lb(ismember(sample_qcp.varnames, 'GROWTH')) = min_obj_pct*tmp_sol.mu_max;
    
    % generate as many random orderings as necessary to avoid using
    % reactions twice for minimization
    all_rxn_idx = 1:n_rxns;
    valid_rxn_idx = all_rxn_idx;
    n_select = max(1, ceil(p_select*numel(valid_rxn_idx)));
    rand_idx = cell2mat(arrayfun(@(i)...
        valid_rxn_idx(randperm(numel(valid_rxn_idx), numel(valid_rxn_idx))),...
        1:ceil(n_select*n_samples/numel(valid_rxn_idx)), 'un', 0));
    
    % solve minimization problems for sampling
    j = 1;
    tries = 0;
    tmp_samples = zeros(size(model.S, 2)+2, n_samples);
    while j <= n_samples && tries < max_tries
        tries = tries + 1;
        % create random flux vector
        rand_fluxes = (curr_max_fluxes-curr_min_fluxes) .* rand(n_rxns, 1) + curr_min_fluxes;
        
        weights = zeros(n_rxns, 1);
        curr_rand_idx = rand_idx((j-1)*n_select+1:j*n_select);
        weights(curr_rand_idx) = 1./curr_max_fluxes(curr_rand_idx);
        weights(isnan(weights)|isinf(weights)) = 0;
        sample_qcp.obj = [zeros(size(tmp_qcp.A, 2), 1); weights; weights];
        
        % add random fluxes to minimization problem (RHS)
        sample_qcp.rhs(startsWith(sample_qcp.constrnames, 'delta_constr_')) = -rand_fluxes;
        
        % solve QCP
        min_sol = gurobi(sample_qcp, solver_params);
        
        if isequal(min_sol.status, 'OPTIMAL')
            tmp_samples(:, j) = min_sol.x(~startsWith(sample_qcp.varnames, 'delta_'));
            j = j + 1;
        else
            n_failed(i) = n_failed(i) + 1;
        end
        
        if j > 1 && mod(j, 1000) == 1
            fprintf('Done with %d samples (%d%%)', j-1, round(100*(j-1)/n_samples))
            t2 = toc;
            fprintf(' (appr. %.f min remaining)\n',...
                (t2-t1)*(n_samples-j+1)/(j-1)/60)
        end
        
    end
    
    % calculate coverage
    coverages(i) = getSamplingCoverage(tmp_samples(rxn_idx, :),...
        curr_min_fluxes, curr_max_fluxes);
    
    samples{i} = sparse(tmp_samples);

    t3 = toc;
    fprintf('Time: %.f min (appr. %.f min remaining)\n\n', t3/60,...
        t3*(numel(temperatures)-i)/i/60)
    save(result_ws, 'min_fluxes', 'max_fluxes', 'samples', 'n_failed', ...
        'rgr', 'tmp_qcp', 'model', 'irradiance', 'min_obj_pct', 'temperatures', ...
        'coverages',...
        '-v7.3');
end
tend=toc;
fprintf('Total sampling time: %.2f min\n', tend/60)

function c = getSamplingCoverage(samples, minflux, maxflux)
%% c = getSamplingCoverage(samples, minflux, maxflux)
% Function to calculate the coverage of flux sampling (taken from GAPSPLIT)
% (Keaty & Jensen (2020), doi: 10.1093/bioinformatics/btz971).
%
% INPUT
% double samples:       matrix that contains the sampled fluxes
%                       (#variables x #samples)
% double minflux:       minima of the flux ranges, respectively
% double maxflux:       maxima of the flux ranges, respectively
%
% OUTPUT
% double c:             average coverage

[n_rxns, n_samples] = size(samples);

X = full([...
    reshape(minflux, 1, n_rxns);...
    reshape(maxflux, 1, n_rxns);...
    samples']);
n = (n_samples+2);
X = sort(X);
gaps = X(2:n,:) - X(1:n-1,:);
width = max(gaps,[],1);
rel = width ./ (X(n,:) - X(1,:));
c = 1 - mean(rel, 'omitnan');

end
