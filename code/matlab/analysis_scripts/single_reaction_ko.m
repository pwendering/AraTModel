% Prediction and validation of temperature-induced lethality and reduced growth
clear;clc
%% 1) Data preparation

% read model
model_file = config('tgemFile');
tmp = load(model_file);
model = tmp.TGEM;
clear tmp

% add RubisCO activase module
model = addRCA(model);

% block uptake of H2S and chloroplastic alanine transaminase
model.ub(findRxnIDs(model, {'Im_H2S', 'Im_H2S_REV'})) = 0;
model.ub(findRxnIDs(model, {'AlaTA_h', 'AlaTA_h_REV'})) = 0;
model.ub(findRxnIDs(model, {'Bio_AA', 'Bio_CLim', 'Bio_NLim'})) = 0;

% correct THF transport reaction
model.S(findMetIDs(model, 'H_h'), findRxnIDs(model, {'Tr_THF_h', 'Tr_THF_h_REV'})) = [-1 1];

% light intensitiy
irradiance = config('I'); % [umol m^-2 s^-1]
fprintf('Irradiance: %.4g umol m^-2 s^-1\n', irradiance)

% reaction IDs
reactions = model.rxns(~startsWith(model.rxns, 'prot_') & ~contains(model.rxns, 'draw_prot_'));
reactions(ismember(reactions, 'Bio_opt')) = {'GROWTH'};
reactions(ismember(reactions, 'Im_hnu')) = {'PHOTON_UPTAKE'};
reactions(ismember(reactions, 'Im_H2O_REV')) = {'TRANSPIRATION'};

% solver parameters
solver_params = config('solver_params');
% solver_params = rmfield(solver_params, {'PoolSearchMode', 'PoolSolutions'});

% simulation temperatures
temperatures = [17 25 27 35 45];

c = clock;
datestr = sprintf('%d%02.0f%02.0f', c(1:3));
result_ws = ['thermosensitivity_reaction_level_I_' num2str(irradiance) '_' datestr];

%% 2) Single knock-outs on protein level

% initialize result arrays
n_temp = numel(temperatures);
n_rxns = numel(reactions);

rgr_opt = zeros(n_temp, 1);
co2_exc_opt = zeros(n_temp, 1);
rgr_kos = nan(n_rxns, n_temp);
co2_exc_kos = nan(n_rxns, n_temp);

red_rgr_ids = cell(n_temp, 1);
red_rgr_50_ids = cell(n_temp, 1);
lethal_ids = cell(n_temp, 1);

qcps = cell(n_temp, 1);

for i = 1:n_temp
    
    % determine optimal solution at current temperature
    [tmp_sol_wt, ~, ~, tmp_qcp_wt] = simulateTempEffects(model, irradiance,...
        'tempRange', celsius2kelvin(temperatures(i)));
    qcps{i} = tmp_qcp_wt;
    rgr_opt(i) = tmp_sol_wt.mu_max;
    co2_exc_idx = ismember(tmp_qcp_wt.varnames, model.C_ID);
    co2_exc_opt(i) = tmp_sol_wt.x_step2(co2_exc_idx);
    
    % find reactions that carry flux in the optimal solution
    rxn_idx = ismember(tmp_qcp_wt.varnames, reactions);
    nz_rxn_idx = tmp_sol_wt.x_step1(rxn_idx) > 0;
    
    % loop over model protein and simulate KOs
    for j = 1:numel(reactions)
        
        if nz_rxn_idx(j)
            % block draw reaction for current protein
            tmp_qcp_ko = tmp_qcp_wt;
            
            tmp_qcp_ko.lb(ismember(tmp_qcp_ko.varnames, reactions(j))) = 0;
            tmp_qcp_ko.ub(ismember(tmp_qcp_ko.varnames, reactions(j))) = 0;
            
            % optimize flux through biomass reaction
            tmp_sol_ko = gurobi(tmp_qcp_ko, solver_params);
            
            if ~isequal(tmp_sol_ko.status, 'OPTIMAL')
                fprintf('Protein %s: %s\n', reactions{j}, tmp_sol_ko.status)
            else
                if isfield(tmp_sol_ko, 'pool')
                    rgr_kos(j, i) = max(vertcat(tmp_sol_ko.pool.objval));
                else
                    rgr_kos(j, i) = tmp_sol_ko.objval;
                end
            end
            
            % minimize CO2 uptake at the optimal objective value
            tmp_qcp_ko.lb(find(tmp_qcp_ko.obj)) = tmp_sol_ko.objval - solver_params.FeasibilityTol;
            tmp_qcp_ko.obj(:) = 0;
            tmp_qcp_ko.obj(co2_exc_idx) = 1;
            tmp_qcp_ko.modelsense = 'min';
            tmp_sol_ko = gurobi(tmp_qcp_ko, solver_params);
            co2_exc_kos(j, i) = tmp_sol_ko.x(co2_exc_idx);
        end
    end
    
    red_rgr_ids{i} = reactions(rgr_kos(:, i)>solver_params.FeasibilityTol & ...
        rgr_kos(:, i)<0.99*rgr_opt(i));
    red_rgr_50_ids{i} = reactions(rgr_kos(:, i)>solver_params.FeasibilityTol & ...
        rgr_kos(:, i)<0.5*rgr_opt(i));
    lethal_ids{i} = reactions(rgr_kos(:, i)<=solver_params.FeasibilityTol);
    
end

save(result_ws)

%% 3) Print results
red_rgr_names = red_rgr_ids;
lethal_names = lethal_ids;
for i = 1:n_temp

    n_r = numel(red_rgr_ids{i});
    for j = 1:n_r
        red_rgr_names{i}{j} = model.rxnNames{ismember(qcps{i}.varnames, red_rgr_ids{i}{j})};
    end

    n_r = numel(lethal_ids{i});
    for j = 1:n_r
        lethal_names{i}{j} = model.rxnNames{ismember(qcps{i}.varnames, lethal_ids{i}{j})};
    end

end

% Reduction at 17 °C, but not at 27 °C
result_file = ['thermosensitivity_reaction_level_I_' num2str(irradiance) '_' datestr '.txt'];
fid = fopen(result_file, 'w');
fprintf(fid, 'Reduction only 17 °C\n');
fprintf(fid, 'Gene\tProtein\tReaction ID\tReaction name\tRGR reduction\n');
for i = 1:numel(red_rgr_names{temperatures==17})
    rxn = red_rgr_ids{temperatures==17}(i);
    if ~ismember(rxn, red_rgr_ids{temperatures==27})
        if startsWith(rxn, 'arm_')
            mets = findMetsFromRxns(model, rxn);
            pmet = mets(startsWith(mets, 'pmet_'));
            rxns = findRxnsFromMets(model, pmet);
            mets = findMetsFromRxns(model, rxns);
        else
            mets = findMetsFromRxns(model, rxn);
        end
        
        mets = mets(startsWith(mets, 'prot_'));
        prots = erase(mets, 'prot_');
        
        up_info = getProtInfoUniProt(prots, [], [], 'xref_araport');
        araport_ids = repmat({'no match'}, size(up_info, 1), 1);
        for j = 1:size(up_info, 1)
            if isfield(up_info.to(j).uniProtKBCrossReferences, 'id')
                araport_ids{j} = up_info.to(j).uniProtKBCrossReferences.id;
            end
        end
        rgr_red = 100-100*rgr_kos(ismember(reactions, rxn),...
            temperatures==17)/rgr_opt(temperatures==17);
        for j = 1:numel(prots)
            fprintf(fid, '%s\t%s\t', araport_ids{j}, prots{j});
            fprintf(fid, '%s\t%s\t', char(rxn), red_rgr_names{temperatures==17}{i});
            fprintf(fid, '%.2f%%\n', rgr_red);
        end
    end
end

fprintf(fid, 'Lethal only 17 °C\n');
fprintf(fid, 'Reaction\tProtein(s)\tGene(s)\n');
for i = 1:numel(red_rgr_names{temperatures==17})
    rxn = lethal_ids{temperatures==17}(i);
    if ~ismember(rxn, lethal_ids{temperatures==27})
        if startsWith(rxn, 'arm_')
            mets = findMetsFromRxns(model, rxn);
            pmet = mets(startsWith(mets, 'pmet_'));
            rxns = findRxnsFromMets(model, pmet);
            mets = findMetsFromRxns(model, rxns);
        else
            mets = findMetsFromRxns(model, rxn);
        end
        mets = mets(startsWith(mets, 'prot_'));
        prots = erase(mets, 'prot_');
        
        up_info = getProtInfoUniProt(prots, [], [], 'xref_araport');
        araport_ids = repmat({'no match'}, size(up_info, 1), 1);
        for j = 1:size(up_info, 1)
            if isfield(up_info.to(j).uniProtKBCrossReferences, 'id')
                araport_ids{j} = up_info.to(j).uniProtKBCrossReferences.id;
            end
        end
        
        fprintf(fid, '%s\t', [lethal_names{temperatures==17}{i} '[' char(rxn) ']']);
        fprintf(fid, '%s\t', strjoin(prots, ', '));
        fprintf(fid, '%s\n', strjoin(araport_ids, ', '));
    end
end

% Reduction at 27 °C, but not at 17 °C
fprintf(fid, 'Reduction only 27 °C\n');
fprintf(fid, 'Reaction\tRGR reduction\tProtein(s)\tGene(s)\n');
for i = 1:numel(red_rgr_names{temperatures==27})
    rxn = red_rgr_ids{temperatures==27}(i);
    if ~ismember(rxn, red_rgr_ids{temperatures==17})
        if startsWith(rxn, 'arm_')
            mets = findMetsFromRxns(model, rxn);
            pmet = mets(startsWith(mets, 'pmet_'));
            rxns = findRxnsFromMets(model, pmet);
            mets = findMetsFromRxns(model, rxns);
        else
            mets = findMetsFromRxns(model, rxn);
        end
        mets = mets(startsWith(mets, 'prot_'));
        prots = erase(mets, 'prot_');
        
        up_info = getProtInfoUniProt(prots, [], [], 'xref_araport');
        araport_ids = repmat({'no match'}, size(up_info, 1), 1);
        for j = 1:size(up_info, 1)
            if isfield(up_info.to(j).uniProtKBCrossReferences, 'id')
                araport_ids{j} = up_info.to(j).uniProtKBCrossReferences.id;
            end
        end
        
        rgr_red = 100-100*rgr_kos(ismember(reactions, rxn),...
            temperatures==27)/rgr_opt(temperatures==27);
        for j = 1:numel(prots)
            fprintf(fid, '%s\t%s\t', araport_ids{j}, prots{j});
            fprintf(fid, '%s\t%s\t', char(rxn), red_rgr_names{temperatures==27}{i});
            fprintf(fid, '%.2f%%\n', rgr_red);
        end
    end
end

fprintf(fid, 'Lethal only 27 °C\n');
fprintf(fid, 'Reaction\tProtein(s)\tGene(s)\n');
for i = 1:numel(lethal_ids{temperatures==27})
    rxn = lethal_ids{temperatures==27}(i);
    if ~ismember(rxn, lethal_ids{temperatures==17})
        if startsWith(rxn, 'arm_')
            mets = findMetsFromRxns(model, rxn);
            pmet = mets(startsWith(mets, 'pmet_'));
            rxns = findRxnsFromMets(model, pmet);
            mets = findMetsFromRxns(model, rxns);
        else
            mets = findMetsFromRxns(model, rxn);
        end
        mets = mets(startsWith(mets, 'prot_'));
        prots = erase(mets, 'prot_');
        
        up_info = getProtInfoUniProt(prots, [], [], 'xref_araport');
        araport_ids = repmat({'no match'}, size(up_info, 1), 1);
        for j = 1:size(up_info, 1)
            if isfield(up_info.to(j).uniProtKBCrossReferences, 'id')
                araport_ids{j} = up_info.to(j).uniProtKBCrossReferences.id;
            end
        end

        fprintf(fid, '%s\t', [lethal_names{temperatures==27}{i} '[' char(rxn) ']']);
        fprintf(fid, '%s\t', strjoin(prots, ', '));
        fprintf(fid, '%s\n', strjoin(araport_ids, ', '));
    end
end
fclose(fid);

%% 4) Find negative controls
rgr_kos_for_ctrl = rgr_kos;
for i = 1:size(rgr_kos_for_ctrl, 1)
    if any(~isnan(rgr_kos_for_ctrl(i, :)))
        rgr_kos_for_ctrl(i, isnan(rgr_kos_for_ctrl(i, :))) = 1;
    end
end
neg_ctrl_idx = all(rgr_kos_for_ctrl./rgr_opt' > 0.99, 2);

rxn_ids = reactions(neg_ctrl_idx);
rxn_names = model.rxnNames(findRxnIDs(model, rxn_ids));

fid = fopen(['neg_ctrl_reaction_level_I_' num2str(irradiance) '_' datestr '.txt'], 'w');
for i = 1:numel(rxn_ids)
    if startsWith(rxn_ids(i), 'arm_')
        mets = findMetsFromRxns(model, rxn_ids(i));
        pmet = mets(startsWith(mets, 'pmet_'));
        rxns = findRxnsFromMets(model, pmet);
        mets = findMetsFromRxns(model, rxns);
    else
        mets = findMetsFromRxns(model, rxn_ids(i));
    end
    mets = mets(startsWith(mets, 'prot_'));
    prots = erase(mets, 'prot_');
    
    up_info = getProtInfoUniProt(prots, [], [], 'xref_araport');
    araport_ids = repmat({'no match'}, size(up_info, 1), 1);
    for j = 1:size(up_info, 1)
        if isfield(up_info.to(j).uniProtKBCrossReferences, 'id')
            araport_ids{j} = up_info.to(j).uniProtKBCrossReferences.id;
        end
    end
    
    for j = 1:numel(prots)
        fprintf(fid, '%s\t%s\t', araport_ids{j}, prots{j});
        fprintf(fid, '%s\t%s\n', rxn_ids{i}, rxn_names{i});
    end
end
fclose(fid);

%% 5) RGR per exchanged CO2
rgr_opt_norm = rgr_opt ./ co2_exc_opt;
rgr_kos_norm = rgr_kos ./ co2_exc_kos;
rgr_ratios_norm = rgr_kos_norm ./ rgr_opt_norm';

red_rgr_norm_ids = cell(n_temp, 1);
red_rgr_norm_50_ids = cell(n_temp, 1);
lethal_norm_ids = cell(n_temp, 1);

for i = 1:n_temp
    red_rgr_norm_ids{i} = reactions(rgr_kos_norm(:, i)>solver_params.FeasibilityTol & ...
        rgr_ratios_norm(:, i)<0.99);
    red_rgr_norm_50_ids{i} = reactions(rgr_kos_norm(:, i)>solver_params.FeasibilityTol & ...
        rgr_ratios_norm(:, i)<0.5);
    lethal_norm_ids{i} = reactions(rgr_kos_norm(:, i)<=solver_params.FeasibilityTol);
end

% Reduction only at 27 °C
t_idx = temperatures==27;
r_only_27 = setdiff(union(red_rgr_norm_ids{t_idx}, red_rgr_ids{t_idx}),...
    union(red_rgr_norm_ids{temperatures==17}, red_rgr_ids{temperatures==17}));
red_idx = cellfun(@(x)find(ismember(reactions, x)), r_only_27);
r_names = model.rxnNames(red_idx);
r_ids = reactions(red_idx);
res_table = table(...
    arrayfun(@(i)[r_names{i} ' [' r_ids{i} ']'], 1:numel(red_idx), 'un', 0)',...
    100-100*rgr_kos(red_idx, t_idx)./rgr_opt(t_idx),...
    100-100*rgr_ratios_norm(red_idx, t_idx),...
    'VariableNames', {'Reaction', 'RGR Reduction', 'RGR/CO2 Reduction'});
res_table_sort = sortrows(res_table, 'RGR Reduction', 'descend');
disp(res_table_sort)

% Reduction only at 17 °C
t_idx = temperatures==17;
r_only_17 = setdiff(union(red_rgr_norm_ids{t_idx}, red_rgr_ids{t_idx}),...
    union(red_rgr_norm_ids{temperatures==27}, red_rgr_ids{temperatures==27}));
red_idx = cellfun(@(x)find(ismember(reactions, x)), r_only_17);
r_names = model.rxnNames(red_idx);
r_ids = reactions(red_idx);
res_table = table(...
    arrayfun(@(i)[r_names{i} ' [' r_ids{i} ']'], 1:numel(red_idx), 'un', 0)',...
    r_ids,...
    100-100*rgr_kos(red_idx, t_idx)./rgr_opt(t_idx),...
    100-100*rgr_ratios_norm(red_idx, t_idx),...
    'VariableNames', {'Reaction', 'ID', 'RGR Reduction', 'RGR/CO2 Reduction'});
res_table_sort = sortrows(res_table, 'RGR Reduction', 'descend');
disp(res_table_sort)

