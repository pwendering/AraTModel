% Fit beta growth function to TPP data for all proteins in the Meltome
% Atlas

meltome_files = dir(fullfile('protein-stability', 'Meltome_Atlas'));

proteinIDs = [];
funParams = [];
rmse = [];
r2adj = [];
T_m = [];
T_m_meltome = [];
T_opt = [];
fN_TH = [];
T_H = [];

org_names = [];
for i=1:numel(meltome_files)
    if contains(meltome_files(i).name, 'protein-stability')
        tmp_file = fullfile(meltome_files(i).folder, meltome_files(i).name);
        tmp_tab = readtable(tmp_file);
        [~,ia] = unique(tmp_tab.Protein_ID,'stable');
        tmp_tm_meltome = tmp_tab.meltPoint(ia);
        try
            [tmp_ids,~, tmp_fit,tmp_rmse,tmp_tm,tmp_topt,tmp_th,tmp_fnth,tmp_r2adj] = fitTPPBetaFcn(tmp_file, [], false);
            proteinIDs = [proteinIDs; tmp_ids];
            funParams = [funParams; tmp_fit];
            rmse = [rmse; tmp_rmse];
            r2adj = [r2adj; tmp_r2adj];
            T_m = [T_m; tmp_tm];
            T_opt = [T_opt; tmp_topt];
            fN_TH = [fN_TH; tmp_fnth];
            T_H = [T_H; tmp_th];
            org_names = [org_names; ...
                repmat(cellstr(regexprep(tmp_file, '.*\\(?<org>([^\\\-]+)).*\.csv', '$<org>')),...
                numel(tmp_ids), 1)];
            T_m_meltome = [T_m_meltome; tmp_tm_meltome];
            clear tmp_tab tmp_tm_meltome ia
        catch ME
            fprintf('Fitting not successful for %s\n', tmp_file)
            disp(ME.message)
            clear tmp_tab tmp_tm_meltome ia
        end
        
        
    end
end

save('meltome_ml_ws')

%remove datasets from intact cells
rem_idx = contains(org_names, 'cells');
org_names(rem_idx) = [];
proteinIDs(rem_idx) = [];
funParams(rem_idx,:) = [];
rmse(rem_idx) = [];
r2adj(rem_idx) = [];
T_m(rem_idx) = [];
fN_TH(rem_idx) = [];
T_H(rem_idx) = [];
T_opt(rem_idx) = [];

%% write results as table
writetable(...
    cell2table(...
    [org_names proteinIDs num2cell([funParams T_m T_opt rmse r2adj fN_TH])],...
    'variablenames', {'dataset', 'protein_id', 'beta_param_1',...
    'beta_param_2', 'beta_param_3', 'melt_point', 't_opt', 'rmse', 'r2_adj', 'fn_th'}),...
    fullfile('protein-stability', 'Meltome_Atlas', 'curve_params.csv'));

writetable(...
    cell2table([proteinIDs num2cell(T_H)], 'VariableNames', {'protein_id' 'T_H'}),...
    fullfile('protein-stability', 'Meltome_Atlas', 't_h.csv'));

%% plot key temperatures with fit quality

% adjusted R2
n_bins = 10;
[bins, edges] = discretize(r2adj(~isnan(T_m)), n_bins);

figure
colormap(parula(n_bins))
scatter(T_m(~isnan(T_m)), T_opt(~isnan(T_m)), 'filled',...
    'CData', bins)
cb = colorbar('Ticks', 1:.9:numel(edges),...
    'TickLabels', arrayfun(@(i)sprintf('%.2f', i), edges, 'un', 0));
cb.Title.String = 'R^2_{adj}';

xlabel('T_m [°C]')
ylabel('T_{opt} [°C]')

set(gca,...
    'Box', 'on',...
    'LineWidth', 1,...
    'FontSize', 12,...
    'Position', [0.1764    0.2033    0.56    0.65])

set(gcf, 'OuterPosition', [477.0000  203.0000  477.3333  436.6667])

exportgraphics(gcf, 'R2_key_temp_meltome_all.png')

% RMSE
n_bins = 10;
[bins, edges] = discretize(rmse(~isnan(T_m)), n_bins);

figure
colormap(parula(n_bins))
scatter(T_m(~isnan(T_m)), T_opt(~isnan(T_m)), 'filled',...
    'CData', bins, 'MarkerFaceAlpha', .1)
cb = colorbar('Ticks', 1:.9:numel(edges),...
    'TickLabels', arrayfun(@(i)sprintf('%.2f', i), edges, 'un', 0));
cb.Title.String = 'RMSE';

xlabel('T_m [°C]')
ylabel('T_{opt} [°C]')

set(gca,...
    'Box', 'on',...
    'LineWidth', 1,...
    'FontSize', 14,...
    'Position', [0.1764    0.2033    0.56    0.65])

set(gcf, 'OuterPosition', [477.0000  203.0000  477.3333  436.6667])

exportgraphics(gcf, 'RMSE_key_temp_meltome_all.png')

