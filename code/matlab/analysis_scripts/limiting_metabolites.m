% Identify lmiting metabolites at different temperatures
clear;clc
%% 1) Load model and specify conditions

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

% set ratio of 3:1 between NO3- and NH4+ uptake (optimal for Arabidopsis)
% https://doi.org/10.1007/s11104-010-0445-8
model = addMetabolite(model, 'nh3_nh4_ratio');
model.S(findMetIDs(model, 'nh3_nh4_ratio'),...
    findRxnIDs(model, {'Im_NO3', 'Im_NH4'})) = [1 -4];

% light intensitiy
irradiance = 400;%config('I'); % [umol m^-2 s^-1]

% temperature range
temperatures = celsius2kelvin(10:5:40);

% whether to adjust the total protein content to temperature
adj_ptot = true;

% metabolite IDs
metabolites = model.mets(~startsWith(model.mets, {'prot_', 'pmet_', 'nh3_nh4_ratio'}));
metabolite_names = model.metNames(findMetIDs(model, metabolites));

% solver parameters
solver_params = config('solver_params');

%% 2) Predict relative growth rates with supplied compounds

% initialize array for relative growth rates
rgr = nan(numel(metabolites), numel(temperatures));
rgr_opt = nan(numel(temperatures), 1);
import_fluxes = nan(sum(startsWith(model.rxns, 'Im_')), numel(temperatures));

for i = 1:numel(temperatures)
    
    [tmp_sol, ~, ~, qcp] = simulateTempEffects(model, irradiance,...
        'tempRange', temperatures(i),...
        'adjPtotFlag', adj_ptot);
    rgr_opt(i) = tmp_sol.mu_max;
    
    % store input fluxes
    import_rxn_idx = startsWith(qcp.varnames, {'Im_', 'TRANSPIRATION', 'PHOTON_UPTAKE'});
    import_fluxes(:, i) = tmp_sol.x_step2(import_rxn_idx);
    
    empty_column = zeros(size(qcp.rhs));
    
    for j = 1:numel(metabolites)
        
        tmp_qcp = qcp;
        
        % add exchange reaction for current metabolite
        tmp_col = empty_column;
        tmp_col(ismember(qcp.constrnames, metabolites(j))) = 1;
        
        tmp_qcp.A = [tmp_qcp.A tmp_col];
        tmp_qcp.lb = [tmp_qcp.lb; 0];
        tmp_qcp.ub = [tmp_qcp.ub; 1];
        tmp_qcp.obj = [tmp_qcp.obj; 0];
        tmp_qcp.quadcon.Qc = [tmp_qcp.quadcon.Qc zeros(size(tmp_qcp.quadcon.Qc, 1), 1);
            zeros(1, size(tmp_qcp.quadcon.Qc, 2) + 1)];
        tmp_qcp.quadcon.q = [tmp_qcp.quadcon.q; 0];
        tmp_qcp.varnames = [tmp_qcp.varnames; {['Im_' metabolites{j}]}];
        tmp_qcp.vtype = [tmp_qcp.vtype; 'C'];
                
        % solve QCP
        solution = gurobi(tmp_qcp, solver_params);
        
        rgr(j, i) = solution.objval;
        
    end
end

save(['supplementation_workspace_I_' num2str(irradiance)]);

%% 3) Clustering of growth responses to metabolite supplementation
rng('default')
min_rel_inc = 1;

% define temperatures that should be used for clustering
t_idx = 1:numel(temperatures);

% select metabolites that give an increased RGR at any of the temperatures
rel_icr = 100*(rgr(:, t_idx) ./ rgr_opt(t_idx)' - 1);
red_idx = any(rel_icr>min_rel_inc, 2);
rgr_red = rgr(red_idx, t_idx);
rel_icr_red = rel_icr(red_idx, :);
metabolite_names_red = metabolite_names(red_idx);
metabolites_red = metabolites(red_idx);

clust_data = rel_icr_red;
clust_data(isinf(clust_data)) = 0;

% K-medoids clustering

% find optimal K
K_array = 5:30;
algorithm = 'pam';
d_fun = 'cosine';
silh = nan(size(K_array));
kmedopt = statset;
kmedopt.Replicates = 10;
for i = 1:numel(K_array)
    % run clustering
    clust = kmedoids(clust_data, K_array(i),...
        'Options', kmedopt,...
        'Distance', d_fun,...
        'Algorithm', algorithm);
    % silhouette index
    silh(i) = median(silhouette(clust_data, clust, d_fun), 'omitnan');
end

% rerun k-medoids clustering with optimal K
K = K_array(silh == max(silh));
clust = kmedoids(clust_data, K,...
        'Options', kmedopt,...
        'Distance', d_fun,...
        'Algorithm', algorithm);
    
% plot clustering
colors = [[215,76,84]; [83,132,51]; [133,42,95]; [96,137,212]; [165,182,65];...
    [146,57,27]; [211,83,135]; [189,128,212]; [85,52,131]; [184,73,164];...
    [203,133,48]; [71,187,138]; [51,212,209]; [166,148,70]; [178,69,85];...
    [102,198,115]; [122,128,234]; [219,124,87]; [215,123,183]; [86,86,186]]/255;
colors = [colors; colors];
figure
t = tiledlayout('flow');
for i = 1:K
    nexttile
    hold on
    Y = clust_data(clust == i, :);
    plot(kelvin2celsius(temperatures(t_idx)),...
        Y,...
        'LineWidth', 1.3,...
        'Color', colors(i,:));
    
    min_y = min(min(Y));
    max_y = max(max(Y));
    ylim([min_y max(1.4*max_y, max_y+5)])
    xlim(kelvin2celsius([min(temperatures(t_idx)) max(temperatures(t_idx))]))
    
    y_limits = get(gca, 'YLim');
    
    x_limits = get(gca, 'XLim');
    rect_y = y_limits(2)-0.2*range(y_limits);
    rectangle(...
        'Position', [x_limits(1) rect_y range(x_limits) y_limits(2)-rect_y],...
        'FaceColor', [.8 .8 .8],...
        'EdgeColor', 'k',...
        'LineWidth', 1.3)
    text(0.5, 0.9, sprintf('Cluster %d (n=%d)', i, sum(clust == i)),...
        'Units', 'normalized',...
        'FontName', 'Arial',...
        'FontSize', 8,...
        'FontWeight', 'bold',...
        'HorizontalAlignment', 'center')

    set(gca,...
        'Box', 'on',...
        'FontName', 'Arial',...
        'LineWidth', 1.3)
end
t.XLabel.String = 'Temperature (°C)';
t.XLabel.FontSize = 14;

t.YLabel.String = 'Scaled growth increment';
t.YLabel.FontSize = 14;

t.Title.String = strrep(char(d_fun), '_', '\_');

set(gcf, 'OuterPosition', 1000*[0.1977    0.1203    1.0513    0.5787])

saveas(gcf, ['supplementation_kmedoids_I_' num2str(irradiance) '.svg'])

% print metabolites per cluster
for i = 1:K
    fprintf('#############\nCluster %i\n#############\n', i)
    arrayfun(@(j)disp([metabolite_names_red{j} ...
        ' (' char(regexprep(metabolites_red{j}, '.*_(?<c>\w)$', '$<c>')) ')']),...
        find(clust==i))
end

% plot Silhouette index at different K
figure
plot(K_array, silh, 'k', 'linewidth', 1.5)
hold on
mc = lines(2);
scatter(K, silh(K_array==K), 'filled', 'MarkerFaceColor', mc(2, :),...
        'MarkerEdgeColor', mc(2, :))
xlabel('K', 'fontsize', 14)
ylabel('Median Silhouette Index', 'fontsize', 14)
h = scatter(NaN, NaN, 'filled', 'MarkerFaceColor', mc(2, :),...
    'MarkerEdgeColor', mc(2, :));
legend(h, {'selected'},...
    'box', 'off',...
    'fontname', 'Arial',...
    'Fontsize', 14,...
    'location', 'se')
set(gca,...
    'Box', 'off',...
    'FontName', 'Arial')
saveas(gcf, ['silhouette_index_supplementation_kmedoids_I_' num2str(irradiance) '.png'])

%% 4) Write results to file
clust_all_mets = nan(size(metabolites));
clust_all_mets(red_idx) = clust;
filename = ['supplementation_I_' num2str(irradiance) '.xlsx'];
writetable(table(metabolites, metabolite_names, clust_all_mets,...
    'VariableNames', {'metabolite ID' 'metabolite name' 'cluster label'}),...
    filename,...
    'Range', 'A3');
letters = char(68:68+23);
for i = 1:numel(temperatures)
    
    % Calculate increase in relative growth rate
    rgr_inc = 100 * (rgr(1:end-1, i) ./ rgr_opt(i) - 1);
        % determine ranks for increases
    [~, sort_idx] = sort(rgr_inc, 'descend');
    [~, sort_idx] = sort(sort_idx, 'ascend');
    
    % write RGR and increases to file
    res_tab = table(rgr(1:end-1, i), rgr_inc, sort_idx,...
        'VariableNames', {'RGR (h^-1)', 'RGR increase (%)', 'rank'});
    n_col = size(res_tab, 2);
    writetable(res_tab, filename,...
        'Range', [letters((i-1)*n_col+1) '3'])
    
    % temperature headers
    writetable(cell2table(cellstr(...
        [num2str(kelvin2celsius(temperatures(i))) ' °C'])),...
        filename,...
        'Range', [letters(1+(i-1)*n_col+1) '2'],...
        'WriteVariableNames', false)
end

%% 5) Venn diagram
% threshold for scaled relative increase in RGR
th_inc = 0.1;
% selected temperatures
t_dc = [10 25 40];
t_idx = find(ismember(kelvin2celsius(temperatures), t_dc));
set_names = arrayfun(@(i)sprintf('%d °C', i), t_dc, 'un', 0);

mets_per_t = cell(numel(t_idx), 1);
for i = 1:numel(t_idx)

    % increase in realtive growth rate
    rgr_icr = 100 * (rgr(:, t_idx(i)) / rgr_opt(t_idx(i)) - 1);
    rgr_icr_rel = rgr_icr / max(rgr_icr);
    sidx = rgr_icr_rel >= th_inc & rgr_icr > min_rel_inc;
    metabolites_red = metabolites(sidx);
    met_names_red = metabolite_names(sidx);
    
    mets_per_t{i} = metabolites_red;
end

if numel(set_names) == 3
    % three sets
sets = {
    setdiff(mets_per_t{1}, vertcat(mets_per_t{2:3}));... A
    setdiff(mets_per_t{2}, vertcat(mets_per_t{[1 3]}));... B
    setdiff(mets_per_t{3}, vertcat(mets_per_t{1:2}));... C
    setdiff(intersect(mets_per_t{1}, mets_per_t{2}), mets_per_t{3});... A&B
    setdiff(intersect(mets_per_t{1}, mets_per_t{3}), mets_per_t{2});... A&C
    setdiff(intersect(mets_per_t{2}, mets_per_t{3}), mets_per_t{1});... B&C
    intersect(intersect(mets_per_t{1}, mets_per_t{2}), mets_per_t{3})... A&B&C
    };
elseif numel(set_names) == 4
    % four sets
sets = {
    setdiff(mets_per_t{1}, vertcat(mets_per_t{2:4}));... A
    setdiff(mets_per_t{2}, vertcat(mets_per_t{[1 3 4]}));... B
    setdiff(mets_per_t{3}, vertcat(mets_per_t{[1 2 4]}));... C
    setdiff(mets_per_t{4}, vertcat(mets_per_t{1:3}));... D
    setdiff(intersect(mets_per_t{1}, mets_per_t{2}), union(mets_per_t{3}, mets_per_t{4}));... A&B
    setdiff(intersect(mets_per_t{1}, mets_per_t{3}), union(mets_per_t{2}, mets_per_t{4}));... A&C
    setdiff(intersect(mets_per_t{1}, mets_per_t{4}), union(mets_per_t{2}, mets_per_t{3}));... A&D
    setdiff(intersect(mets_per_t{2}, mets_per_t{3}), union(mets_per_t{1}, mets_per_t{4}));... B&C
    setdiff(intersect(mets_per_t{2}, mets_per_t{4}), union(mets_per_t{1}, mets_per_t{3}));... B&D
    setdiff(intersect(mets_per_t{3}, mets_per_t{4}), union(mets_per_t{1}, mets_per_t{2}));... C&D
    setdiff(intersect(intersect(mets_per_t{1}, mets_per_t{2}), mets_per_t{3}), mets_per_t{4});... A&B&C
    setdiff(intersect(intersect(mets_per_t{1}, mets_per_t{2}), mets_per_t{4}), mets_per_t{3});... A&B&D
    setdiff(intersect(intersect(mets_per_t{1}, mets_per_t{3}), mets_per_t{4}), mets_per_t{2});... A&C&D
    setdiff(intersect(intersect(mets_per_t{2}, mets_per_t{3}), mets_per_t{4}), mets_per_t{1});... B&C&D
    intersect(intersect(intersect(mets_per_t{1}, mets_per_t{2}), mets_per_t{3}), mets_per_t{4});... A&B&C&D
    };
end
venn_labels = cellfun(@numel, sets);

venn(numel(set_names),...
    'sets', set_names,...
    'labels', venn_labels,...
    'colors', zeros(numel(set_names), 3),...
    'alpha', 0,...
    'edgeC', [0 0 0],...
    'edgeW', 1);
set(gcf, 'OuterPosition', [649.0000  443.6667  260.6667  265.3333])
saveas(gcf, ['supplementation_venn_I_' num2str(irradiance) '_t_' num2str(th_inc) '.svg'])

%% 6) Barplot
% threshold for scaled relative increase in RGR
th_inc = 0.1;
min_rel_inc = 1;
t_idx = find(ismember(kelvin2celsius(temperatures), kelvin2celsius(temperatures)));
biomass_mets = ismember(metabolites, [findMetsFromRxns(model,...
    {'Bio_opt', 'protein', 'DNA', 'RNA', 'carbohydrate', 'lipid'});...
    {'starch1_h', 'starch3_h', 'starch5_h', 'cellulose1_c', 'RBC_inact'}']);
mets_per_t = cell(numel(t_idx), 1);
incr_per_t = cell(numel(t_idx), 1);
top_x_mets_per_t = cell(numel(t_idx), 1);
top_x_incr_per_t = cell(numel(t_idx), 1);
for i = 1:numel(t_idx)

    % increase in relative growth rate
    rgr_icr = 100 * (rgr(:, t_idx(i)) / rgr_opt(t_idx(i)) - 1);
    rgr_icr_rel = rgr_icr / max(rgr_icr);
    % select those with sufficient scaled in unscaled increase in RGR
    % exclude biomass reaction metabolites
    sidx = rgr_icr_rel >= th_inc & rgr_icr > min_rel_inc & ~biomass_mets;
    
    rgr_icr_rel_red = rgr_icr_rel(sidx);
    rgr_icr_red = rgr_icr(sidx);
    [~, sort_idx] = sort(rgr_icr_red, 'descend');
    
    rgr_icr_red_sort = rgr_icr_red(sort_idx);
    rgr_icr_rel_red_sort = rgr_icr_rel_red(sort_idx);
    
    metabolites_red = metabolites(sidx);
    metabolites_red_sort = metabolites_red(sort_idx);
    
    met_names_red = metabolite_names(sidx);
    met_names_red_sort = met_names_red(sort_idx);
    
    mets_per_t{i} = metabolites_red_sort;
    incr_per_t{i} = rgr_icr_rel_red_sort;
    top_x_mets_per_t{i} = metabolites_red_sort(1:10);
    top_x_incr_per_t{i} = rgr_icr_rel_red_sort(1:10);
end

all_mets = vertcat(top_x_mets_per_t{:});
mets_uniq = unique(all_mets);

plot_mat = zeros(numel(t_idx), numel(mets_uniq));
for i = 1:numel(t_idx)
    plot_mat(i,:) = cellfun(@(x)ismember(x, top_x_mets_per_t{i}), mets_uniq);
end

% count the number of times each metabolite is among the top-limiting
% metabolites at each temperature
n_per_met = cellfun(@(x)sum(ismember(all_mets, x)), mets_uniq);
[n_per_met_sort, sort_idx] = sort(n_per_met, 'ascend');

% create barplot
bp = bar(plot_mat(:, sort_idx)', 'stacked',...
    'FaceColor', 'flat',...
    'EdgeColor', 'none',...
    'Horizontal', 'on');
colors = [
    0.6118    0.7804    0.8510
    0.6046    0.6791    0.7497
    0.5974    0.5778    0.6484
    0.5902    0.4765    0.5471
    0.5830    0.3752    0.4458
    0.5758    0.2739    0.3444
    0.5686    0.1725    0.2431];
bar_colors =  mat2cell(colors(t_idx,:), ones(numel(bp),1), 3); 
set(bp, {'CData'}, bar_colors)
yticks(1:numel(mets_uniq))
y_lab_met_names = metabolite_names(cellfun(@(x)find(ismember(metabolites, x)),...
    mets_uniq(sort_idx)));
y_lab_comps = regexprep(mets_uniq(sort_idx), '.*_(?<c>\w)', ' ($<c>)');
yticklabels(strcat(y_lab_met_names, y_lab_comps))
xlabel('limitation occurrence',...
    'horizontalalignment', 'center')
xlim([0 numel(t_idx)])
legend(strcat(strsplit(num2str(kelvin2celsius(temperatures(t_idx)))), ' °C'),...
    'Box', 'off')

set(gca,...
    'FontName', 'Arial',...
    'FontSize', 10,...
    'TickLength', [0 0],...
    'box', 'off')
set(gcf, 'OuterPosition', [350.3333  261.6667  354.6667  437.3333])

saveas(gcf, ['supplementation_bar_I_' num2str(irradiance) '_t_' num2str(th_inc) '.svg'])

%% 7) Similarity between growth response vectors
rgr_icr = 100 * (rgr ./ rgr_opt' - 1);
lab = arrayfun(@(i)sprintf('%d °C', i), kelvin2celsius(temperatures), 'un', 0);

% Spearman correlation
[rho_s, p_corr_s] = corr(rgr_icr, 'type', 'Spearman');
figure
imagesc(rho_s)
xticks(1:size(rho_s, 1))
xticklabels(lab)
xtickangle(90)
yticks(1:size(rho_s, 1))
yticklabels(lab)
colormap('summer')
cb = colorbar;
cb.Title.String = '\rho_S';
set(gcf, 'OuterPosition', [11.6667  369.6667  324.6667  340.6667])

saveas(gcf, ['supplementation_heatmap_rho_s_I_' num2str(irradiance) '.svg'])

% Jaccard Index
rgr_icr_rel = rgr_icr ./ max(rgr_icr);
rgr_icr_rel_th = rgr_icr_rel >= th_inc & rgr_icr > 1e-6;
n = size(rgr_icr_rel_th, 2);
n_mets = size(rgr_icr_rel_th, 1);
ji = ones(n);
for i = 1:n-1
    for j = i+1:n
        ji(i, j) = sum(rgr_icr_rel_th(:, i)&rgr_icr_rel_th(:, j))/...
            sum(rgr_icr_rel_th(:, i)|rgr_icr_rel_th(:, j));
        ji(j, i) = ji(i, j);
    end
end

figure
imagesc(ji)
xticks(1:size(ji, 1))
xticklabels(lab)
xtickangle(90)
yticks(1:size(ji, 1))
yticklabels(lab)
colormap('summer')
cb = colorbar;
cb.Title.String = sprintf('JI (t=%.2g)', th_inc);

set(gcf, 'OuterPosition', [328.3333  369.6667  324.6667  340.6667])
saveas(gcf, ['supplementation_heatmap_ji_I_' num2str(irradiance) '_t_' num2str(th_inc) '.svg'])