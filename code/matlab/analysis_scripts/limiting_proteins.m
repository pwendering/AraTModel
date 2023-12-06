% Find limiting proteins at different temperatures
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

% light intensitiy
irradiance = config('I'); % [umol m^-2 s^-1]

% temperature range
temperatures = celsius2kelvin(10:5:45);

% whether to adjust the total protein content to temperature
adj_ptot = true;

% protein IDs
protein_mets = model.mets(startsWith(model.mets, 'prot_') & ~contains(model.mets, 'pool'));
proteins = erase(protein_mets, 'prot_');
protein_names = cellfun(@(x)model.enzNames(ismember(model.enzymes, x)), proteins);

% solver parameters
solver_params = config('solver_params');

%% 2) Predict relative growth rates with thermostable proteins

% initialize array for relative growth rates
rgr = nan(numel(proteins), numel(temperatures));
rgr_opt = nan(numel(temperatures), 1);

for i = 1:numel(temperatures)
    
    [tmp_sol, ~, ~, qcp] = simulateTempEffects(model, irradiance,...
        'tempRange', temperatures(i),...
        'adjPtotFlag', adj_ptot);
    rgr_opt(i) = tmp_sol.mu_max;
    
    for j = 1:numel(proteins)
        
        tmp_qcp = qcp;
        
        % get original kcat entries
        row_idx = ismember(model.mets, ['prot_' proteins{j}]);
        kcats_orig = model.S(row_idx, model.S(row_idx, :) < 0);
        
        % replace adjusted kcats for current protein with original kcats
        row_idx = ismember(tmp_qcp.constrnames, ['prot_' proteins{j}]);
        
        if sum(tmp_qcp.A(row_idx, :) < 0) == numel(kcats_orig)
            tmp_qcp.A(row_idx, tmp_qcp.A(row_idx, :) < 0) = kcats_orig / config('kcat_scaling');
        else
            error('Numbers of kcats in original and adjusted model are not equal.')
        end
        
        % solve QCP
        solution = gurobi(tmp_qcp, solver_params);
        
        rgr(j, i) = solution.objval;
        
    end    
end

% get protein descriptions from UniProt
up_prot_info = getProtInfoUniProt(proteins, [], [], 'protein_name');
protein_descriptions = repmat({''}, size(proteins));
for i = 1:numel(proteins)
    if ~isfield(up_prot_info.to(i).proteinDescription, 'recommendedName')
        protein_descriptions{i} = up_prot_info.to(i).proteinDescription.submissionNames(1).fullName.value;
    else
        protein_descriptions{i} = up_prot_info.to(i).proteinDescription.recommendedName.fullName.value;
    end
end
clear up_prot_info

save(['engineering_targets_workspace_I_' num2str(irradiance)])

%% 3) Write results to file
filename = ['engineering_targets_I_' num2str(irradiance) '.xlsx'];
writetable(cell2table([proteins protein_names protein_descriptions num2cell(model.T_opt)],...
    'VariableNames', {'protein ID' 'protein name' 'protein description' 'T_opt'}),...
    filename,...
    'Range', 'A3');
comb_letters = char(combvec(char(65:65+25), char(65:65+25)));
letters = [cellstr(char(69:69+21)'); cellstr(comb_letters(2:-1:1, :)')];

for i = 1:numel(temperatures)
    
    % Calculate increase in relative growth rate
    rgr_inc = 100 * (rgr(:, i) ./ rgr_opt(i) - 1);
    % determine ranks for increases
    [~, sort_idx] = sort(rgr_inc, 'descend');
    [~, sort_idx] = sort(sort_idx, 'ascend');
    
    % write RGR and increases to file
    res_tab = table(rgr(:, i), rgr_inc, sort_idx,...
        'VariableNames', {'RGR (h^-1)', 'RGR increase (%)', 'rank'});
    n_col = size(res_tab, 2);
    writetable(res_tab, filename,...
        'Range', [letters{(i-1)*n_col+1} '3'])
    
    % temperature headers
    writetable(cell2table(cellstr(...
        [num2str(kelvin2celsius(temperatures(i))) ' °C'])),...
        filename,...
        'Range', [letters{1+(i-1)*n_col+1} '2'],...
        'WriteVariableNames', false)
end

%% 4) Clustering of growth responses to increased protein stability
rng('default')

min_rel_inc = 0.01;

% define temperatures that should be used for clustering
t_idx = 1:numel(temperatures);

% select protein that give an increased RGR at any of the temperatures
rel_icr = 100*(rgr(:, t_idx) ./ rgr_opt(t_idx)' - 1);
red_idx = any(rel_icr>min_rel_inc, 2);
rgr_red = rgr(red_idx, t_idx);
rel_icr_red = rel_icr(red_idx, :);
protein_names_red = protein_names(red_idx);
proteins_red = proteins(red_idx);

clust_data = rel_icr_red ./ max(rel_icr_red);
clust_data(isinf(clust_data)) = 0;

% K-medoids clustering

% find optimal K
K_array = 2:30;
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
    silh(i) = mean(silhouette(clust_data, clust,d_fun), 'omitnan');
end

% rerun k-medoids clustering with optimal K
K = K_array(silh == max(silh));
clust = kmedoids(clust_data, K,...
        'Options', kmedopt,...
        'Distance', d_fun,...
        'Algorithm', algorithm);

% plot clustering
colors = lines(K);
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
    
    min_rgr_pct = min(min(Y));
    max_rgr_pct = max(max(Y));
    ylim([min_rgr_pct max(1.4*max_rgr_pct, 2)])
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

saveas(gcf, ['engineering_targets_I_' num2str(irradiance) '.svg'])

% print proteins per cluster
for i = 1:K
    fprintf('#############\nCluster %i\n#############\n', i)
    arrayfun(@(j)disp([protein_names_red{j} ...
        ' (' char(regexprep(proteins_red{j}, '.*_(?<c>\w)$', '$<c>')) ')']),...
        find(clust==i))
end

% plot Silhouette index per K
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
saveas(gcf, ['silhouette_index_engineering_targets_kmedoids_I_' num2str(irradiance) '.png'])


%% 5) Venn diagram
% threshold for scaled relative increase in RGR
th_inc = 0.05;
% selected temperatures
t_dc = [10 25 40];
t_idx = find(ismember(kelvin2celsius(temperatures), t_dc));
set_names = arrayfun(@(i)sprintf('%d °C', i), t_dc, 'un', 0);

prots_per_t = cell(numel(t_idx), 1);
for i = 1:numel(t_idx)

    % increase in realtive growth rate
    rgr_icr = 100 * (rgr(:, t_idx(i)) / rgr_opt(t_idx(i)) - 1);
    rgr_icr_rel = rgr_icr / max(rgr_icr);
    
    % select those proteins with sufficient scaled and unscaled increase in
    % RGR
    sidx = rgr_icr_rel >= th_inc & rgr_icr > min_rel_inc;
    proteins_red = proteins(sidx);
    prot_names_red = protein_names(sidx);
    
    prots_per_t{i} = proteins_red;
end

if numel(set_names) == 3
    % three sets
sets = {
    setdiff(prots_per_t{1}, vertcat(prots_per_t{2:3}));... A
    setdiff(prots_per_t{2}, vertcat(prots_per_t{[1 3]}));... B
    setdiff(prots_per_t{3}, vertcat(prots_per_t{1:2}));... C
    setdiff(intersect(prots_per_t{1}, prots_per_t{2}), prots_per_t{3});... A&B
    setdiff(intersect(prots_per_t{1}, prots_per_t{3}), prots_per_t{2});... A&C
    setdiff(intersect(prots_per_t{2}, prots_per_t{3}), prots_per_t{1});... B&C
    intersect(intersect(prots_per_t{1}, prots_per_t{2}), prots_per_t{3})... A&B&C
    };
elseif numel(set_names) == 4
    % four sets
sets = {
    setdiff(prots_per_t{1}, vertcat(prots_per_t{2:4}));... A
    setdiff(prots_per_t{2}, vertcat(prots_per_t{[1 3 4]}));... B
    setdiff(prots_per_t{3}, vertcat(prots_per_t{[1 2 4]}));... C
    setdiff(prots_per_t{4}, vertcat(prots_per_t{1:3}));... D
    setdiff(intersect(prots_per_t{1}, prots_per_t{2}), union(prots_per_t{3}, prots_per_t{4}));... A&B
    setdiff(intersect(prots_per_t{1}, prots_per_t{3}), union(prots_per_t{2}, prots_per_t{4}));... A&C
    setdiff(intersect(prots_per_t{1}, prots_per_t{4}), union(prots_per_t{2}, prots_per_t{3}));... A&D
    setdiff(intersect(prots_per_t{2}, prots_per_t{3}), union(prots_per_t{1}, prots_per_t{4}));... B&C
    setdiff(intersect(prots_per_t{2}, prots_per_t{4}), union(prots_per_t{1}, prots_per_t{3}));... B&D
    setdiff(intersect(prots_per_t{3}, prots_per_t{4}), union(prots_per_t{1}, prots_per_t{2}));... C&D
    setdiff(intersect(intersect(prots_per_t{1}, prots_per_t{2}), prots_per_t{3}), prots_per_t{4});... A&B&C
    setdiff(intersect(intersect(prots_per_t{1}, prots_per_t{2}), prots_per_t{4}), prots_per_t{3});... A&B&D
    setdiff(intersect(intersect(prots_per_t{1}, prots_per_t{3}), prots_per_t{4}), prots_per_t{2});... A&C&D
    setdiff(intersect(intersect(prots_per_t{2}, prots_per_t{3}), prots_per_t{4}), prots_per_t{1});... B&C&D
    intersect(intersect(intersect(prots_per_t{1}, prots_per_t{2}), prots_per_t{3}), prots_per_t{4});... A&B&C&D
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
saveas(gcf, ['engineering_targets_venn_I_' num2str(irradiance) '_t_' num2str(th_inc) '.svg'])

%% 6) Barplot
% threshold for scaled relative increase in RGR
th_inc = 0.05;

% 400 umol m-2 s-1: thermal limitations only found at 10 °C and 40 °C
t_idx = find(ismember(kelvin2celsius(temperatures), [10 40]));
prots_per_t = cell(numel(t_idx), 1);
top_20_prots_per_t = cell(numel(t_idx), 1);
for i = 1:numel(t_idx)

    % increase in realtive growth rate
    rgr_icr = 100 * (rgr(:, t_idx(i)) / rgr_opt(t_idx(i)) - 1);
    rgr_icr_rel = rgr_icr / max(rgr_icr);
    sidx = rgr_icr_rel >= th_inc & rgr_icr > min_rel_inc;
    [rgr_icr_rel_sort, sort_idx] = sort(rgr_icr_rel(sidx), 'descend');
    prots_red = proteins(sidx);
    prots_red_sort = prots_red(sort_idx);
    
    prot_names_red = protein_names(sidx);
    prot_names_red_sort = prot_names_red(sort_idx);
    
    prot_descr_red = protein_descriptions(sidx);
    prot_descr_red_sort = prot_descr_red(sort_idx);
    
    prots_per_t{i} = prots_red_sort;
    top_20_prots_per_t{i} = prot_names_red_sort(1:min(numel(prots_red_sort),20));
end

all_prots = vertcat(top_20_prots_per_t{:});
prots_uniq = unique(all_prots);

plot_mat = zeros(numel(t_idx), numel(prots_uniq));
for i = 1:numel(t_idx)
    plot_mat(i,:) = cellfun(@(x)ismember(x, top_20_prots_per_t{i}), prots_uniq);
end

% count the number of times each of the top-limiting proteins is limiting
% across all tested temperatures
n_per_prot = cellfun(@(x)sum(ismember(all_prots, x)), prots_uniq);
[n_per_prot_sort, sort_idx] = sort(n_per_prot, 'descend');

% create barplot
bp = bar(plot_mat(:, sort_idx)', 'stacked',...
    'FaceColor', 'flat', 'EdgeColor', 'none');
colors = cool(numel(bp));
bar_colors =  mat2cell(colors, ones(numel(bp),1), 3); 
set(bp, {'CData'}, bar_colors)
xticks(1:numel(prots_uniq))
xticklabels(prots_uniq(sort_idx))
xtickangle(45)
yl = ylabel({'Number of temperatures', 'with thermal limitation'},...
    'horizontalalignment', 'left');
yl.Position(2) = 0;

legend(strsplit(num2str(kelvin2celsius(temperatures(t_idx)))),...
    'Box', 'off')

set(gca, 'FontName', 'Arial')
set(gcf, 'OuterPosition', [86.3333  301.0000  848.6667  398.0000])

saveas(gcf, ['engineering_targets_bar_I_' num2str(irradiance) '_t_' num2str(th_inc) '.svg'])

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

saveas(gcf, ['engineering_targets_heatmap_rho_s_I_' num2str(irradiance) '.svg'])

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
saveas(gcf, ['engineering_targets_heatmap_ji_I_' num2str(irradiance) '_t_' num2str(th_inc) '.svg'])
