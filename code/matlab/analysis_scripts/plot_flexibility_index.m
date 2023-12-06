% Calculate the flexibility index per reaction as the ratio between the
% interquatile range and the median flux across all sampled flux
% distributions

%% Read input workspace
clear; clc

% load sampling workspace
load sampling_workspace_I_150_20230817_constr_sum_of_flux.mat

model.subSystems(cellfun(@isempty, model.subSystems)) = {''};
n_asm_idx = cellfun(@(x)contains(x,...
    {'nitrate assimilation', 'glutamate synthesis', 'glutamine synthesis'}), model.subSystems);
model.subSystems(n_asm_idx) = strcat(model.subSystems(n_asm_idx), ', N assimilation');

% reaction IDs
reactions = model.rxns(~startsWith(model.rxns, {'draw_', 'prot_'}));
subsystems = model.subSystems(~startsWith(model.rxns, {'draw_', 'prot_'}));
reactions(ismember(reactions, 'Bio_opt')) = {'GROWTH'};
reactions(ismember(reactions, 'Im_hnu')) = {'PHOTON_UPTAKE'};
reactions(ismember(reactions, 'Im_H2O_REV')) = {'TRANSPIRATION'};
n_rxns = numel(reactions);

% Determine statistics of sampled reaction fluxes
iqr_fluxes = zeros(numel(reactions), numel(temperatures));
median_fluxes = zeros(numel(reactions), numel(temperatures));
for i = 1:numel(reactions)
    
    rxn_idx = ismember(tmp_qcp.varnames, reactions(i));
    flux_samples = cell2mat(arrayfun(@(tp)samples{tp}(rxn_idx, :),...
        1:numel(temperatures), 'un', 0)')';
    
    q25 = quantile(flux_samples, 0.25, 1);
    q75 = quantile(flux_samples, 0.75, 1);
    iqrange = q75-q25;
    
    iqr_fluxes(i, :) = iqrange;
    median_fluxes(i, :) = median(flux_samples, 1);
end
clear q25 q75 iqrange flux_samples rxn_idx

save('metabolic_flexibility_clustering_rxn_fluxes', 'iqr_fluxes', 'median_fluxes');

%% Plot differences in IQR/Sum(v) across all reactions

% remove ratios that are either Inf or NaN or all-zero
ratios = iqr_fluxes./median_fluxes;
keep_idx = ~any(isnan(ratios)|isinf(ratios), 2) & ~all(ratios<1e-6, 2);
rxns_red = reactions(keep_idx);
subs_red = subsystems(keep_idx);
ratios = ratios(keep_idx, :);

% remove ratios that do not have a range of at least 10% variation from the
% mean
keep_idx = any(abs(1-ratios./mean(ratios, 2)) >= 0.1, 2);
rxns_red = rxns_red(keep_idx);
subs_red = subs_red(keep_idx);
ratios = ratios(keep_idx, :);

% cluster temperature responses by K-medoids clustering
rng('default')

clust_data = ratios;
for i = 1:size(clust_data, 1)
    clust_data(i,:) = zscore(smooth(clust_data(i,:), 3));
end

kmedopt = statset;
kmedopt.Replicates = 1000;
algorithm = 'pam';

K=5;
% rerun k-medoids clustering with chosen K and distance metric
[clust, MedoidsBest, sumDbest, Dbest, Midx, info] = kmedoids(clust_data, K,...
        'Options', kmedopt,...
        'Distance', 'sqEuclidean',...
        'Algorithm', algorithm);

figure
hold on
colors = [[106, 140, 205]; [98, 182, 80]; [138, 99, 202];...
    [185, 177, 66]; [205, 86, 168]; [77, 181, 152]; [204, 77, 62];...
    [101, 125, 57]; [195, 99, 127]; [196, 130, 65]]/255;

for i = 1:K
    plot(MedoidsBest(i, :), 'Color', colors(i, :), 'LineWidth', 2)
end
xticks(1:numel(temperatures))
xticklabels(round(kelvin2celsius(temperatures)))
xlabel('Temperature (°C)')
ylabel({'scaled flexibility index', 'of reaction fluxes'})
xlim([1 numel(temperatures)])

set(gca, 'FontName', 'Arial', 'FontSize', 10)
set(gcf, 'OuterPosition', [891.6667  339.6667  361.3333  356.6667],...
    'Renderer', 'painters')

saveas(gcf, 'flexibility_index_rxn_fluxes.svg')

%% print clusters

for i = 1:K
    fprintf('Cluster %i\n', i)
    r = rxns_red(clust==i);
    s = subs_red(clust==i);
    for j = 1:numel(r)
        fprintf('\t%s\t%s\n', r{j}, char(s{j}))
    end
end

%% t-SNE
rng default
close all
subs_red = cellfun(@(x)char(x), subs_red, 'un', 0);

[Y, loss] = tsne(ratios, 'Distance', 'correlation');

subs_red_mod = cellfun(@(x)strtok(x, ','), subs_red, 'un', 0);
subs_red_mod_red = subs_red_mod;
subs_red_mod_red(cellfun(@isempty, subs_red_mod_red)) = {'NA'};

% find optimal K
kmedopt = statset;
kmedopt.Replicates = 1000;

% for each distance measure, find K that minimizes the maximum correlation
% between median flexibility index profiles

d_fun = {'sqEuclidean', 'Euclidean', 'seuclidean', 'cityblock', 'minkowski',...
    'chebychev'};

min_max_corr = ones(size(d_fun));
k_opt = nan(size(d_fun));
min_opt = optimset;
min_opt.MaxFunEvals = 1000;
for i = 1:numel(d_fun)
    fun = @(k)kmedoids_eval(Y, clust_data, k, d_fun{i}, kmedopt);
    for j = 1:10
        [tmp_k, tmp_c] = fminbnd(fun, 2, 15);
        if tmp_c < min_max_corr(i)
            k_opt(i) = tmp_k;
            min_max_corr(i) = tmp_c;
        end
    end
end

best_k = k_opt(min_max_corr==min(min_max_corr));
best_max_corr = min_max_corr(min_max_corr==min(min_max_corr));
best_dist = d_fun(min_max_corr==min(min_max_corr));

disp(best_k)
disp(best_max_corr)
disp(best_dist)

K = round(best_k(1));
[clust, medoids] = kmedoids(Y, K, 'Options', kmedopt,...
    'Distance', best_dist{1});

% color by cluster
figure
hold on

colors = [[106, 140, 205]; [98, 182, 80]; [138, 99, 202];...
    [185, 177, 66]; [205, 86, 168]; [77, 181, 152]; [204, 77, 62];...
    [101, 125, 57]; [195, 99, 127]; [196, 130, 65]]/255;

for i = 1:K
    
    % find most-abundant pathway per cluster
    tmp_subs = cellfun(@(x)strtrim(strsplit(x, ',')), subs_red(clust==i), 'un', 0);
    tmp_subs = [tmp_subs{:}]';
    tmp_subs_uniq = setdiff(tmp_subs, {'transport', 'import', 'export'});
    n_per_sub = cellfun(@(x)sum(ismember(tmp_subs, x)), tmp_subs_uniq);
    [n_per_sub_sort, ia] = sort(n_per_sub, 'descend');
    tmp_subs_uniq_sort = tmp_subs_uniq(ia);
    
    % top three most represented subsystems
    top_3_subs = tmp_subs_uniq_sort(1:3);
    
    tmp_idx = find(clust==i);
    for j = 1:numel(tmp_idx)
        if contains(subs_red(tmp_idx(j)), top_3_subs(1))
            alpha = 1;
            f = -0.2;
        elseif contains(subs_red(tmp_idx(j)), top_3_subs(2))
            alpha = 1;
            f = 0;
        elseif contains(subs_red(tmp_idx(j)), top_3_subs(3))
            alpha = 1;
            f = 0.2;
        else
            alpha = 0;
            f = 0;
        end
        scatter(Y(tmp_idx(j), 1), Y(tmp_idx(j), 2), 70, min(max(colors(i, :)+f, 0), 1),...
            'o', 'filled',...
            'MarkerFaceAlpha', alpha,...
            'MarkerEdgeColor', colors(i, :),...
            'LineWidth', 1)
    end
    
    
    % plot most abundant subsystem names
    y_pos =  medoids(i, 2) + [0.5 0 -0.5];
    f = [-0.2 0 0.2];
    for j = 1:numel(top_3_subs)
        text(medoids(i, 1), y_pos(j), top_3_subs{j},...
            'FontSize', 8,...
            'FontName', 'Arial',...
            'HorizontalAlignment', 'left',...
            'Color', min(max(colors(i, :)+f(j), 0), 1))
    end
end
xlabel('tSNE 1')
ylabel('tSNE 2')
set(gca, 'FontName', 'Arial')
saveas(gcf, 'flux_flexibility_tnse_by_cluster.svg')

% color by subsystem
figure
hold on
markers = {'v', 's', '^', 'd', 'o'};

subs_red_mod_red_uniq = unique(subs_red_mod_red, 'stable');
colors = colorcube(numel(subs_red_mod_red_uniq)+1);
colors(end, :) = [];

for i = 1:K
    tmp_colors = colors(cellfun(@(x)find(ismember(subs_red_mod_red_uniq,x)),...
        subs_red_mod_red(clust==i)), :);
    scatter(Y(clust==i, 1), Y(clust==i, 2), 70,...
        markers{i}, 'filled',...
        'CData', tmp_colors)
end

h = zeros(numel(subs_red_mod_red_uniq), 1);
for i = 1:numel(subs_red_mod_red_uniq)
    h(i) = scatter(NaN, NaN, 100, colors(i, :), 's', 'filled');
end
legend(h, subs_red_mod_red_uniq, 'box', 'off')
saveas(gcf, 'flux_flexibility_tnse_by_subsystem.svg')
        
% plot subsystems names
figure
hold on
scatter(Y(:, 1), Y(:, 2), 'MarkerFaceColor', 'none')
for i = 1:size(Y, 1)
    text(Y(i, 1), Y(i, 2), subs_red_mod_red{i},...
        'FontSize', 3,...
        'FontName', 'Arial',...
        'HorizontalAlignment', 'center');
end
saveas(gcf, 'flux_flexibility_tnse_subs_text.svg')

% plot cluster medians
figure
hold on
colors = [[106, 140, 205]; [98, 182, 80]; [138, 99, 202];...
    [185, 177, 66]; [205, 86, 168]; [77, 181, 152]; [204, 77, 62];...
    [101, 125, 57]; [195, 99, 127]; [196, 130, 65]]/255;

for i = 1:K
    md = median(clust_data(clust==i, :));
    sd = std(clust_data(clust==i, :));
    min_err = md - sd;
    max_err = md + sd;
    x = kelvin2celsius([temperatures(1); temperatures(1); temperatures(2:end)';...
        temperatures(end:-1:2)']);
    y = [max_err(1); min_err(1); min_err(2:end)'; max_err(end:-1:2)'];
    F = fill(x, y, colors(i, :), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    plot(kelvin2celsius(temperatures), md,...
        'Color', colors(i, :),...
        'LineWidth', 3)
end

h = zeros(K, 1);
for i = 1:K
    h(i) = plot(NaN, NaN,...
        'Color', colors(i, :),...
        'LineWidth', 3);
end
legend(h, strsplit(num2str(1:K)), 'box', 'off')
xlabel('Temperature (°C)')
ylabel('scaled flexibility index of reaction flux')
set(gca, 'LineWidth', 1.3,...
    'FontSize', 10,...
    'FontName', 'Arial')
set(gcf, 'OuterPosition', [965.0000  379.0000  302.6667  322.6667]) 
saveas(gcf, 'flux_flexibility_tnse_medoids.svg')

    