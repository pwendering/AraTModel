%% Correlation between operational ranges and RGR
clear; clc
load sampling_workspace_I_150_20230817_constr_sum_of_flux.mat
op_ranges = max_fluxes - min_fluxes;
nrxns = size(op_ranges, 1);

% --- median operational ranges vs RGR ---
[r, p] = corr(median(op_ranges)', rgr', 'type', 'Pearson');
fprintf('Pearson correlation between median operational range and RGR: r = %.3g (P=%.3g)\n',...
    r, p)
[r, p] = corr(median(op_ranges)', rgr', 'type', 'Spearman');
fprintf('Spearman correlation between median operational range and RGR: r = %.3g (P=%.3g)\n',...
    r, p)

% --- operational ranges vs RGR ---
r = nan(nrxns, 1);
p = nan(nrxns, 1);
for i = 1:nrxns
    [r(i), p(i)] = corr(op_ranges(i, :)', rgr', 'type', 'Pearson');
end
fprintf('Average Pearson correlation between median operational range and RGR: r = %.3g\n',...
    mean(r, 'omitnan'))
fprintf('Median Pearson correlation between median operational range and RGR: r = %.3g\n',...
    median(r, 'omitnan'))

[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(p);
disp(sum(adj_p<0.05&r>=0.8)/numel(adj_p))

r = nan(nrxns, 1);
p = nan(nrxns, 1);
for i = 1:nrxns
    [r(i), p(i)] = corr(op_ranges(i, :)', rgr', 'type', 'Spearman');
end
fprintf('Average Spearman correlation between median operational range and RGR: r = %.3g\n',...
    mean(r, 'omitnan'))
fprintf('Median Spearman correlation between median operational range and RGR: r = %.3g\n',...
    median(r, 'omitnan'))
