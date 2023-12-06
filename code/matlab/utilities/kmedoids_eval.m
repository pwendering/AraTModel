function max_corr = kmedoids_eval(X, raw, K, dist_measure, options)
K = round(K);
clust = kmedoids(X, K,...
    'Options', options,...
    'Distance', dist_measure);
if numel(unique(clust)) ~= K
    fprintf('Number of required clusters not reached\n')
    max_corr = 1;
else
    median_val = cell2mat(arrayfun(@(i)median(raw(clust==i,:), 1), 1:K, 'un', 0)');
    corr_mat = corr(median_val');
    corr_mat(diag(true(K, 1))) = nan;
    max_corr = max(max(corr_mat));
end

end