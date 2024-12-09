% Test different ratios of chloroplastic O2 and CO2

%% load temperature model
tmp = load(config('tgemFile'));
model = tmp.TGEM;

%% Predict RGR and Anet temperature responses with different gamma values

% temperature range
T = celsius2kelvin(linspace(10, 40, 10));
% light intensity
I = 800;

% === gamma must be changed in the config file! ===

% gamma = Oa / Ca
simResults_1 = simulateTempEffects(model, I, 'tempRange', T);
% gamma = Oa / Ca / 0.75
simResults_2 = simulateTempEffects(model, I, 'tempRange', T);
% gamma = Oa / Ca / 0.75^2
simResults_3 = simulateTempEffects(model, I, 'tempRange', T);
% gamma = Oa / Ca / 0.75^3
simResults_4 = simulateTempEffects(model, I, 'tempRange', T);

save(sprintf('test_gamma_workspace_I%i', I), 'model', 'T', 'I', ...
    'simResults_1', 'simResults_2', 'simResults_3', 'simResults_4')

%% Plot results

% RGR
fig1 = figure;
hold on
plot(simResults_1.mu, 'Color', 'k')
plot(simResults_2.mu, 'Color', [0 0.4470 0.7410])
plot(simResults_3.mu, 'Color', [0.8500 0.3250 0.0980])
plot(simResults_4.mu, 'Color', [0.9290 0.6940 0.1250])
xlabel('Temperature (°C)')
ylabel('RGR (h^{-1})')
xticklabels(round(kelvin2celsius(T)))
xlim([1 numel(T)])
legend({'\gamma_1 = O_a/C_a', '\gamma_2 = 1.33*O_a/C_a',...
    '\gamma_3 = 1.78*O_a/C_a', '\gamma_4 = 2.37*O_a/C_a'}, ...
    'Box', 'off', 'Location', 'best')
text(0.1, 0.8, 'O_c = \gamma * C_c', 'Units', 'normalized', 'FontSize', 14)
set(gca, 'FontSize', 14)
% saveas(fig1, fullfile(outdir, sprintf('gamma_test_rgr_I%i.pdf', I)), 'pdf')
saveas(fig1, fullfile(outdir, sprintf('gamma_test_rgr_I%i.svg', I)), 'svg')

% Anet
fig2 = figure;
hold on
plot(simResults_1.A_max, 'Color', 'k')
plot(simResults_2.A_max, 'Color', [0 0.4470 0.7410])
plot(simResults_3.A_max, 'Color', [0.8500 0.3250 0.0980])
plot(simResults_4.A_max, 'Color', [0.9290 0.6940 0.1250])
xlabel('Temperature (°C)')
ylabel('A_{net} (mmol m^{-2} s^{-1})')
xticklabels(round(kelvin2celsius(T)))
xlim([1 numel(T)])
set(gca, 'FontSize', 14)
legend({'\gamma_1 = O_a/C_a', '\gamma_2 = 1.33*O_a/C_a',...
    '\gamma_3 = 1.78*O_a/C_a', '\gamma_4 = 2.37*O_a/C_a'}, ...
    'Box', 'off', 'Location', 'best')
text(0.3, 0.2, 'O_c = \gamma * C_c', 'Units', 'normalized', 'FontSize', 14)
set(gca, 'FontSize', 14)
% saveas(fig2, fullfile(outdir, sprintf('gamma_test_anet_I%i.pdf', I)), 'pdf')
saveas(fig2, fullfile(outdir, sprintf('gamma_test_anet_I%i.svg', I)), 'svg')

% Fluxes
histogram((simResults_1.x_step2./simResults_1.mu) - (simResults_2.x_step2./simResults_2.mu))

% congruence coefficient
v1 = reshape(simResults_1.x_step2, numel(simResults_1.x_step2), []);
v2 = reshape(simResults_2.x_step2, numel(simResults_2.x_step2), []);
v3 = reshape(simResults_3.x_step2, numel(simResults_3.x_step2), []);
cos_v1v2 = v1'*v2 / sqrt((v1'*v1)*(v2'*v2));
cos_v1v3 = v1'*v3 / sqrt((v1'*v1)*(v3'*v3));

% Rv coefficient
rvfun = @(X,Y)trace(X*X'*Y*Y')/sqrt(trace(X*X'*X*X')*trace(Y*Y'*Y*Y'));
rv_12 = rvfun(simResults_1.x_step2, simResults_2.x_step2);
rv_13 = rvfun(simResults_1.x_step2, simResults_3.x_step2);

% simulation 1 vs. simulation 2
fig3 = figure;
corr_mat = corr(simResults_1.x_step2, simResults_2.x_step2, 'Type', 'Spearman');
imagesc(corr_mat)
xticklabels(round(kelvin2celsius(T)))
yticklabels(round(kelvin2celsius(T)))
xlabel('Temperature (°C)')
ylabel('Temperature (°C)')
cb = colorbar;
cb.Title.String = 'r_S';
colormap('summer')
yyaxis('right')
yticklabels([])
ylabel('\gamma = O_a/C_a', 'Color', 'k')
title('\gamma = 1.33*O_a/C_a', 'FontWeight', 'normal')
set(gca, 'FontSize', 14)
set(fig3, 'OuterPosition', [13.6667  295.0000  488.0000  420.0000])
saveas(fig3, fullfile(outdir, 'heatmap_gamma_flux_step2_v1_v2.pdf'), 'pdf')

% simulation 1 vs. simulation 2
fig4 = figure;
corr_mat = corr(simResults_1.x_step2, simResults_3.x_step2, 'Type', 'Spearman');
imagesc(corr_mat)
xticklabels(round(kelvin2celsius(T)))
yticklabels(round(kelvin2celsius(T)))
xlabel('Temperature (°C)')
ylabel('Temperature (°C)')
cb = colorbar;
cb.Title.String = 'r_S';
colormap('summer')
yyaxis('right')
yticklabels([])
ylabel('\gamma = O_a/C_a', 'Color', 'k')
title('\gamma = 1.78*O_a/C_a', 'FontWeight', 'normal')
set(gca, 'FontSize', 14)
set(fig4, 'OuterPosition', [489.0000  295.0000  488.0000  420.0000])
saveas(fig4, fullfile(outdir, 'heatmap_gamma_flux_step2_v1_v3.pdf'), 'pdf')

% boxplots of correlations of temperature response patterns across
% reactions
corrtype = 'Spearman';
corr_flux_over_temp_g1_g2 = diag(corr(simResults_1.x_step2', simResults_2.x_step2', 'Type', corrtype));
corr_flux_over_temp_g1_g3 = diag(corr(simResults_1.x_step2', simResults_3.x_step2', 'Type', corrtype));
corr_flux_over_temp_g1_g4 = diag(corr(simResults_1.x_step2', simResults_4.x_step2', 'Type', corrtype));

fig5 = figure;
hold on

corr_flux_over_temp_g1_g2_sort = sort(corr_flux_over_temp_g1_g2);
n_corr_sort_g1_g2 = sum(~isnan(corr_flux_over_temp_g1_g2_sort));
corr_flux_over_temp_g1_g3_sort = sort(corr_flux_over_temp_g1_g3);
n_corr_sort_g1_g3 = sum(~isnan(corr_flux_over_temp_g1_g3_sort));
corr_flux_over_temp_g1_g4_sort = sort(corr_flux_over_temp_g1_g4);
n_corr_sort_g1_g4 = sum(~isnan(corr_flux_over_temp_g1_g4_sort));
n_rxns_max = max([n_corr_sort_g1_g2 n_corr_sort_g1_g3 n_corr_sort_g1_g4]);


s1 = scatter((1:numel(corr_flux_over_temp_g1_g2)), corr_flux_over_temp_g1_g2_sort, 25);
s2 = scatter((1:numel(corr_flux_over_temp_g1_g3)), corr_flux_over_temp_g1_g3_sort, 25);
s3 = scatter((1:numel(corr_flux_over_temp_g1_g4)), corr_flux_over_temp_g1_g4_sort, 25);

q10_g1_g2 = quantile(corr_flux_over_temp_g1_g2, 0.1);
q10_g1_g3 = quantile(corr_flux_over_temp_g1_g3, 0.1);
q10_g1_g4 = quantile(corr_flux_over_temp_g1_g4, 0.1);

x_lim = get(gca, 'XLim');

line([x_lim(1) n_rxns_max], [q10_g1_g2 q10_g1_g2], 'Color', s1.CData)
line([x_lim(1) n_rxns_max], [q10_g1_g3 q10_g1_g3], 'Color', s2.CData)
line([x_lim(1) n_rxns_max], [q10_g1_g4 q10_g1_g4], 'Color', s3.CData)
xlim(x_lim)
ylim([-1 1.1])

xlabel('Reactions')
ylabel(sprintf('%s correlation', corrtype))
legend({...
    sprintf('r_S(\\gamma_1, \\gamma_2), n=%i', n_corr_sort_g1_g2), ...
    sprintf('r_S(\\gamma_1, \\gamma_3), n=%i', n_corr_sort_g1_g3), ...
    sprintf('r_S(\\gamma_1, \\gamma_4), n=%i', n_corr_sort_g1_g4) ...
    }, ...
    'Box', 'off', 'Location', 'best')
set(gca, 'FontSize', 14)
set(fig5, 'OuterPosition', [353.0000  225.6667  375.3333  372.0000])

saveas(fig5, 'flux_corr_gamma.svg', 'svg')