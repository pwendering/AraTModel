% Figure 2 A-C
% Comparison of predicted and experimental relative growth rates
% Comparison of predicted and experimental net CO2 assimilation rate with
%   the dataset from Weston et al. 2011
% Boxplots of feasible flux ranges at different temperatures
clear; clc; close all

fig2 = figure;
pnxbr = 0.06;
pnybr = 0.06;
pnxpad = 0.07;
pnypad = 0.2;
pnxwd_3 = (1-2*pnxbr-2*pnxpad)/3;
pnywd_3 = (1-2*pnybr-pnxpad-pnypad)/3;
pnxwd_4 = (1-3*pnxbr-2*pnxpad)/4;
fs = 10;
bx_lw = 1.3;

p_idx = 0;
letters = char(65:65+25);
letter_pos = [-0.23 1.05];
letter_fs = 12;

% temperature range
temperatures = celsius2kelvin(linspace(10, 40, 50));

% load colormap
tmp = load('blue_red_colormap.mat');
colors = tmp.blueredcm;
min_temp = 6;
max_temp = 42;
n_colors = max_temp-min_temp+1;
clear tmp

%% Panels A and B (relative growth rate and net CO2 assimilation rate)

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

% -------------------- Panel A
panel_a = axes(fig2, 'Units', 'normalized',...
    'Position', [pnxpad pnypad+2*(pnywd_3+pnybr) pnxwd_3 pnywd_3]);

% scatter plot experimental vs predicted RGR
t_exp = [12 16 24 20 21 20 20 20 20 20 16 6 16 22 28 22 22 22];
t_exp_uniq = unique(t_exp);
rgr_exp = [0.0063 0.0083 0.0100 0.0082 0.0071 0.0092 0.0077 0.0077 0.0104 ...
     0.0076 0.0055 0.0020 0.0124 0.0124 0.0132 0.0141 0.0119 0.0077];
 
% light intensities
i_array = 150;

rgr_pred = zeros(numel(rgr_exp), numel(i_array));
for i = 1:numel(i_array)
    tmp_sol = simulateTempEffects(model, i_array(i),...
        'tempRange', celsius2kelvin(t_exp_uniq));
    
    for j = 1:numel(t_exp_uniq)
        exp_idx = ismember(t_exp, t_exp_uniq(j));
        rgr_pred(exp_idx, i) = repelem(tmp_sol.mu_max(j), sum(exp_idx));
    end
end

[~,~,bins] = histcounts(t_exp, n_colors, 'BinLimits', [min_temp max_temp]);

% create scatter plot of experimental vs predicted RGR
mrk = 'os^*';
pearson_r = zeros(size(i_array));
corr_pval = zeros(size(i_array));
ptcolors = colors(bins, :);
for i = 1:numel(i_array)
    lm = fitlm(rgr_exp, rgr_pred(:, i));
    line(panel_a, [min(rgr_exp) max(rgr_exp)],...
        [lm.predict(min(rgr_exp)) lm.predict(max(rgr_exp))],...
        'Color', [.1 .1 .1],...
        'LineWidth', 1.5)
    hold on
    scatter(panel_a, rgr_exp, rgr_pred(:, i),...
        60,...
        'filled',...
        'Marker', mrk(i),...
        'CData', ptcolors,...
        'MarkerEdgeColor', [.4 .4 .4])
    [pearson_r(i), corr_pval(i)] = corr(rgr_exp', rgr_pred(:, i));
    
end

text(panel_a, 0.2, 0.5, sprintf('r = %.2f', pearson_r(1)), 'Units', 'normalized',...
    'FontSize', fs-2, 'FontName', 'Arial')

h = zeros(size(i_array));
for i = 1:numel(h)
    h(i) = scatter(panel_a, NaN, NaN,...
        100,...
        'filled',...
        'Marker', mrk(i),...
        'MarkerFaceColor', [.7 .7 .7],...
        'MarkerEdgeColor', [.4 .4 .4]);
end

ax = gca;
ax.XAxis.Exponent = round(mean(log10(rgr_exp)));

legend(h,...
    arrayfun(@(i)sprintf('%i \\mumol m^{-2} s^{-1}', i), i_array, 'un', 0),...
    'box', 'off',...
    'fontname', 'Arial',...
    'fontsize', fs-2,...
    'position', [0.15 0.91 0.001 0.02])
xlabel(panel_a, 'RGR measured (h^{-1})');
ylabel(panel_a, 'RGR predicted (h^{-1})');

p_idx = p_idx + 1;
text(panel_a, letter_pos(1), letter_pos(2), letters(p_idx),...
    'FontSize', letter_fs,...
    'FontName', 'Arial',...
    'FontWeight', 'bold',...
    'Units', 'normalized')

set(panel_a,...
    'box', 'off',...
    'LineWidth', bx_lw,...
    'FontSize', fs,...
    'FontName', 'Arial',...
    'XColor', 'k',...
    'YColor', 'k')

% -------------------- Panel B
panel_b = axes(fig2, 'Units', 'normalized',...
    'Position', [pnxpad+pnxwd_3+pnxbr pnypad+2*(pnywd_3+pnybr) 0.8*pnxwd_3 pnywd_3]);

% load A_net predictions from workspace
% (generated using script 'correlation_a_net_predictions.m' in
% code/matlab/analysis_scripts directory)
load a_net_pred_meas_workspace.mat a_net_tab no_imp_ca_idx a_cbm

filter_idx = no_imp_ca_idx;

A_exp = a_net_tab.A_mumol_m2_s(filter_idx);
A_cbm = a_cbm(filter_idx);
T_exp = a_net_tab.Temperature_C(filter_idx);

% colors
[~,~,bins] = histcounts(T_exp, n_colors, 'BinLimits', [min_temp max_temp]);
ptcolors = colors(bins, :);

% line indicating perfect correlation
lm = fitlm([0 1], [0 1]);
line(panel_b, [min(A_exp) max(A_exp)],...
        [lm.predict(min(A_exp)) lm.predict(max(A_exp))],...
        'Color', [.1 .1 .1],...
        'LineWidth', 1.5)
hold on
% scatter plot measured vs. predicted
scatter(panel_b, A_exp, A_cbm, 60,...
    'filled',...
    'CData',  ptcolors,...
    'MarkerEdgeColor', [.4 .4 .4])

xlabel('A measured (\mumol m^{-2} s^{-1})')
ylabel('A predicted (\mumol m^{-2} s^{-1})')

pearson_r = corr(A_exp, A_cbm, 'type', 'Pearson');
text(panel_b, 0.8, 0.2, sprintf('r = %.2f', pearson_r), 'Units', 'normalized',...
    'FontSize', fs-2, 'FontName', 'Arial')

p_idx = p_idx + 1;
text(letter_pos(1), letter_pos(2), letters(p_idx),...
    'FontSize', letter_fs,...
    'FontName', 'Arial',...
    'FontWeight', 'bold',...
    'Units', 'normalized')

set(gca,...
    'box', 'off',...
    'LineWidth', bx_lw,...
    'FontSize', fs,...
    'FontName', 'Arial',...
    'XColor', 'k',...
    'YColor', 'k')

%% Panel C (feasible flux ranges)
panel_c = axes(fig2, 'Units', 'normalized',...
    'Position', [pnxpad+2*(pnxwd_3+pnxbr) pnypad+2*(pnywd_3+pnybr) 0.8*pnxwd_3 pnywd_3]);

load flux_ranges_I_150_20230829.mat

fs_ranges_log10 = fs_ranges;
fs_ranges_log10(fs_ranges_log10<=0) = NaN;
fs_ranges_log10 = log10(fs_ranges_log10);

y_min = floor(min(min(fs_ranges_log10)));
y_max = ceil(max(max(fs_ranges_log10)));

% feasible ranges
[~,~,bins] = histcounts(kelvin2celsius(temperatures), n_colors, 'BinLimits', [min_temp max_temp]);
bxcolor = colors(bins, :);
for i = 1:size(fs_ranges_log10, 2)
    bc = boxchart(panel_c, categorical(repelem(i, size(fs_ranges_log10, 1), 1)), fs_ranges_log10(:, i),...
        'MarkerStyle', '.',...
        'JitterOutliers', 'on',...
        'MarkerColor', 'k',...
        'BoxFaceAlpha', .5);
    hold on
    bc.BoxFaceColor = bxcolor(i, :);
end
xticklabels(panel_c, round(kelvin2celsius(temperatures)))
xtickangle(panel_c, 90)
xlabel(panel_c, 'Temperature (°C)')

ylabel(panel_c,{'log_{10} feasible flux ranges', '(mmol gDW^{-1} h^{-1})'})

p_idx = p_idx + 1;
text(panel_c, letter_pos(1), letter_pos(2), letters(p_idx),...
    'FontSize', letter_fs,...
    'FontName', 'Arial',...
    'FontWeight', 'bold',...
    'Units', 'normalized')

set(panel_c,...
    'Box', 'off',...
    'LineWidth', bx_lw,...
    'FontName', 'Arial',...
    'FontSize', fs,...
    'YLim', [y_min y_max],...
    'XColor', 'k',...
    'YColor', 'k')

%% overall colorbar
ncbtcks = 5;
cbtck = linspace(0, ncbtcks, ncbtcks)/ncbtcks;
cbtcklab = round(linspace(min_temp, max_temp, ncbtcks));
cb = colorbar(panel_c, 'Ticks', cbtck,...
    'TickLabels', cbtcklab,...
    'Units', 'normalized',...
    'Position', [0.76900673779618,0.12,0.168318402136886,0.018606701940035],...
    'Location', 'southoutside');
colormap(colors)

cb.Title.String = 'Temperature (°C)';
cb.Title.Units = 'normalized';
cb.Title.Position = [0.5 -2 0];

cb.Title.VerticalAlignment = 'bottom';
cb.Title.HorizontalAlignment = 'center';
cb.Title.FontName = 'Arial';
cb.Title.FontSize = fs;

%% Save figure
set(fig2,...
    'OuterPosition', 1000*[1.6757   -0.2063    1.1106    0.9273],...
    'Renderer', 'painters')

saveas(fig2, 'figure_2_20240430.svg')
