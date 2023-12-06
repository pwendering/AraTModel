% Create boxplots of operational flux ranges

clear;clc;close all

% load colormap
tmp = load('blue_red_colormap.mat');
colors = tmp.blueredcm;
min_temp = 6;
max_temp = 42;
n_colors = max_temp-min_temp+1;
clear tmp

fs = 10;
bx_lw = 1.3;

load sampling_workspace_I_150_20230817_constr_sum_of_flux.mat min_fluxes ...
    max_fluxes temperatures irradiance

op_ranges = max_fluxes-min_fluxes;
op_ranges_log10 = log10(op_ranges);
op_ranges_log10(op_ranges_log10<=0) = NaN;
op_ranges_log10 = log10(op_ranges_log10);

y_min = floor(min(min(op_ranges_log10)));
y_max = ceil(max(max(op_ranges_log10)));

% operational ranges
[~,~,bins] = histcounts(kelvin2celsius(temperatures), n_colors, 'BinLimits', [min_temp max_temp]);
bxcolor = colors(bins, :);
for i = 1:size(op_ranges_log10, 2)
    bc = boxchart(categorical(repelem(i, size(op_ranges_log10, 1), 1)), op_ranges_log10(:, i),...
        'MarkerStyle', '.',...
        'JitterOutliers', 'on',...
        'MarkerColor', 'k',...
        'BoxFaceAlpha', .5);
    hold on
    bc.BoxFaceColor = bxcolor(i, :);
end
xticklabels(round(kelvin2celsius(temperatures)))
xtickangle(90)
xlabel('Temperature (°C)')

ylabel({'log_{10} operational flux ranges', '(mmol gDW^{-1} h^{-1})'})

set(gca,...
    'Box', 'off',...
    'LineWidth', bx_lw,...
    'FontName', 'Arial',...
    'FontSize', fs,...
    'YLim', [y_min y_max],...
    'XColor', 'k',...
    'YColor', 'k')
set(gcf, 'OuterPosition', [353.6667  327.6667  376.0000  371.3333])
exportgraphics(gcf, ['operational_ranges_I_' num2str(irradiance) '.png'],...
    'Resolution', 300)
