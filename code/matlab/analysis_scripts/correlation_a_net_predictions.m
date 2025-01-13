% Compare predictions of the net CO2 assimilation rate with experimental data

clear;clc;close all

% read model
model_file = config('tgemFile');
tmp = load(model_file);
TGEM = tmp.TGEM;
clear tmp

% add RubisCO activase module
TGEM = addRCA(TGEM);

% block uptake of H2S and chloroplastic alanine transaminase
TGEM.ub(findRxnIDs(TGEM, {'Im_H2S', 'Im_H2S_REV'})) = 0;
TGEM.ub(findRxnIDs(TGEM, {'AlaTA_h', 'AlaTA_h_REV'})) = 0;
TGEM.ub(findRxnIDs(TGEM, {'Bio_AA', 'Bio_CLim', 'Bio_NLim'})) = 0;

% correct THF transport reaction
TGEM.S(findMetIDs(TGEM, 'H_h'), findRxnIDs(TGEM, {'Tr_THF_h', 'Tr_THF_h_REV'})) = [-1 1];

% load colormap
tmp = load('blue_red_colormap.mat');
colors = tmp.blueredcm;
min_temp = 5;
max_temp = 42;
n_colors = max_temp-min_temp+1;
clear tmp

tempRange = celsius2kelvin(linspace(min_temp, max_temp, 10));

%% Correlation of experimental vs predicted A

a_net_tab = readtable(fullfile(config('pathToTGEM'),...
    'photosynthesis-parameters', 'a_net_dataset.xlsx'));

% only select rows with first author IDs
a_net_tab = a_net_tab(~cellfun(@isempty, a_net_tab.Author), :);
no_imp_ca_idx = ~isnan(a_net_tab.pCa_ubar);
no_imp_ci_idx = ~isnan(a_net_tab.pCi_ubar);

% impute missing Ci and Ca values assuming that Ci = 0.75*C_a
% C_a = 1.33 * C_i
a_net_tab.pCa_ubar(isnan(a_net_tab.pCa_ubar)) = ...
    4/3*a_net_tab.pCi_ubar(isnan(a_net_tab.pCa_ubar));
% C_i = 0.75 * C_a
a_net_tab.pCi_ubar(isnan(a_net_tab.pCi_ubar)) = ...
    0.75*a_net_tab.pCa_ubar(isnan(a_net_tab.pCi_ubar));

% predict A
a_cbm = nan(size(a_net_tab, 1), 1);
a_far = nan(size(a_net_tab, 1), 1);
fvcb_params = getFvCBParams;
for i = 1:size(a_net_tab, 1)

    % constraint-based prediction
    tmp_model = TGEM;
    I = a_net_tab.Irradiance_mumol_m2_s(i);
    tmp_model.C_a = a_net_tab.pCa_ubar(i);
    T = celsius2kelvin(a_net_tab.Temperature_C(i));

    if ~isnan(a_net_tab.pO_mbar(i))
        tmp_model.O_a = a_net_tab.pO_mbar(i)*1000;
    end

    if ~isnan(I) && ~isnan(T) && ~isnan(tmp_model.C_a)
        sol = simulateTempEffects(tmp_model, I, 'tempRange', T);
        a_cbm(i) = sol.A_max;

        % FvCB model prediction
        a_far(i) = farquhar(...
            'T', T,...
            'C', a_net_tab.pCi_ubar(i),...
            'O', tmp_model.O_a,...
            'I', I,...
            'kc', fvcb_params.k_c,...
            'ko', fvcb_params.k_o,...
            'K_c', fvcb_params.K_c,...
            'K_o', fvcb_params.K_o,...
            'E_a', fvcb_params.E_a,...
            'V_c_max', fvcb_params.V_c_max,...
            'J_max', fvcb_params.J_max);
    end
end

save('a_net_pred_meas_workspace')

%% Plot comparison

clear; clc; close all

load a_net_pred_meas_workspace.mat

filter_idx = no_imp_ca_idx;

A_exp = a_net_tab.A_mumol_m2_s(filter_idx);
A_cbm = a_cbm(filter_idx);
A_far = a_far(filter_idx);
a_net_tab_red = a_net_tab(filter_idx, :);

a_net_tab_red.I_Ca_ratio = sqrt(a_net_tab_red.Irradiance_mumol_m2_s./a_net_tab_red.pCa_ubar);
a_net_tab_red.Ci_C_a_ratio = a_net_tab_red.pCi_ubar./a_net_tab_red.pCa_ubar;
a_net_tab_red.Study = strcat(a_net_tab_red.Author,...
    arrayfun(@(y)cellstr([' (' num2str(y) ')']), a_net_tab_red.Year));

fig = figure('Renderer', 'painters');
t = tiledlayout('flow', 'TileSpacing', 'compact');
cats = {'publication', 'T (°C)', {'C_a', '(\mubar)'}, {'irradiance', '(\mumol/m^2/s)'}};

headers = {'Study', 'Temperature_C', 'pCa_ubar', 'Irradiance_mumol_m2_s'};
msz = 20;
letters = char(65:65+25);
g_array = gobjects(numel(cats), 1);
cb_array = gobjects(numel(cats), 1);

% color blind friendly palette
cb_colors = [[140 168 239]; [120 94 240]; [220 38 127]; [254 97 0]; ...
    [255 176 0]; [159 194 106]]/255;
markers = 'os^v';
c = 0;
for i = 1:numel(cats)

    if isequal(headers{i}, 'Study')
        studies = categorical(a_net_tab_red.(headers{i}));
        X = grp2idx(studies);
        n_bins = numel(unique(X));
        colors = repmat(cb_colors, ceil(n_bins/size(cb_colors, 1)), 1);
    else
        X = a_net_tab_red.(headers{i});
        X(isnan(X)) = mean(X, 'omitnan');
        n_bins = 7;
        colors = cool(n_bins);
        markers = repelem('o', size(colors, 1), 1);
    end
    [~,edges,bins] = histcounts(X,'NumBins', n_bins);

    c = c + 1;
    g_array(c) = nexttile;
    lm = fitlm([0 1], [0 1]);
    x_limits = [min(A_exp) max(A_exp)];
    line(g_array(c), x_limits,...
        lm.predict(x_limits'),...
        'Color', 'k',...
        'LineWidth', 2)
    hold on

    if isequal(headers{i}, 'Study')
        plt1 = zeros(ceil(n_bins/size(cb_colors, 1)), 1);
        for j = 1:numel(plt1)
            idx = (j-1)*size(cb_colors, 1)+1:j*size(cb_colors, 1);
            plt1 = scatter(g_array(c), A_exp(ismember(bins, idx)), ...
                A_cbm(ismember(bins, idx)), msz, 'filled',...
                'CData', colors(bins(ismember(bins, idx)), :),...
                'MarkerEdgeColor', 'k', ...
                'Marker', markers(j));
        end
    else
        plt1 = scatter(g_array(c), A_exp, A_cbm, msz, 'filled',...
            'CData', colors(bins,:),...
            'MarkerEdgeColor', 'k');
    end
    text(g_array(c), 0.04, 0.92, letters(c),...
        'units', 'normalized',...
        'fontname', 'Arial',...
        'fontsize', 14,...
        'fontweight', 'bold')
    set(g_array(c), 'FontSize', 10, 'FontName', 'Arial')

    c = c + 1;
    g_array(c) = nexttile;
    lm = fitlm([0 1], [0 1]);
    x_limits = [min(A_exp) max(A_exp)];
    line(g_array(c), x_limits,...
        lm.predict(x_limits'),...
        'Color', 'k',...
        'LineWidth', 2)
    hold on
    if isequal(headers{i}, 'Study')
        plt2 = zeros(ceil(n_bins/size(cb_colors, 1)), 1);
        for j = 1:numel(plt2)
            idx = (j-1)*size(cb_colors, 1)+1:j*size(cb_colors, 1);
            plt2 = scatter(g_array(c), A_exp(ismember(bins, idx)), ...
                A_far(ismember(bins, idx)), msz, 'filled',...
                'CData', colors(bins(ismember(bins, idx)), :),...
                'MarkerEdgeColor', 'k', ...
                'Marker', markers(j));
        end
    else
        plt2 = scatter(g_array(c), A_exp, A_far, msz, 'filled',...
            'CData', colors(bins,:),...
            'MarkerEdgeColor', 'k');
    end
    text(g_array(c), 0.04, 0.92, letters(c),...
        'units', 'normalized',...
        'fontname', 'Arial',...
        'fontsize', 14,...
        'fontweight', 'bold')
    set(g_array(c), 'FontSize', 10, 'FontName', 'Arial')

    if ~isequal(headers{i}, 'Study')
        if isequal(headers{i}, 'Year')
            edges = round(edges);
        end
        cb_array(c) = colorbar(g_array(c), 'Ticks', linspace(0, 1, numel(edges)));
        colormap(g_array(c), colors)
        cb_array(c).TickLabels = edges;
        cb_array(c).Title.String = cats{i};
    else
        c = c + 1;
        g_array(c) = nexttile;

        h = zeros(n_bins, 1);
        for j = 1:numel(h)
            h(j) = plot(NaN, NaN, 's', 'MarkerFaceColor', colors(j,:),...
                'MarkerEdgeColor', 'k', 'MarkerSize', msz/3, ...
                'Marker', markers(ceil(j/size(cb_colors, 1))));
            hold on
        end
        studies_uniq = unique(a_net_tab_red.Study, 'stable');
        labels = studies_uniq;
        labels(unique(bins, 'stable')) = studies_uniq;
        legend(h, labels,...
            'box', 'off',...
            'numcolumns', 4,...
            'location', 'w')
        set(g_array(c), 'Color', 'none', 'XColor', 'none', 'YColor', 'none')

        c = c + 1;
        g_array(c) = nexttile;
        set(g_array(c), 'Color', 'none', 'XColor', 'none', 'YColor', 'none')
    end

    if i==1
        [rho_cbm, p_cbm] = corr(A_exp, A_cbm,...
            'rows', 'complete');
        [rho_far, p_far] = corr(A_exp, A_far,...
            'rows', 'complete');
        legend([plt1 plt2], ['ecAraCore r = ' num2str(rho_cbm, '%.2f')],...
            ['FvCB r = ' num2str(rho_far, '%.2f')],...
            'Location', 'se',...
            'box', 'off')
    end

end

t.XLabel.String = 'A measured (\mumol m^{-2} s^{-1})';
t.YLabel.String = 'A predicted (\mumol m^{-2} s^{-1})';

set(gcf, 'OuterPosition', 1000*[-1.4057   -0.2250    0.6960    1.0534])

for i = 1:numel(g_array)
    set(g_array(i), 'ylim', [-10 30])
end

saveas(gcf, 'a_net_exp_scatter.svg')

%% CV and correlation

% Pearson correlation
[r,p] = corr(a_net_tab_red.A_mumol_m2_s, a_net_tab_red.Irradiance_mumol_m2_s);
fprintf('Correlation between A and irradiance: %.3g (P=%.3g)\n', r, p)
[r,p] = corr(a_net_tab_red.A_mumol_m2_s, a_net_tab_red.Temperature_C);
fprintf('Correlation between A and temperature: %.3g (P=%.3g)\n', r, p)
[r,p] = corr(a_net_tab_red.A_mumol_m2_s, a_net_tab_red.pCa_ubar);
fprintf('Correlation between A and ambient p(CO2): %.3g (P=%.3g)\n', r, p)
[r,p] = corr(a_net_tab_red.A_mumol_m2_s, a_net_tab_red.pCi_ubar);
fprintf('Correlation between A and intercellular p(CO2): %.3g (P=%.3g)\n', r, p)
[r,p] = corr(a_net_tab_red.A_mumol_m2_s, a_net_tab_red.pO_mbar,...
    'rows', 'complete');
fprintf('Correlation between A and ambient p(O2): %.3g (P=%.3g)\n', r, p)
[r,p] = corr(a_net_tab_red.A_mumol_m2_s, a_net_tab_red.RH_pct,...
    'rows', 'complete');
fprintf('Correlation between A and relative humidity: %.3g (P=%.3g)\n', r, p)
[r,p] = corr(a_net_tab_red.A_mumol_m2_s, a_net_tab_red.Age_dag,...
    'rows', 'complete');
fprintf('Correlation between A and plant age: %.3g (P=%.3g)\n', r, p)

% coefficient of variation
cv_fun = @(x)std(x, 'omitnan')/mean(x, 'omitnan');
fprintf('A: %.3g\n', cv_fun(a_net_tab_red.A_mumol_m2_s))
fprintf('I: %.3g\n', cv_fun(a_net_tab_red.Irradiance_mumol_m2_s))
fprintf('pCa: %.3g\n', cv_fun(a_net_tab_red.pCa_ubar))
fprintf('T: %.3g\n', cv_fun(a_net_tab_red.Temperature_C))
fprintf('Age: %.3g\n', cv_fun(a_net_tab_red.Age_dag))
fprintf('RH: %.3g\n', cv_fun(a_net_tab_red.RH_pct))
fprintf('Ci: %.3g\n', cv_fun(a_net_tab_red.pCi_ubar))
fprintf('pO: %.3g\n', cv_fun(a_net_tab_red.pO_mbar))

%% specific differences in prediction quality for both models

% calculate residuals between predicted and experimental data
rmse_cbm = abs(A_exp-A_cbm);
rmse_far = abs(A_exp-A_far);

% compare differences to exp. data at high temperatures (> 30 °C)
high_t_idx = a_net_tab_red.Temperature_C > 30;

p = ranksum(rmse_cbm, rmse_far, 'tail', 'left');
fprintf('All data points: Ha: |A(cbm)-A(exp)| < |A(fvcb)-A(exp)| => P = %.3g\n',...
    p)
fprintf('Median RMSE CBM: %.3g, median FvCB: %.3g\n',...
    median(rmse_cbm), median(rmse_far))
p = ranksum(rmse_cbm(high_t_idx), rmse_far(high_t_idx), 'tail', 'left');
fprintf('T > 25 °C: Ha: |A(cbm)-A(exp)| < |A(fvcb)-A(exp)| => P = %.3g\n',...
    p)
fprintf('Median RMSE CBM: %.3g, median FvCB: %.3g\n',...
    median(rmse_cbm(high_t_idx)), median(rmse_far(high_t_idx)))

% compare predictions at high CO2 and high irradiance
high_co2_idx = a_net_tab_red.pCa_ubar > 600;
high_light_idx = a_net_tab_red.Irradiance_mumol_m2_s > 1000;

p = ranksum(rmse_cbm(high_co2_idx|high_light_idx),...
    rmse_far(high_co2_idx|high_light_idx), 'tail', 'left');
fprintf('p(CO2)>600, I>1000: Ha: |A(cbm)-A(exp)| < |A(fvcb)-A(exp)| => P = %.3g\n',...
    p)
fprintf('Median RMSE CBM: %.3g, median FvCB: %.3g\n',...
    median(rmse_cbm(high_co2_idx|high_light_idx)),...
    median(rmse_far(high_co2_idx|high_light_idx)))