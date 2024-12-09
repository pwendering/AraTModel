% Perform sensitivity analysis

%% Load model
tmp = load(config('tgemFile'));
model = tmp.TGEM;

%% Define variables

sa_params = {'V_c_max', 'K_c', 'K_o', 'k_c', 'k_o', 'phi', 'J_max', ...
    'g_s', 'g_m'};

%% Sensitivity analysis I at a specific temperature
T = celsius2kelvin([10 25 40]);
I = 150;

sol_default = simulateTempEffects(model, I, 'tempRange', T);
rgr_default = sol_default.mu;
anet_default = sol_default.A_max;

sa_pct = 0.1;
rgr_decr = zeros(numel(sa_params), numel(T));
rgr_incr = zeros(numel(sa_params), numel(T));
anet_decr = zeros(numel(sa_params), numel(T));
anet_incr = zeros(numel(sa_params), numel(T));

for i = 1:numel(sa_params)

    for j = 1:numel(T)

        % -10%
        sol = simulateTempEffects(model, I, 'tempRange', T(j), ...
            'saBool', true, 'saParameter', sa_params{i}, 'saPercentage', -sa_pct);
        rgr_decr(i, j) = sol.mu / rgr_default(j);
        anet_decr(i, j) = sol.A_max / anet_default(j);

        % +10%
        sol = simulateTempEffects(model, I, 'tempRange', T(j), ...
            'saBool', true, 'saParameter', sa_params{i}, 'saPercentage', sa_pct);
        rgr_incr(i, j) = sol.mu / rgr_default(j);
        anet_incr(i, j) = sol.A_max / anet_default(j);
    end
end

%% plot results

tick_labels = regexprep(sa_params, '_(?<p>([\w_]+))', '_{$<p>}');
tick_labels = regexprep(tick_labels, '(?<=\{[\w_]+)_', '');

c = 0;
letters = 'ABCDEF';

fig_file_name = 'rgr_anet_sensitivity.png';
fig1 = figure;
t = tiledlayout(2, 3, 'TileSpacing', 'compact');

% RGR
for i = 1:numel(T)
    nexttile

    line([1 numel(sa_params)], [0 0], 'LineStyle', '--')

    hold on
    arrayfun(@(j)line([j j], [0 100*(rgr_decr(j, i)-1)], ...
        'LineWidth', 6, 'Color', 'r'), 1:numel(sa_params))

    hold on
    arrayfun(@(j)line([j j], [0 100*(rgr_incr(j, i)-1)], ...
        'LineWidth', 6, 'Color', 'k'), 1:numel(sa_params))

    if i == 1
        ylabel('Change in RGR (%)')
    end

    xlim([0.5 numel(sa_params)+0.5])
    ylim([-6 6])
    
    xticks([])

    c = c + 1;
    text(0.05, 0.93, letters(c), 'FontSize', 12, 'FontWeight', 'bold', ...
        'Units', 'normalized')

    text(0.05, 0.1, sprintf('%d °C', kelvin2celsius(T(i))), 'FontSize', 12, 'FontWeight', 'bold', ...
        'Units', 'normalized')

    set(gca, 'FontSize', 12)
end

% Anet
for i = 1:numel(T)
    nexttile

    line([1 numel(sa_params)], [0 0], 'LineStyle', '--')

    hold on
    arrayfun(@(j)line([j j], [0 100*(anet_decr(j, i)-1)], ...
        'LineWidth', 6, 'Color', 'r'), 1:numel(sa_params))

    hold on
    arrayfun(@(j)line([j j], [0 100*(anet_incr(j, i)-1)], ...
        'LineWidth', 6, 'Color', 'k'), 1:numel(sa_params))

    if i == 1
        ylabel('Change in A (%)')
    end

    xlim([0.5 numel(sa_params)+0.5])
    ylim([-6 6])

    c = c + 1;
    text(0.05, 0.93, letters(c), 'FontSize', 12, 'FontWeight', 'bold', ...
        'Units', 'normalized')

    text(0.05, 0.1, sprintf('%d °C', kelvin2celsius(T(i))), 'FontSize', 12, 'FontWeight', 'bold', ...
        'Units', 'normalized')

    set(gca, 'FontSize', 12)
    
    xticks(1:numel(sa_params))
    xticklabels(tick_labels)
    xtickangle(90)
end

h = zeros(2, 1);
h(1) = scatter(NaN, NaN, 80, 'r', 'filled', 'Marker', 's');
h(2) = scatter(NaN, NaN, 80, 'k', 'filled', 'Marker', 's');
legend(h, {'10% decrease', '10% increase'}, ...
    'box', 'off',...
    'position', [0.682976190476191 0.0113095249448503 0.256547619047619 0.0964285714285714])

t.XLabel.String = 'Parameter';
t.XLabel.FontSize = 14;

set(fig1, 'OuterPosition', [170.3333   90.3333  859.3333  507.3333])
saveas(fig1, fig_file_name, 'png')

%% Robustness analysis across temperature range between 10 °C and 40 °C

% temperature range
T = celsius2kelvin(linspace(10, 40, 10));
% light intensity
I = 150;

% percentage deviation for simulating standard deviation around the mean
sa_pct = 0.05;

% number of random samples
n_samples = 1000;

% define parameters
sa_params = {'V_c_max', 'K_c', 'K_o', 'k_c', 'k_o', 'J_max'};

% simulated parameter distributions
params_mean = cellfun(@(x)config(x), sa_params);
params_log_mean = log(params_mean);
params_log_sd = sqrt(log(1+(sa_pct*params_mean).^2./params_mean.^2));

param_samples_log = cell2mat(arrayfun(@(i) ...
    normrnd(params_log_mean(i), params_log_sd(i), n_samples, 1), ...
    1:numel(sa_params), 'un', 0));

param_samples_lin = exp(param_samples_log);

% Simulate default temperature response
sol = simulateTempEffects(model, I, 'tempRange', T);
rgr_default = sol.mu;
anet_default = sol.A_max;

% Simulate temperature responses over the temperature range for each random
% parameter point
rgr_sol = nan(n_samples, numel(T));
anet_sol = nan(n_samples, numel(T));

model_tmp = model;
for i = 1:n_samples

    % set parameter values
    for j = 1:numel(sa_params)
        model_tmp.(sa_params{j}) = param_samples_lin(i, j);
    end

    % simulate temperature response
    try
        sol_tmp = simulateTempEffects(model_tmp, I, 'tempRange', T);
        rgr_sol(i, :) = sol_tmp.mu;
        anet_sol(i, :) = sol_tmp.A_max;
    catch
        fprintf('Parameter point %i is invalid.\n', i)
    end
    
    if mod(i, 100) == 0
        save('sensitivity_analysis_original_gamma.mat')
    end
end

%% Plot results

fig2 = figure;
t = tiledlayout(1, 2, 'TileSpacing', 'compact');
fig_file_name = 'anet_rgr_robustness.png';

% RGR
nexttile
hold on
arrayfun(@(i)plot(kelvin2celsius(T), rgr_sol(i, :), 'Color', [.8 .8 .8]), 1:n_samples)
plot(kelvin2celsius(T), rgr_default, ...
    'LineWidth', 1.5, 'Color', 'k')
ylabel('RGR (h^{-1})')
text(0.05, 0.93, 'A', 'FontSize', 14, 'FontWeight', 'bold', 'Units', 'normalized')
set(gca, 'FontSize', 14)

errorbar(kelvin2celsius(T), rgr_default, std(rgr_sol), std(rgr_sol), ...
    'Color', 'k')

% Anet
nexttile
hold on
arrayfun(@(i)plot(kelvin2celsius(T), anet_sol(i, :), 'Color', [.8 .8 .8]), 1:n_samples)
plot(kelvin2celsius(T), anet_default, ...
    'LineWidth', 1.5, 'Color', 'k')
ylabel({'A (\mumol m^{-2} s^{-1})'})
text(0.05, 0.93, 'B', 'FontSize', 14, 'FontWeight', 'bold', 'Units', 'normalized')
set(gca, 'FontSize', 14)

errorbar(kelvin2celsius(T), anet_default, std(anet_sol), std(anet_sol), ...
    'Color', 'k')

% common x-label
t.XLabel.String = 'Temperature (°C)';
t.XLabel.FontSize = 14;

set(fig2, 'OuterPosition', [59.0000  252.3333  765.3333  388.0000])
saveas(fig2, fig_file_name, 'png')
