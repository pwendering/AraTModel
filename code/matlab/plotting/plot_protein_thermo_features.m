% Plot general statistics of protein thermostability

%% 1) Read the model

% read model
model_file = config('tgemFile');
tmp = load(model_file);
model = tmp.TGEM;
clear tmp model_file

% add RubisCO activase module
model = addRCA(model);

% protein UniProt IDs
proteins = erase(...
    model.mets(...
    startsWith(model.mets,'prot_') &...
    ~contains(model.mets,'pool')),...
    'prot_');

%% 2) Melting temperatures
figure
t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile(1)

T_m = model.T_m;
fprintf('Melting temperatures from Meltome Atlas for %.2f%% (%d/%d)\n',...
    100*sum(~isnan(T_m))/numel(T_m), sum(~isnan(T_m)), numel(T_m))

histogram(T_m(~isnan(T_m)))
xlabel('T_m (°C)')
ylabel('Count')

text(0.08, 0.9, '(a)',...
    'Fontname', 'Arial', 'FontSize', 14, 'fontweight', 'bold', 'Units', 'normalized')

set(gca,...
    'FontSize', 12,...
    'Box', 'on',...
    'LineWidth', 1)

%% 3) Optimal temperatures
nexttile(2)

T_opt = model.T_opt;
meltome_idx = model.T_H==70.4;
fprintf('Optimal temperatures from Meltome Atlas for %.2f%% (%d/%d)\n',...
    100*sum(meltome_idx)/numel(T_opt), sum(meltome_idx), numel(T_opt))

histogram(T_opt(meltome_idx),...
    'NumBins', 20,...
    'BinLimits', [0 55])
hold on
histogram(T_opt(~meltome_idx),...
    'NumBins', 20,...
    'BinLimits', [0 55])

xticks([0 10 20 30 40 50])
xlabel('T_{opt} (°C)')
ylabel('Count')

text(0.08, 0.9, '(b)',...
    'Fontname', 'Arial', 'FontSize', 14, 'fontweight', 'bold', 'Units', 'normalized')

set(gca,...
    'FontSize', 12,...
    'Box', 'on',...
    'LineWidth', 1)

legend({'fit', 'pred.'},...
    'Location', 'sw',...
    'Box', 'off',...
    'FontSize', 12)

%% 4) Key temperatures with fit quality
meltingCurves = struct;
[proteinIDs, meltingCurves.f, estimates, rmse, T_m, T_opt, T_H, fN_TH, r2adj]= ...
    fitTPPBetaFcn(config('meltomeFile'), proteins, false);

% by R2
nexttile(4)

r2adj(isnan(T_m)|isnan(T_opt)) = NaN;

n_bins = 5;
[bins, edges] = discretize(r2adj, n_bins);

colormap(parula(n_bins))
scatter(T_m, T_opt, 'filled', 'CData', bins)
cb = colorbar('Ticks', 1:.799:numel(edges),...
    'TickLabels', arrayfun(@(i)sprintf('%.2f', i), edges, 'un', 0));
cb.Title.String = 'R^2_{adj}';

xlabel('T_m (°C)')
ylabel('T_{opt} (°C)')

text(0.08, 0.9, '(c)',...
    'Fontname', 'Arial', 'FontSize', 14, 'fontweight', 'bold', 'Units', 'normalized')

set(gca,...
    'FontSize', 12,...
    'Box', 'on',...
    'LineWidth', 1)

% whole meltome
nexttile(5)
load meltome_ml_ws_20230531.mat

[bins, edges] = discretize(r2adj(~isnan(T_m)), n_bins);

colormap(parula(n_bins))
scatter(T_m(~isnan(T_m)), T_opt(~isnan(T_m)), 'filled',...
    'CData', bins)
cb = colorbar('Ticks', 1:.799:numel(edges),...
    'TickLabels', arrayfun(@(i)sprintf('%.2f', i), edges, 'un', 0));
cb.Title.String = 'R^2_{adj}';

xlabel('T_m (°C)')
ylabel('T_{opt} (°C)')

text(0.08, 0.9, '(d)',...
    'Fontname', 'Arial', 'FontSize', 14, 'fontweight', 'bold', 'Units', 'normalized')

set(gca,...
    'Box', 'on',...
    'LineWidth', 1,...
    'FontSize', 12)

%% 5) kcats
tempRange = linspace(10, 40, 10);

kcat_mmrt = cell(1,numel(tempRange));

enzyme_row_idx = find(startsWith(model.mets,'prot_')&~contains(model.mets,'prot_pool'));

for i=1:numel(tempRange)
    tmpModel = adjustKcatsMMRTGECKO(model,celsius2kelvin(tempRange(i)));
    adj_kcats = arrayfun(@(i)-1./tmpModel.S(i,model.S(i,:)<0)/3600,enzyme_row_idx,'un',0);
    kcat_mmrt{i} = [adj_kcats{:}]';
end
clear tmpModel adj_kcats

kcat_mmrt = cell2mat(kcat_mmrt);

nexttile(3, [2 1])

rng(42)
x = 0.3*rand(1, numel(tempRange)*size(kcat_mmrt,1))-0.15 + repmat(1:numel(tempRange),1,size(kcat_mmrt,1));
y = reshape(kcat_mmrt',numel(kcat_mmrt),1);

scatter(x, y, 8, 'filled',...
    'MarkerFaceColor', [.6 .6 .6],...
    'MarkerEdgeColor', [.2 .2 .2])
hold on
boxchart(repmat(1:numel(tempRange),1,size(kcat_mmrt,1)), y,...
    'BoxFaceColor', [1 .5 .5],...
    'WhiskerLineColor', [1 .5 .5],...
    'LineWidth', 2,...
    'MarkerStyle', 'none')

% x-axis
xticks(1:numel(tempRange))
xticklabels(round(tempRange))
xtickangle(90)
xlabel('Temperature (°C)')

% y-axis
yticks(logspace(-5, 5, 6))
ylabel('log_{10} k_{cat} (s^{-1})')

text(0.08, 0.96, '(e)',...
    'Fontname', 'Arial', 'FontSize', 14, 'fontweight', 'bold', 'Units', 'normalized')

set(gca,...
    'LineWidth', 1,...
    'FontSize', 12,...
    'YScale','log',...
    'Box', 'on')

set(gcf, 'OuterPosition', 1000*[-1.4650    0.1390    1.1433    0.6987])

exportgraphics(gcf, 'aracore_tm_topt_kcat.png')

%% 6) Plot MMRT parameters
DCp = model.DCp;
DH = model.DH;
DS = model.DS;

tiledlayout(2,3)

letter_fs = 14;
letter_xpos = -0.05;
letter_ypos = 1.1;

nexttile
histogram(DCp(DCp~=0), 'FaceColor', [.9 .9 .9], 'NumBins', 20)
ylabel('Count', 'Interpreter', 'latex')
text(letter_xpos, letter_ypos, '(a)', 'FontName', 'Arial', 'FontWeight', 'bold',...
    'FontSize', letter_fs, 'Units', 'normalized')
set(gca, 'Box', 'off', 'FontSize', 12)

nexttile
histogram(DH(DH~=0), 'FaceColor', [.9 .9 .9], 'NumBins', 20)
text(letter_xpos, letter_ypos, '(b)', 'FontName', 'Arial', 'FontWeight', 'bold',...
    'FontSize', letter_fs, 'Units', 'normalized')
set(gca, 'Box', 'off', 'FontSize', 12)

nexttile
histogram(DS(DS~=0), 'FaceColor', [.9 .9 .9], 'NumBins', 20)
text(letter_xpos, letter_ypos, '(c)', 'FontName', 'Arial', 'FontWeight', 'bold',...
    'FontSize', letter_fs, 'Units', 'normalized')
set(gca, 'Box', 'off', 'FontSize', 12)

nexttile
scatter(DCp(DCp~=0), DS(DS~=0), 30, 'filled', 'MarkerFaceColor', [.9 .9 .9],...
    'MarkerEdgeColor', [.4 .4 .4])
xlabel('${\Delta}C_p^{\ddagger}$', 'Interpreter', 'latex')
ylabel('${\Delta}S_{T_0}^{\ddagger}$', 'Interpreter', 'latex')
text(letter_xpos, letter_ypos, '(d)', 'FontName', 'Arial', 'FontWeight', 'bold',...
    'FontSize', letter_fs, 'Units', 'normalized')
set(gca,'Box', 'off', 'FontSize', 12)

nexttile
scatter(DH(DH~=0), DCp(DCp~=0), 30, 'filled', 'MarkerFaceColor', [.9 .9 .9],...
    'MarkerEdgeColor', [.4 .4 .4])
xlabel('${\Delta}H_{T_0}^{\ddagger}$', 'Interpreter', 'latex')
ylabel('${\Delta}C_p^{\ddagger}$', 'Interpreter', 'latex')
text(letter_xpos, letter_ypos, '(e)', 'FontName', 'Arial', 'FontWeight', 'bold',...
    'FontSize', letter_fs, 'Units', 'normalized')
set(gca, 'Box', 'off', 'FontSize', 12)

nexttile
scatter(DS(DS~=0), DH(DH~=0), 30, 'filled', 'MarkerFaceColor', [.9 .9 .9],...
    'MarkerEdgeColor', [.4 .4 .4])
xlabel('${\Delta}S_{T_0}^{\ddagger}$', 'Interpreter', 'latex')
ylabel('${\Delta}H_{T_0}^{\ddagger}$', 'Interpreter', 'latex')
text(letter_xpos, letter_ypos, '(f)', 'FontName', 'Arial', 'FontWeight', 'bold',...
    'FontSize', letter_fs, 'Units', 'normalized')
set(gca, 'Box', 'off', 'FontSize', 12)

print('mmrt_parameters_a_th_gecko_raw_model', '-dpng', '-painters')
