%% Read input workspace
clear; clc
load sampling_workspace_I_150_20230817_constr_sum_of_flux.mat

% reaction IDs
reactions = model.rxns(~startsWith(model.rxns, {'draw_', 'prot_'}));
reactions(ismember(reactions, 'Bio_opt')) = {'GROWTH'};
reactions(ismember(reactions, 'Im_hnu')) = {'PHOTON_UPTAKE'};
reactions(ismember(reactions, 'Im_H2O_REV')) = {'TRANSPIRATION'};
n_rxns = numel(reactions);

% name of workspace to store sampling results
c = clock;
datestring = sprintf('%d%02.0f%02.0f', c(1:3));

%% 1) Plot pathway fluxes

model.subSystems(cellfun(@isempty, model.subSystems)) = {''};
n_asm_idx = cellfun(@(x)contains(x,...
    {'nitrate assimilation', 'glutamate synthesis', 'glutamine synthesis'}), model.subSystems);
model.subSystems(n_asm_idx) = strcat(model.subSystems(n_asm_idx), ', N assimilation');

% define subsystems that should be plotted
subsystems = {
 'Calvin-Benson cycle'
    'sucrose synthesis'
    'starch synthesis'
    'N assimilation'
    'light reactions'
    'photorespiration'
    'starch degradation'
    'pentose phosphate pathway'
    'gluconeogenesis'
    'fatty acid synthesis'
    'glycolysis'
    'tricarboxylic acid cycle'
    'Shikimate pathway'
    'nucleotide metabolism'
    'oxidative phosphorylation'
    'sulfur assimilation'
    };

% subsystems = {
%     'glutamate synthesis'
%     'tyrosine synthesis'
%     'aspartate synthesis'
%     'glycine synthesis'
%     'leucine synthesis'
%     'lysine synthesis'
%     'methionine synthesis'
%     'serine synthesis'
%     'threonine synthesis'
%     'valine synthesis'
%     'alanine synthesis'
%     'arginine synthesis'
%     'asparagine synthesis'
%     'cysteine synthesis'
%     'glutamine synthesis'
%     'histidine synthesis'
%     'isoleucine synthesis'
%     'phenylalanine synthesis'
%     'proline synthesis'
%     'tryptophan synthesis'
%     };

% Determine statistics of pathway flux
iqr_pf = zeros(numel(subsystems), numel(temperatures));
median_pf = zeros(numel(subsystems), numel(temperatures));
av_pf = zeros(numel(subsystems), numel(temperatures));
w_min_pf = zeros(numel(subsystems), numel(temperatures));
w_max_pf = zeros(numel(subsystems), numel(temperatures));
q25_pf = zeros(numel(subsystems), numel(temperatures));
q75_pf = zeros(numel(subsystems), numel(temperatures));
for i = 1:numel(subsystems)
    
    rxns = model.rxns(cellfun(@(x)contains(char(x), subsystems{i}), model.subSystems));
    
    flux_sums = cell2mat(arrayfun(@(tp)sum(full(samples{tp}(ismember(tmp_qcp.varnames, rxns), :))),...
        1:numel(temperatures), 'un', 0)')';
    
    q25 = quantile(flux_sums, 0.25, 1);
    q75 = quantile(flux_sums, 0.75, 1);
    iqrange = q75-q25;
    
    iqr_pf(i, :) = iqrange;
    q25_pf(i, :) = q25;
    q75_pf(i, :) = q75;
    w_max_pf(i, :) = min(q75+1.5*iqrange, max(flux_sums));
    w_min_pf(i, :) = max(q25-1.5*iqrange, min(flux_sums));
    median_pf(i, :) = median(flux_sums, 1);
    av_pf(i, :) = mean(flux_sums, 1);
end
clear q25 q75 iqrange flux_sums rxns

% Plot distributions of pathway flux per subsystem
figure
t = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
line_color = [0 0 0];
for i = 1:numel(subsystems)
    
    nexttile
    hold on
    
    % "whiskers"
    x = kelvin2celsius([temperatures(1); temperatures(1); temperatures(2:end)';...
        temperatures(end:-1:2)']);
    y = [w_max_pf(i, 1); w_min_pf(i, 1); w_min_pf(i, 2:end)'; w_max_pf(i, end:-1:2)'];
    F = fill(x, y, 'k--');
    F.FaceColor = [.8 .8 .8];
    F.EdgeColor = [.4 .4 .4];
    
    % quartiles
    x = kelvin2celsius([temperatures(1); temperatures(1); temperatures(2:end)';...
        temperatures(end:-1:2)']);
    y = [q75_pf(i, 1); q25_pf(i, 1); q25_pf(i, 2:end)'; q75_pf(i, end:-1:2)'];
    F = fill(x, y, 'k');
    F.FaceColor = [.6 .6 .6];
    F.EdgeColor = [.4 .4 .4];
    
    % median
    plot(kelvin2celsius(temperatures),...
        median_pf(i, :),...
        'Color', line_color,...
        'LineWidth', 2)
    
    set(gca, 'YScale', 'linear',...
        'FontName', 'Arial',...
        'XColor', 'k',...
        'YColor', 'k')
       
    title(subsystems{i}, 'FontSize', 10)
    
end

t.XLabel.String = 'Temperature (°C)';
t.XLabel.FontName = 'Arial';
t.XLabel.FontSize = 16;

t.YLabel.String = 'Pathway flux (mmol gDW^{-1} h^{-1})';
t.YLabel.FontName = 'Arial';
t.YLabel.FontSize = 16;
