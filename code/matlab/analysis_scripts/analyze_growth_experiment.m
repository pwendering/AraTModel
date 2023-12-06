% Read experimental data from growth experiment with T-DNA lines

file_path = 'dry-weights-tdna-lines';
file_name = fullfile(file_path, '17_leaf_development_DATA_GA_combined.xlsx');

sheets = sheetnames(file_name);
sheets = cellstr(sheets(contains(sheets, 'dataRAW')));
numhead = 10;
sdt_str = 'Normal conditions';
edt_str = 'Move to';
date_pt = '\d+\.\d+(\.\d+)?';
outfmt = '%d\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%s\t%s\n';
rgr_fun = 'a*exp(b*t)';
raw_tab = importdata(file_name);

result_file = fullfile(file_path, 'rgr_measurements_20231005.txt');
fid = fopen(result_file, 'w');
fprintf(fid, 'Set\tGene\tRGR\tR2adj\tRMSE\tDW_final\tFW_final\teffect\tcomment\n');
for i = 1:numel(sheets)
    
    textdata = raw_tab.textdata.(sheets{i});
    start_date = regexp(textdata(any(contains(textdata, sdt_str), 2), :),...
        date_pt, 'match');
    start_date = [start_date{:}];
    
    tab = readtable(file_name, 'Sheet', sheets{i}, 'NumHeaderLines', numhead);
    
    rdcol = startsWith(tab.Properties.VariableNames, 'RD');
    dwcol = startsWith(tab.Properties.VariableNames, {'DW', 'DryWeight'});
    fwcol = startsWith(tab.Properties.VariableNames, {'FW', 'FreshWeight'});
    for j = 1:size(tab, 1)
        
        if ~isempty(tab.Gene{j})
            gene = tab.Gene{j};
            dw_meas = table2array(tab(j, dwcol));
            if ~isnumeric(dw_meas)
                dw_meas = str2double(dw_meas);
            end
            fw_meas = table2array(tab(j, fwcol));
            if ~isnumeric(fw_meas)
                fw_meas = str2double(fw_meas);
            end
            comment = tab.Comments{j};
            effect = tab.Effect{j};
        end
        
        rd_meas = table2array(tab(j, rdcol));
        if ~isnumeric(rd_meas)
            rd_meas = str2double(rd_meas);
        end
        
        % get final age in hours
        agestr = regexp(tab.Properties.VariableNames(rdcol),...
            strrep(date_pt, '\.', '_'), 'match');
        agestr = strrep([agestr{:}]', '_', '.');
        agestr(~endsWith(agestr, {'.23', '.2023'})) = strcat(agestr(~endsWith(agestr, {'.23', '.2023'})), '.23');
        agestr = strrep(agestr, '.23', '.2023');
        age = datenum(datetime(agestr) - start_date)*24;
        
        if all(rd_meas) && ~any(isnan(rd_meas))
            % fit exponential function to RD measurements
            [ffit, gof, output] = fit(age, rd_meas', fittype(rgr_fun, 'Independent', 't'),...
                'StartPoint', [1 0.01]);
            r2adj = gof.adjrsquare;
            rmse = gof.rmse;
            
            rgr = ffit.b;
        else
            r2adj = NaN;
            rmse = NaN;
            rgr = NaN;
        end
        
        fprintf(fid, outfmt, i, gene, rgr, r2adj, rmse, dw_meas,...
            fw_meas, effect, comment);
    end
end
fclose(fid);

%% Statistics
clear;clc;close all

file_path = 'dry-weights-tdna-lines';
file_name = fullfile(file_path, '17_leaf_development_DATA_GA_combined.xlsx');

% read clean table
rgr_tab = readtable(fullfile(file_path, 'rgr_measurements_20231005.txt'),...
    'Delimiter', '\t');

% remove samples with comment
rem_col_idx = ~cellfun(@isempty, rgr_tab.comment);
rgr_tab(rem_col_idx, :) = [];

% read SI table to match expected phenotypes
lookup_tab = readtable(fillfile(file_path, 'si_table_ko_predictions.xlsx'),...
    'Sheet', 'SI Table');
effects = cellfun(@(x)unique(lookup_tab.expectedPhenotype(ismember(lookup_tab.GeneID, x))),...
    rgr_tab.Gene, 'un', 0);
effects(ismember(rgr_tab.Gene, 'WT')) = {'WT'};

% fix cell/cellstr
for i = 1:numel(effects)
    % if both reduced growth and negative control are assigned to a gene,
    % assign only reduced growth
    if any(contains(effects{i}, 'Reduced growth at 17 °C')) && ...
            any(contains(effects{i}, 'negative control'))
        effects{i} = 'Reduced growth at 17 °C';
        warning('Multiple effects found for gene %s', rgr_tab.Gene{i})
    end
    
    if ~iscellstr(effects(i))
        if isempty(effects{i})
            effects(i) = {''};
        else
            effects(i) = effects{i};
        end
    end
end
rgr_tab.effect = effects;

% build unique table based on columns "Set", "Gene", and "DW_final"
col_idx = ismember(rgr_tab.Properties.VariableNames,...
    {'Set', 'Gene', 'DW_final'});
[~, ia] = unique(rgr_tab(:, col_idx), 'stable');
rgr_tab_uniq = sortrows(rgr_tab(ia, :), 'effect');

% print number of samples per category
unique_effects = unique(rgr_tab_uniq.effect);
for i = 1:numel(unique_effects)
    genes = rgr_tab_uniq.Gene(ismember(rgr_tab_uniq.effect, unique_effects(i)));
    fprintf('%d samples with category "%s" (%d genotypes)\n',...
        numel(genes),...
        unique_effects{i},...
        numel(unique(genes)))
end

% number of replicates per gene
n_reps = cellfun(@(x)sum(ismember(rgr_tab_uniq.Gene, x)), rgr_tab_uniq.Gene);

% remove genes with 2 or less replicates
rgr_tab_uniq(n_reps<=2, :) = [];

% remove genes without effect
rgr_tab_uniq(cellfun(@isempty, rgr_tab_uniq.effect), :) = [];

rgr_tab_uniq = sortrows(rgr_tab_uniq, 'effect', 'ascend');

%% Plot dry weights per genotype and experiment / batch
effect_uniq = unique(rgr_tab_uniq.effect);
ko_uniq = setdiff(rgr_tab_uniq.Gene, 'WT', 'stable');

% colors for predicted effects
colors = [[205, 86, 168]; [98, 182, 80]; [106, 140, 205]]/255;

% create a joint matrix for plotting, padded with nan
plot_matrix = nan(10, numel(ko_uniq)*4);
col_idx = 0;
box_colors = zeros(size(plot_matrix, 2), 3);
x_tick_labels = repmat({''}, size(plot_matrix, 2), 1);

for i = 1:numel(ko_uniq)
    
    % index of current KO line gene ID
    ko_idx = ismember(rgr_tab_uniq.Gene, ko_uniq(i));
    
    % sets / batches that the current line was grown in
    sets = unique(rgr_tab_uniq.Set(ko_idx));
    
    % predicted effect of current KO line
    tmp_effect = unique(rgr_tab_uniq.effect(ko_idx));
    effect_idx = ismember(effect_uniq, tmp_effect);
    
    % create separate columns for each set per gene
    for j = 1:numel(sets)
        
        % index of current set across all data
        set_idx = ismember(rgr_tab_uniq.Set, sets(j));
        
        % wild type values in current set
        tmp_wt_idx = ismember(rgr_tab_uniq.Gene, 'WT') & set_idx;
        tmp_wt_val = rgr_tab_uniq.DW_final(tmp_wt_idx);
        
        % indices of current KO line gene ID and set combination
        tmp_ko_idx = ko_idx & set_idx;
        tmp_ko_val = rgr_tab_uniq.DW_final(tmp_ko_idx);
        
        % update arrays with wild type values
        col_idx = col_idx + 1;
        plot_matrix(1:numel(tmp_wt_val), col_idx) = tmp_wt_val;
        box_colors(col_idx, :) = colors(ismember(effect_uniq, 'WT'), :);
        
        % update arrays with KO values
        col_idx = col_idx + 1;
        plot_matrix(1:numel(tmp_ko_val), col_idx) = tmp_ko_val;
        box_colors(col_idx, :) = colors(effect_idx, :);
        x_tick_labels{col_idx} = ko_uniq{i};
    end
    
    % increase column counter after each KO line to generate breaks between
    % boxes
    col_idx = col_idx + 1;
end

% remove empty columns on the right side of the matrix, keep empty columns
% in between for breaks between boxes
rem_col = find(any(~isnan(plot_matrix)), 1, 'last') + 1;
box_colors(rem_col:end, :) = [];
x_tick_labels(rem_col:end) = [];
plot_matrix(:, rem_col:end) = [];

% initialize figure
fig = figure;
ax = axes;

% create boxplot
boxplot(ax, plot_matrix,...
    'Colors', box_colors,...
    'Labels', x_tick_labels,...
    'Widths', .8)
% rotate x-labels
xtickangle(45)
% update x tick label positions
xtpos = get(ax, 'XTick');
set(ax, 'XTick', xtpos-0.5)
% add gene IDs in the middle between KO and WT boxes
% color boxed according to the predicted effect
h = findobj(ax, 'Tag', 'Box');
box_colors_sort = box_colors(size(box_colors, 1):-1:1, :);
for j = 1:numel(h)
    patch(get(h(j), 'XData'), get(h(j), 'YData'), box_colors_sort(j, :),...
        'FaceAlpha', .5);
end
hold(ax, 'on')
% add data points on top of boxes
plot(ax, plot_matrix',...
    'o',...
    'MarkerFaceColor', [.7 .7 .7],...
    'MarkerEdgeColor', [.4 .4 .4],...
    'MarkerSize', 3)
hold on

% create dummy line graphs with unique colors associated with predicted
% effects
h1 = plot(ax, NaN, NaN, 'LineWidth', 6, 'Color', [colors(1, :) 0.5]);
h2 = plot(ax, NaN, NaN, 'LineWidth', 5, 'Color', [colors(2, :) 0.5]);
h3 = plot(ax, NaN, NaN, 'LineWidth', 5, 'Color', [colors(3, :) 0.5]);

% add legend
legend([h1 h2 h3],...
    {'predicted effect', 'wild type', 'predicted no effect'},...
    'box', 'off', 'location', 'northoutside', 'NumColumns', numel(effect_uniq),...
    'AutoUpdate', 'off')

% label and axis format
ylabel(ax, 'dry weight (g)')
set(ax, 'FontSize', 9,...
    'FontName', 'Arial',...
    'TickLength', [0 0])

% set figure shape and size
set(fig, 'OuterPosition', [688.3333  219.0000  467.3333  445.3333])

% save figure
saveas(fig, 'dw_ko_lines_17deg.svg')

%{

% matthews correlation coefficient
mccfun = @(tp,tn,fp,fn)(tp*tn-fp*fn) / sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));

%}

%% write SI table with FW and DW
wt_idx = find(ismember(rgr_tab_uniq.Gene, 'WT'));
rgr_tab_uniq.Gene(wt_idx) = arrayfun(@(i)sprintf('WT_%i',...
    rgr_tab_uniq.Set(i)), wt_idx, 'un', 0);
ko_uniq = unique(rgr_tab_uniq.Gene, 'stable');

fw_dw_mat = nan(numel(ko_uniq), 2*7);
set_id = zeros(numel(ko_uniq), 1);
for i = 1:numel(ko_uniq)
    gene_idx = ismember(rgr_tab_uniq.Gene, ko_uniq(i));
    fw_dw_mat(i, 1:sum(gene_idx)) = rgr_tab_uniq.FW_final(gene_idx);
    fw_dw_mat(i, 8:7+sum(gene_idx)) = rgr_tab_uniq.DW_final(gene_idx);
    set_id(i) = unique(rgr_tab_uniq.Set(gene_idx));
end
    
fw_dw_table = array2table(fw_dw_mat, 'RowNames', ko_uniq,...
    'VariableNames', [arrayfun(@(i)sprintf('FW_R%i', i), 1:7, 'un', 0);...
    arrayfun(@(i)sprintf('DW_R%i', i), 1:7, 'un', 0)]);

writetable(fw_dw_table, 'Table_SX_ko_dw.xlsx',...
    'WriteRowNames', true)
