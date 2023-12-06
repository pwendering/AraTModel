function [proteinIDs,beta_fun,mCurveFit,rmse,T_m,T_opt,T_H,fN_TH,r2adj] = fitTPPBetaFcn(denaturationDataFile, proteins, plotting)
%% [proteinIDs,beta_fun,mCurveFit,rmse,T_m,T_opt,T_H,fN_TH,r2adj] = fitTPPBetaFcn(denaturationDataFile, proteins, plotting)
% Reads the formatted denaturation data obtained from the Meltome Atlas
% (Jarzab et al. 2020, Nat. Methods.) and fits melting curves
% using the beta growth function (Yin et al. 2013, doi: 10.1093/aob/mcg029):
% The function was mirrored by the y-axis to model sigmoidal decay:
%
%         f(T) = a * (1 + (b + T) / (b - c)) * (-T/b) ^ (b / (b-c)).
%
% Input:
%       char denaturationDataFile:          path to file containing protein
%                                           stability data with columns:
%                                           Protein_ID;gene_name;meltPoint;
%                                           fold_change;temperature
%       cellstr proteins:                   (optional) array of protein IDs
%                                           for which melting curves should
%                                           be fitted
%                                           default: all protein IDs
%                                           contained in the provided file
%       logical plotting:                   (optional) if true, saves PDF
%                                           files containing diagnostic 
%                                           plots of at maximum 50 random
%                                           curve fits
%                                           default: false
%
% Output:
%       cellstr proteinIDs:                 array containing the unique set
%                                           of protein IDs found the in
%                                           given data file or only for the
%                                           protein IDs in the "proteins"
%                                           array
%       function_handle f_curve:            function which coefficients
%                                           were fitted to
%       double mCurveFit:                   n x 3 matrix containing the coefficients
%                                           estimated for the sigmoid fit,
%                                           where n is the number of
%                                           protein IDs, the coefficients
%                                           are a,b,c for the function
%                                           given above
%       double rmse:                        root mean square error of
%                                           melting curve fits
%       double T_m:                         estimated melting point for
%                                           each protein, which is the
%                                           temperature [°C] at which the
%                                           function above has a value of
%                                           0.5
%       double T_opt:                       estimated optimal temperature
%                                           with respect to protein
%                                           stability (function parameter b)
%       double T_H:                         heat denaturation temperature
%                                           (maximum of measured
%                                           temperature range)
%       double fN_TH:                       fraction of non-denatured protein
%                                           at T_H as predicted by the
%                                           function above
%       double r2adj:                       adjusted R2 value

if nargin > 0 && exist(denaturationDataFile, 'file')
    % read data from file
    denaturationTable = readtable(denaturationDataFile,...
        'ReadVariableNames',true);
elseif nargin == 0
    fprintf(['\nUSAGE:\n\n[proteinIDs,beta_fun,mCurveFit,rmse,T_m,T_opt,T_H,fN_TH]',...
        '= ...\n\tfitTPPBetaFcn(denaturationDataFile, proteins, plotting)\n\n']);
    return
else
    error('The given denaturation data file at %s does not exist.\n',...
        denaturationDataFile)
end

if nargin < 2 || isempty(proteins)
    % get a unique list of all proteins (first column)
    proteinIDs = unique(denaturationTable.(1), 'stable');
elseif ~iscellstr(proteins)
    error('The given protein ID array is not of type cellstr.')
else
    proteinIDs = unique(...
        proteins(ismember(proteins, strtok(denaturationTable.(1), '_'))),...
        'stable');
    denaturationTable = denaturationTable(...
        ismember(strtok(denaturationTable.(1), '_'), proteinIDs),:);
end

if nargin < 3
    plotting = false;
end

% initialize coefficients matrix and error
mCurveFit   = nan(numel(proteinIDs), 3);
rmse        = nan(numel(proteinIDs), 1);
r2adj       = nan(numel(proteinIDs), 1);
T_m         = nan(numel(proteinIDs), 1);
T_opt       = nan(numel(proteinIDs), 1);
T_H         = nan(numel(proteinIDs), 1);
% fraction of native protein at T_H according to fitted melting curve
fN_TH       = nan(numel(proteinIDs), 1);

% Beta growth function used to fit the melting curve
% F(1): a, F(2): b, F(3): c
beta_fun = @(F,T) F(1)*(1+(F(2)+T)./(F(2)-F(3))).*(-T./F(2)).^(F(2)/(F(2)-F(3)));
ft = fittype('a*(1+(b+T)./(b-c)).*(-T./b).^(b/(b-c))', 'independent', 'T');

% check if every protein has the same number of entries i.e. temperature
% steps
assumedTempSteps = size(denaturationTable, 1)/numel(proteinIDs);
tmp_ids_1 = strtok(denaturationTable.(1), '_');
tmp_ids_2 = strtok(proteinIDs, '_');
checkTempSteps = all(cellfun(@(x)...
    sum(ismember(tmp_ids_1, x)),...
    tmp_ids_2) == assumedTempSteps);
clear assumedTempSteps tmp_ids_1 tmp_ids_2

if ~checkTempSteps
    error('The given file does not contain the same number of measurements for each protein!')
else
    X = unique(denaturationTable.(5), 'stable');
    X(isnan(X)) = [];
    [X, order] = sort(X);
end
clear checkTempSteps assumedTempSteps

fprintf('\nFitting melting curves for %d proteins\n', numel(proteinIDs))
disp(repmat('_', 1, 40))

% initialize counter for the number of plots
c = 0;
for i = 1:numel(proteinIDs)
    
    if mod(i,500)==0
        fprintf('Done with %d proteins...\n', i)
    end
    
    % find row indices for the current protein
    rowIdx = ismember(strtok(denaturationTable.(1), '_'), strtok(proteinIDs(i), '_'));
    
    % find fold changes for each temperature and order them numerically by
    % associated temperature
    Y = denaturationTable.(4)(rowIdx)';
    Y = Y(order)';
    try
        [cf, gof, output] = fit(X, Y, ft, 'StartPoint', [1 -25 -50]);
    catch ME
        fprintf('Protein #%d: %s\n', i, ME.message)
        continue
    end
    
    mCurveFit(i, :) = [cf.a cf.b cf.c];
    rmse(i) = gof.rmse;
    r2adj(i) = gof.adjrsquare;
    
    if output.exitflag > 0
        curve_valid = true;
    else
        curve_valid = false;
    end
    
    if curve_valid
        T_opt(i) = -mCurveFit(i, 2);
        
        % find all points where the value of the fitted function reaches
        % half of the amplitude
        start_points = linspace(-mCurveFit(i, 2), 300, 5);
        roots = nan(numel(start_points), 1);
        options = optimset('Display', 'off');
        for j=1:numel(start_points)
            roots(j) = fzero(...
                @(x)beta_fun(mCurveFit(i, :), x)-0.5*mCurveFit(i, 1),...
                start_points(j),...
                options...
                );
            %             roots(j) = fzero(...
            %                 @(x)beta_fun(mCurveFit(i, :), x)-0.5,...
            %                 start_points(j),...
            %                 options...
            %                 );
        end
        
        roots_uniq = unique(round(roots, 4));
        roots_clean = roots_uniq(~isnan(roots_uniq));
        if numel(roots_clean) == 1
            fprintf('%s: Only one root found for f(x)-0.5*f(x): %.2g\n',...
                strtok(proteinIDs{i},'_'), roots_clean)
        elseif numel(roots_clean) == 2
            T_m(i) = roots_clean(2);
        else
            fprintf('%s: Roots for f(x)-0.5*f(x) could not be determined.\n',...
                strtok(proteinIDs{i}, '_'))
        end
        
        if isnan(T_m(i)) || T_m(i) <= T_opt(i) || T_m(i) < 25
            T_opt(i) = NaN;
            T_m(i) = NaN;
        else
            T_H(i) = max(X);
            fN_TH(i) = beta_fun(mCurveFit(i, :),T_H(i));
        end
        
        if plotting && rand <= 0.5 && c <= 50
            c = c + 1;
            fig = figure('Visible', 'off');
            
            % temperature range [°C]
            tempRange = 20:70;
            
            % plot original data
            plot(X, Y, 'k.', 'MarkerSize', 15)
            hold on
            % plot fitted sigmoid function
            plot(tempRange,...
                beta_fun(mCurveFit(i, :), tempRange),...
                'k',...
                'LineWidth', 2)
            
            % plot vertial line for optimal temperature
            line(...
                [T_opt(i) T_opt(i)],...
                [0 1.5],...
                'Color', 'b',...
                'LineWidth', 1,...
                'LineStyle', '--')
            
            % plot vertial line for melting temperature
            line(...
                [T_m(i) T_m(i)],...
                [0 1.5],...
                'Color', 'k',...
                'LineWidth', 1,...
                'LineStyle', '--')
            hold off
            
            % configure plot
            ylim([0 1.5])
            xlabel('T [°C]',...
                'FontSize',14)
            ylabel('Relative soluble fraction',...
                'FontSize',14)
            text(.75, .7,...
                sprintf('RMSE = %.4g', rmse(i)),...
                'Units', 'normalized')
            text(.75, .65,...
                sprintf('R^2_{adj} = %.4g', gof.adjrsquare),...
                'Units', 'normalized')
            legend({'data', 'fit', 'T_{opt}', 'T_m'},...
                'Box','off',...
                'Location','northeast')
            
            title(strrep(proteinIDs{i}, '_', '\_'))
            
            if mCurveFit(i, 1) > 3 && ~isnan(T_m(i))
                disp('High amplitude')
            end
            exportgraphics(fig,...
                fullfile('protein-stability', 'curve_fits',  [proteinIDs{i} '.pdf']),...
                'Resolution', 100)
        end
    end
end

end
