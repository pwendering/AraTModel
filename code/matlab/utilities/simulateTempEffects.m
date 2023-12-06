function [simulationResult, photParams, fvcbParams, gurobiProblem] = simulateTempEffects(varargin)
%% [simulationResult, photParams, fvcbParams, gurobiProblem] = simulateTempEffects(varargin)
% Simluation of temperature effects on the distribution of fluxes through a
% metabolic network considering enzyme kinetics, protein stability, and CO2
% uptake / photosynthesis constraints.
% INPUT
%   struct TGEM:                metabolic model with additional fields
%                               added using the 'createTGEM' function
%   double I:                   irradiance (umol/m2/s)
%   double tempRange:           (optional) temperatures in Kelvin at which 
%                               the flux
%                               default: 10 temperature steps between 10 °C
%                               and 40 °C
%                               distribution should be predicted
%   logical adjPtotFlag:        (optional) whether to adjust the total 
%                               protein content to different temperatures
%                               default: true
%   logical adjKcatFlag:        (optional) whether to adjust the kcat
%                               values to different temperatures
%                               default: true
%   double nRxns:               (optional) number of reaction flux variables in the
%                               problem; it is assumed that the reaction
%                               indices go from 1:nRxns
%                               default: column dimension of S
%   logical saBool:             (optional) specifies whether or not a
%                               sensitivity analysis is performed
%                               default: false
%   double saPercentage:        (optional) if saBool is true, specifies the
%                               percentage for increase of the parameter
%                               given in saParameter at the given
%                               temperature (between 0 and 1)
%   char saParameter:           (optional) name of the parameter that
%                               should be increased by saPercentage
%   logical visualize:          (optional) boolean value that indicates
%                               whether the results should be visualized in
%                               plots for
%                               (1) predicted relative growth rate, net CO2
%                               assimilation rate predicted by the
%                               metabolic model and the FvCB model
%                               (2) temperature dependences of parameters
%                               that enter the CO2 uptake and
%                               photosynthesis constraints as well as the
%                               FvCB model prediciton
%                               default: false
%   struct solverParams:        Parameters that should be passed to the
%                               Gurobi solver
%                               default: []
% OUTPUT
%   struct simulationResult:    results from predicted flux distribution
%                               .v      matrix with flux distributions as
%                                       column vectors (#reactions x #temperatures)
%                               .mu     average predicted relative growth rate [h^-1]
%                               .mu_min minimum predicted relative growth rate [h^-1]
%                                       (obtained from pool solutions)
%                               .mu_max maximum predicted relative growth rate [h^-1]
%                                       (obtained from pool solutions)
%                               .A      average predicted net CO2 assimilation rate
%                                       (A = v_c - 0.5v_o - Rd; Rd ~ export
%                                       of CO2 from mitochondrion) [umol m^-2 s^-1]
%                               .A_min  minimum predicted net CO2 assimilation rate
%                                       (obtained from pool solutions) [umol m^-2 s^-1]
%                               .A_max  maximum predicted net CO2 assimilation rate
%                                       (obtained from pool solutions) [umol m^-2 s^-1]
%                               .A_net  net CO2 assimilation rate
%                                       calculated by
%                                       sum(fluxes cons. CO2) - sum(fluxes prod. CO2)
%                               .objVal average objective value; if the
%                                       biomass reaction is the objective,
%                                       this field contains the same value
%                                       as the "mu" field [mmol gDW^-1 h^-1]
%                               .objVal_min
%                                       minimum objective value [mmol gDW^-1 h^-1]
%                               .objVal_max
%                                       maximum objective value [mmol gDW^-1 h^-1]
%   struct photParams:          temperature-adjusted photosynthesis/CO2
%                               uptake parameters
%                               .Kc     Michaelis constant for the RuBisCO
%                                       carboxylation reaction [ubar]
%                               .Ko     Michaelis constant for the RuBisCO
%                                       oxygenation reaction [ubar]
%                               .V_c_max
%                                       maximum velocity of the
%                                       carbobylation reaction [umol m^-2 s^-1]
%                               .phi    ratio of oxygenation to
%                                       carboxylation flux
%                               .S_co   specificity of carboxylation over
%                                       oxygenation [bar bar^-1]
%                               .g_m    mesophyll conductance [mol bar^-1 m^-2 s^-1]
%                               .g_s    stomatal conductance [mol m^-2 s^-1]
%                               .LMA    leaf mass per area [g m^-2]
%   struct fvcbParams:          (only non-empty if visualization is true)
%                               Parameters of the FvCB model
%                               .k_c    turnover number (kcat) of the RuBisCO
%                                       carboxylation reaction [s^-1]
%                               .k_o    turnover number (kcat) of the RuBisCO
%                                       oxygenation reaction [s^-1]
%                               .K_c    Michaelis constant for the RuBisCO
%                                       carboxylation reaction [ubar]
%                               .K_o    Michaelis constant for the RuBisCO
%                                       oxygenation reaction [ubar]
%                               .E_a    activation energies for
%                                       K_c (.K_c)
%                                       K_o (.K_o)
%                                       V_c_max (.V_c_max) (not required
%                                           for FvCB model)
%   struct gurobiProblem:       optimization problem that was solved to
%                               obtain the LAST solution in gurobiSolution

p = parseInput(varargin);

TGEM = p.Results.TGEM;
I = p.Results.I;
tempRange = p.Results.tempRange;
adjPtotFlag = p.Results.adjPtotFlag;
adjKcatFlag = p.Results.adjKcatFlag;
nRxns = p.Results.nRxns;
saBool = p.Results.saBool;
saPercentage = p.Results.saPercentage;
saParameter = p.Results.saParameter;
visualize = p.Results.visualize;
solverParams = p.Results.solverParams;

% add sink reaction for water (for modelling transpiration)
TGEM.TRANS_ID = 'Sink_H2O';
if ismember('H2O[e]', TGEM.mets)
    water_id = 'H2O[e]';
elseif ismember('H2O[c]', TGEM.mets)
    water_id = 'H2O[c]';
elseif ismember('H2O_e', TGEM.mets)
    % for GECKO-formatted model
    water_id = 'H2O_e';
elseif ismember('H2O_c', TGEM.mets)
    % for GECKO-formatted model
    water_id = 'H2O_c';
else
    water_idx = find(ismember(TGEM.metFormulas, 'H2O'));
    water_idx_ext = ~cellfun(@isempty, regexp(TGEM.mets(water_idx), '[\[_]e\]?$'));
    water_idx_cyt = ~cellfun(@isempty, regexp(TGEM.mets(water_idx), '[\[_]c\]?$'));
    if any(water_idx_ext)
        water_id = TGEM.mets{water_idx(water_idx_ext)};
    elseif any(water_idx_cyt)
        water_id = TGEM.mets{water_idx(water_idx_cyt)};
    else
        error('Could not find water metabolite ID.')
    end
end

% attempt to find an existing sink reaction for water
water_sink_idx = TGEM.S(findMetIDs(TGEM, water_id),:)<0 & sum(TGEM.S~=0)==1;
add_water_sink_bool = false;
if sum(water_sink_idx) == 1
    fprintf('Using the following reaction for transpiration modelling: %s\n',...
        TGEM.rxns{water_sink_idx})
    TGEM.TRANS_ID = TGEM.rxns(water_sink_idx);
elseif sum(water_sink_idx) > 1
    fprintf('Multiple sink reactions found for water:\n')
    disp(TGEM.rxns(water_sink_idx))
    water_sink_idx = find(water_sink_idx);
    fprintf('Choosing reaction the following reaction for transpiration modelling: %s\n',...
        TGEM.rxns{water_sink_idx(1)})
    TGEM.TRANS_ID = TGEM.rxns(water_sink_idx(1));
    fprintf('Blocking remaining water sink reactions\n')
    TGEM.ub(water_sink_idx(2:end)) = 0;
else
    fprintf('No water sink reaction was found, adding sink reaction\n')
    TGEM = addReaction(TGEM, TGEM.TRANS_ID,...
        'reactionFormula', [water_id ' ->']);
    add_water_sink_bool = true;
end

% in case the H2O sink reaction is blocked, change upper bound to 1000
TGEM.ub(findRxnIDs(TGEM, TGEM.TRANS_ID)) = 1000;

% block water sink reaction in another compartments
water_idx = find(ismember(TGEM.metFormulas, 'H2O'));
water_sink_idx = any(TGEM.S(water_idx,:)<0) & sum(TGEM.S~=0)==1 & ~ismember(TGEM.rxns, TGEM.TRANS_ID)';
fprintf('Blocking remaining water sink reactions:\n')
disp(TGEM.rxns(water_sink_idx))
TGEM.ub(water_sink_idx) = 0;

if isfield(TGEM, 'vartype')
    % only TFA
    TGEM.vartype(end) = 'C';
end

if ~exist('nRxns','var') || isempty(nRxns)
    fprintf('Number of reaction variables not given, using column dimension of S.\n')
    nRxns = size(TGEM.S, 2);
elseif add_water_sink_bool
    nRxns = nRxns + 1;
end

% convert stoichiometric matrix from sparse to full
TGEM.S = full(TGEM.S);

% increase the upper bounds of CO2 and photon import reactions
TGEM.ub(findRxnIDs(TGEM, TGEM.CO2_IMP_ID)) = 1e4;
TGEM.ub(findRxnIDs(TGEM, TGEM.ABS_ID)) = 1e4;

%% Simulation
fprintf('SELECTED OBJECTIVE: %s\n',TGEM.rxnNames{TGEM.c==1})
% initialize result structures for plotting
[simulationResult,photParams] = initResultStruct(tempRange, nRxns);

for i=1:numel(tempRange)
    
    fprintf('\nCURRENT TEMPERATURE: %3.2f K (%.2f °C)\n\n',...
        tempRange(i),kelvin2celsius(tempRange(i)))
    
    tmpModel = TGEM;

    % limit protein pool if present
    if ismember('prot_pool_exchange', tmpModel.rxns)
        if adjPtotFlag
            P = tmpModel.proteinModel(kelvin2celsius(tempRange(i)));
        else
            P = 1e4;
        end
        fprintf('Total protein content: %.2f g/gDW\n', P)
        tmpModel.ub(findRxnIDs(tmpModel, 'prot_pool_exchange')) = P;
    end
    
    if adjKcatFlag
        disp('// Kcat adjustment using MMRT')
        kcat_scaling = config('kcat_scaling');
        kcatAdjModel = adjustKcatsMMRTGECKO(tmpModel, tempRange(i), [],...
            kcat_scaling);
    else
        kcatAdjModel = tmpModel;
    end
    clear tmpModel
    
    disp('// Growth simulation with FvCB constraints')
    [solution,adjParams,~,gurobiProblem] = modelFvCBConstraints(...
        kcatAdjModel,...
        tempRange(i),...
        I,...
        solverParams,...
        nRxns,...
        saBool,...
        saPercentage,...
        saParameter);
    
    if isequal(solution.status,'OPTIMAL')
        fprintf('v_c:\t%.6g mmol/gDW/h\n', solution.x_step2(findRxnIDs(kcatAdjModel, kcatAdjModel.C_ID)))
        fprintf('v_o:\t%.6g mmol/gDW/h\n', solution.x_step2(findRxnIDs(kcatAdjModel, kcatAdjModel.O_ID)))
        fprintf('CO2 (m->c):\t%.6g mmol/gDW/h\n', solution.x_step2(findRxnIDs(kcatAdjModel, kcatAdjModel.RESP_ID)))
    end
    % fill in arrays for adjusted parameters
    fields = setdiff(fieldnames(solution),{'status'});
    for j=1:numel(fields)
        if ~isempty(solution.(fields{j}))
            simulationResult.(fields{j})(:,i) = solution.(fields{j});
        end
    end
    
    % fill in arrays for adjusted parameters
    fields = fieldnames(adjParams);
    for j=1:numel(fields)
        photParams.(fields{j})(i) = adjParams.(fields{j});
    end
    fprintf('\n---\n')
end

%% Visualization
if visualize
    close all
    fvcbParams = getFvCBParams();
    plotPhotSimResults(tempRange,simulationResult,photParams,fvcbParams,I);
else
    fvcbParams = struct;
end

    function [predFig, adjPhotParamFig] = plotPhotSimResults(tempRange,solution,photParams,fvcbParams,I)
        %% [predFig, adjPhotParamFig] = plotPhotSimResults(tempRange,solution,photParams,I)
        % Plot results from simulations of net CO2 assimilation rate and growth
        % with a metabolic model augmented with constraints on CO2 uptake.
        % Input:
        %       double tempRange:           temperature interval for which the
        %                                   simulations were carried out
        %       struct solution:            contains minimum, maximum and average values
        %                                   for growth and net CO2 assimilation rate;
        %                                   all array fields must have the same
        %                                   dimensions as tempRange
        %       struct photParams:          contains temperature-adjusted
        %                                   photosynthesis parameters;
        %                                   all array fields must have the same
        %                                   dimensions as tempRange
        %       double fvcbParams:          parameters for the FvCB model
        %       double I:                   value for the irradiance [umol m^-2 s^-1]
        % Output:
        %       Figure predFig:             figure handle for simulation results
        %       Figure adjPhotParamFig:     figure handle for adjusted parameters
        
        %% ~~~~~~~~~~~~~~~~~~~~~ Plot simulation results ~~~~~~~~~~~~~~~~~~~~~~~~ %
        predFig = figure;
        
        % ~~~~~~~ Predicted growth rate ~~~~~~~ %
        subplot(1,3,1)
        plot(kelvin2celsius(tempRange), solution.mu, 'LineWidth', 1.5)
        hold on
        % minimum and maximum of pool solutions
        plot(kelvin2celsius(tempRange), solution.mu_min, 'LineWidth', 1.5,...
            'Color', [1 .5 .5], 'LineStyle', '--')
        plot(kelvin2celsius(tempRange), solution.mu_max, 'LineWidth', 1.5,...
            'Color', [1 .5 .5], 'LineStyle', '--')
        % add optimal temperature
        growthOptTemp = tempRange(solution.mu==max(solution.mu));
        yLim = get(gca,'YLim');
        ylim([0 max(yLim)])
        text(0.3,0.1,['T_{opt}=' num2str(kelvin2celsius(growthOptTemp), '%.1f') ' °C'],...
            'Units', 'normalized');
        % add panel letter
        text(.05,.95,'A', 'FontSize', 14, 'FontWeight', 'bold', 'Units', 'normalized')
        % format axes
        xlabel('Temperature [°C]')
        ylabel('predicted growth [h{^1}]')
        xlim([min(kelvin2celsius(tempRange)) max(kelvin2celsius(tempRange))])
        box on
        set(gca,'LineWidth', 1.5, 'FontSize', 10)
        
        % ~~~~~~~ Predicted net CO2 assimilation rate (metabolic model) ~~~~~~~ %
        subplot(1,3,2)
        plot(kelvin2celsius(tempRange), solution.A, 'LineWidth', 1.5)
        hold on
        % plot(kelvin2celsius(tempRange), solution.A_net, 'LineWidth', 1.5)
        % hold on
        % minimum and maximum of pool solutions
        plot(kelvin2celsius(tempRange), solution.A_min, 'LineWidth', 1.5,...
            'Color', [1 .5 .5], 'LineStyle', '--')
        plot(kelvin2celsius(tempRange), solution.A_max, 'LineWidth', 1.5,...
            'Color', [1 .5 .5], 'LineStyle', '--')
        % add optimal temperature
        AOptTemp = tempRange(solution.A==max(solution.A));
        % AOptTemp = tempRange(solution.A_net==max(solution.A_net));
        yLim = get(gca,'YLim');
        ymax = max(yLim);
        text(0.3,0.1,['T_{opt}=' num2str(kelvin2celsius(AOptTemp), '%.1f') ' °C'],...
            'units','normalized');
        % add panel letter
        text(.05,.95,'B', 'FontSize', 14, 'FontWeight', 'bold', 'Units', 'normalized')
        
        % format axes
        xlabel('Temperature [°C]')
        ylabel('Net CO_2 assimilation rate [\mumol m^{-2} s^{-1}]')
        xlim([min(kelvin2celsius(tempRange)) max(kelvin2celsius(tempRange))])
        box on
        set(gca,'LineWidth', 1.5, 'FontSize', 10)
        
        % ~~~~~~~ Predicted net CO2 assimilation rate (FvcB model) ~~~~~~~ %
        % CO2 partial pressues
        pCO2 = [55 110 165 230 330 480];
        % initialize assimulation rate matrix (T x p(CO2))
        A_far = zeros(numel(tempRange),numel(pCO2));
        
        for i=1:numel(tempRange)
            for j=1:numel(pCO2)
                % parameterized for A. thaliana Col-0
                A_far(i,j) = farquhar('T', tempRange(i),'C', pCO2(j),'I', I,...
                    'kc',fvcbParams.k_c, 'ko', fvcbParams.k_o,'K_c',fvcbParams.K_c,...
                    'K_o',fvcbParams.K_o, 'E_a',fvcbParams.E_a);
            end
        end
        
        subplot(1,3,3)
        plot(kelvin2celsius(tempRange),A_far,'Color', lines(1),'LineWidth',1.5)
        
        
        % add annotation for CO2 partial pressures
        C_units = repmat({''},1,numel(pCO2));
        C_units(1) = {' \mubar'};
        C_prefix = repmat({''},1,numel(pCO2));
        C_prefix(1) = {'C='};
        text([30 44 42 41 39 37],max(A_far)+[.7 .1 .1 .1 .1 .1],...
            arrayfun(@(i)[C_prefix{i}, num2str(pCO2(i)), C_units{i}],1:numel(pCO2),'un',0),...
            'FontSize',8)
        % add annotation for irradiance
        text(13,-7,['I = ', num2str(I), ' \mumol m^{-2} s^{-1}'],'FontSize',8)
        % add panel letter
        text(.05,.95,'C', 'FontSize', 14, 'FontWeight', 'bold', 'Units', 'normalized')
        % format axes
        xlim([min(kelvin2celsius(tempRange)) max(kelvin2celsius(tempRange))])
        ylabel('Net CO_2 assimilation rate [\mumol m^{-2} s^{-1}]')
        xlabel('Temperature [°C]')
        box on
        set(gca,'LineWidth', 1.5, 'FontSize',10)
        
        % save figure as PNG image
        set(gcf, 'Position', [100 200 800 400])
        print(['ath_growth_A_obj_1_max_bio_2_min_v_abs_I_' num2str(I) '.png'],...
            '-painters', '-dpng')
        
        %% ~~~~~~~~~~~~ Plot T-adjusted photosynthesis parameters ~~~~~~~~~~~~~~~ %
        adjPhotParamFig = figure;
        % Y-axes labels
        yLabels = {
            'K_c (\mubar)', 'K_o (\mubar)', 'k_c (s^{-1})', 'k_o (s^{-1})',...
            'V_{cmax} (\mumol m^{-2} s^{-1})', '\phi', 'S_{c/o} (bar bar^{-1})',...
            'g_m (mol bar^{-1} m^{-2} s^{-1})', 'g_s (mol bar^{-1} m^{-2} s^{-1})',...
            'J_{max} (\mumol m^{-2} s^{-1})', 'LMA (g m^{-2})'...
            };
        fNames = fieldnames(photParams);
        % define plot dimensions (eight panels) and panel letters
        plotXDim = floor(sqrt(numel(yLabels)));
        plotYDim = ceil(numel(yLabels)/plotXDim);
        plotDim = [plotXDim, plotYDim];
        letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
        
        % plot panels with adjsted parameters of the considered temperature interval
        c=0;
        for i=1:numel(fNames)
            c=c+1;
            % define position and add plot
            subplot(plotDim(1),plotDim(2),c)
            plot(kelvin2celsius(tempRange),photParams.(fNames{i}),'k', 'LineWidth',1.5)
            % add letter
            text(.05, .9, letters(c),...
                'FontSize', 12,...
                'FontWeight', 'bold',...
                'FontName', 'Arial',...
                'Units', 'normalized')
            % format axes
            ylabel(yLabels{i})
            set(gca,'LineWidth', 1.2, 'FontSize', 10)
            
            % set x ticks
            xticks(kelvin2celsius(min(tempRange):10:max(tempRange)))
            
            box on
        end
        % replace individual X-axis labels by one common label
        axHandle = axes(adjPhotParamFig, 'visible', 'off');
        axHandle.XLabel.Visible='on';
        xlabel(axHandle, 'Temperature (°C)', 'FontSize', 12)

        % save figure as PNG
        set(gcf, 'Position', [120 51 1000 550])
        print('ath_phot_params_T_dependence.png', '-painters', '-dpng')
    end
    function [simulationResult,photParams] = initResultStruct(tempRange,nRxns)
        % initialize result structures for plotting
        photParams = struct;
        photParams.Kc = nan(size(tempRange));
        photParams.Ko = nan(size(tempRange));
        photParams.kc = nan(size(tempRange));
        photParams.ko = nan(size(tempRange));
        photParams.V_c_max = nan(size(tempRange));
        photParams.phi = nan(size(tempRange));
        photParams.S_co = nan(size(tempRange));
        photParams.g_m = nan(size(tempRange));
        photParams.g_s = nan(size(tempRange));
        photParams.J_max = nan(size(tempRange));
        photParams.LMA = nan(size(tempRange));
        
        simulationResult = struct;
        simulationResult.x_step1 = nan(nRxns+2,numel(tempRange));
        simulationResult.x_step2 = nan(nRxns+2,numel(tempRange));
        simulationResult.mu = nan(size(tempRange));
        simulationResult.mu_min = nan(size(tempRange));
        simulationResult.mu_max = nan(size(tempRange));
        simulationResult.A = nan(size(tempRange));
        simulationResult.A_min = nan(size(tempRange));
        simulationResult.A_max = nan(size(tempRange));
        simulationResult.A_net = nan(size(tempRange));
        simulationResult.objVal = nan(size(tempRange));
        simulationResult.objVal_min = nan(size(tempRange));
        simulationResult.objVal_max = nan(size(tempRange));
    end

    function p = parseInput(args)
        
        DEFAULT_T_RANGE = celsius2kelvin(linspace(10, 40, 10));
        
        validTemp = @(T) isnumeric(T) && all(T >= 273.15) && all(T <= 373.15);
        validNrRxns = @(n) n <= numel(args{1}.rxns);
        validSaPerc = @(p) 0 < p && p <= 1;
        
        p = inputParser;
        
        addRequired(p, 'TGEM', @validateTGEM);
        addRequired(p, 'I', @isnumeric)
        addOptional(p, 'tempRange', DEFAULT_T_RANGE, validTemp);
        addOptional(p, 'adjPtotFlag', true, @islogical);
        addOptional(p, 'adjKcatFlag', true, @islogical);
        addOptional(p, 'nRxns', [], validNrRxns);
        addOptional(p, 'saBool', false, @islogical);
        addOptional(p, 'saPercentage', 0, validSaPerc);
        addOptional(p, 'saParameter', '', @ischar);
        addOptional(p, 'visualize', false, @islogical);
        addOptional(p, 'solverParams', [], @isstruct);
        
        parse(p, args{:});
        
        if p.Results.saBool && p.Results.visualize
            fprintf('Visualization is disables with sensitivity analysis\n')
            p.Results.visualize = false;
        end
        
        if p.Results.saBool && (p.Results.saPercentage==0 || isempty(p.Results.saParameter))
            error('Parameter or percentage missing for sensitivity analysis')
        end
    end


end