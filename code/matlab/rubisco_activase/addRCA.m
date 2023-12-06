function new_model = addRCA(model)
% add RubisCO activase (RCA) mechanism to GECKO model
%
% RuBisCO(inactive) + n RCA + ATP ---> RuBisCO(active) + ADP + Pi
%
% This reaction includes both the re-activation mechanism of RuBisCO after
% RuBP binding as well as the re-activation cost for RCA (ATP -> ADP + Pi).
% 
% The coeffiecient n is determined by the scaling factor used for kcat 
% values (see config.m).
% 
% Further, the original draw reaction for the RuBisCO large subunit is 
% replaced with the inactive RuBisCO form. 
% A draw reaction for the RCA protein is also added.


% add metabolites for inactive RuBisCO and RuBisCO activase
new_model = addMultipleMetabolites(model,...
    {'RBC_inact', 'prot_P10896'});

% change metabolite in RuBisCO draw reaction to inactive protein
new_model = changeRxnMets(new_model, 'prot_O03042', 'RBC_inact', 'draw_prot_O03042');

% add draw reaction for RCA
new_model = addReaction(new_model,...
    'draw_prot_P10896',...
    'reactionName','draw_prot_P10896',...
    'reactionFormula','51.981 prot_pool -> prot_P10896',...
    'reversible', false,...
    'upperBound',Inf ...
    );

% add RCA reaction
new_model = addReaction(new_model, 'RCA_h',...
    'reactionName', 'RuBisCO activation',...
    'reactionFormula', [...
    'RBC_inact + ' num2str(config('kcat_scaling')) ' prot_P10896 + ATP_h + H2O_h + CO2_h -> ',...
    'prot_O03042 + ADP_h + Pi_h']);

% add melting curve fit parameters for RCA to model
[~, ~, ~, ~, T_m, T_opt, T_H, fN_TH] = ...
    fitTPPBetaFcn(config('meltomeFile'), {'P10896'});
new_model.T_opt(end+1) = T_opt;
new_model.T_m(end+1) = T_m;
new_model.T_H(end+1) = T_H;
new_model.meltingCurves.fN_TH(end+1) = fN_TH;
new_model.enzymes{end+1} = 'P10896';
new_model.enzNames{end+1} = 'RCA';

% obtain MMRT fit for RCA
[new_model.DH, new_model.DS, new_model.DCp] = fitMMRTGECKOModel(new_model);


end