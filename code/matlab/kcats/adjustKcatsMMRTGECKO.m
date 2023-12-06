function model_adj = adjustKcatsMMRTGECKO(model, T, T0, kcat_scaling)
%% model_adj = adjustKcatsMMRTGECKO(model, T, T0, kcat_scaling)
% Adjust turnover numbers to a given temperature using the MMRT model
% described by Hobbs et al. (2013)
% 
% By default, a reference temperature of 20 °C is used.
% 
% The temperature adjustment of kcats does not rely on the kcats
% contained in the model (only uses MMRT fit). 
% HOWEVER, if a scaling factor is used, the kcats of enzymes for which no 
% MMRT fit is available will be altered, based on the kcat entries in the 
% stoichiometric matrix.
%
% Input:
%       struct model:       GECKO-formatted metabolic model with additional fields:
%                           DH, DS, DCp (obtained by running fitMMRTGECKOModel)
%       scalar T:           temperature [K]
%       scalar T0:          (opt) reference temperature [K] (default: 20 °C)
%       double kcat_scaling:(opt) scaling factor for multiplication of 
%                           adjusted kcat values (if not given, will be
%                           read from config file)
% Output:
%       struct model_adj:   temperature-adjusted model
% Philipp Wendering, University of Potsdam (philipp.wendering@gmail.com) 24/06/2023

% constants
NPROT = numel(model.enzymes);
ENZ_PFX = 'prot_';

% set reference temperature
if nargin < 3 || isempty(T0); T0 = celsius2kelvin(20); end
% set kcat scaling factor
if nargin < 4 || isempty(kcat_scaling)
    kcat_scaling = config('kcat_scaling');
    warning('Multiplying kcats with %d.', kcat_scaling)
end

% check some fields to find out if we are dealing with a TFA model (matTFA)
if any(startsWith(model.rxns, 'DGo_')) && any(startsWith(model.rxns, 'NF_'))
    fprintf('Detected TFA model, temporarily updating model fields.\n')
    % temporarily change update fields
    mets = model.mets;
    % temporarily remove 'M_' before protein metabolites
    model.mets = regexprep(model.mets, '^M_prot_', 'prot_');
    
    tfa_flag = 1;
else
    tfa_flag = 0;
end
    
% find row indices of enzymes
enzRowIdx = find(startsWith(model.mets,ENZ_PFX));
enzRowIdx = setdiff(enzRowIdx,find(ismember(model.mets,'prot_pool')));

if NPROT ~= numel(enzRowIdx)
    error('Numbers of proteins in S and in ''enzymes'' field are not equal')
end

model_adj = model;

% initilize counter for low kcats (defined by t)
t = 1e-5;
counter = 0;
for i=1:NPROT
    
    % find reactions associated with the current protein and are not draw
    % reactions (not a problem in a normal GECKO EC model but the matTFA
    % toolbox creates forward and reverse reactions for every reaction, so
    % the reverse draw reactions have negative coefficients)
    protRxnsIdx = find(model.S(enzRowIdx(i),:)<0 & ~contains(model.rxns, 'draw_prot_')');
    
    for j=1:numel(protRxnsIdx)
        rxn_idx = protRxnsIdx(j);
        if model.DH(i,rxn_idx)~=0 && model.DS(i,rxn_idx)~=0 && model.DCp(i,rxn_idx)~=0
            ln_k_adj = mmrtFun(T,model.DH(i,rxn_idx),model.DS(i,rxn_idx),model.DCp(i,rxn_idx),T0);
            k_adj = kcat_scaling*exp(ln_k_adj); % adjusted kcat [s^-1]
            if k_adj<=t;counter=counter+1;end % increase counter
            model_adj.S(enzRowIdx(i),protRxnsIdx(j)) = -1./3600/k_adj; % assing adjusted kcat
        else
            kcat = -1/model_adj.S(enzRowIdx(i),protRxnsIdx(j))/3600;
            k_scaled = kcat_scaling*kcat;
            model_adj.S(enzRowIdx(i),protRxnsIdx(j)) = -1./3600/k_scaled;
        end
    end
end

if tfa_flag
    % if TFA model, revert model fields
    model_adj.mets = mets;
end

if counter > 0
    fprintf('Number of low kcats after adjustment (<= %.2g):\t%d\n',t,counter);
end
end