function TGEM = createTGEM(model,is_ec)
%% TGEM = createTGEM(model)
% Add fields that allow for temperature adaptation.
% Input
%   struct model:           metabolic model with all irreversible reactions
%   logical is_ec:          true if model is a GECKO-formatted ecModel, 0
%                           otherwise; default: true
% Output
%   struct TGEM:            metabolic model with additional field required
%                           for temperature adjustment
% Philipp Wendering, University of Potsdam (philipp.wendering@gmail.com) 31/05/2022
if nargin < 2 || ~islogical(is_ec)
    is_ec = true;
end

% Configure MATLAB path and COBRA toolbox
fprintf('Configuring MATLAB path and COBRA solver...\n')
config
fprintf('Done\n')

%% Adding field with biomass reaction ID
model.BIO_ID = config('BIO_ID');

%% Protein melting curves
fprintf('Fitting melting curves...\n')
meltingCurves = struct;
meltomeFile = config('meltomeFile');

% set threshold of 0.6 for R2_adj (melting curve fits)
r2_t = 0.6;
    
if is_ec
    % get protein UniProt IDs
    ec_model_prot_ids = erase(...
        model.mets(...
        startsWith(model.mets,'prot_') &...
        ~contains(model.mets,'pool')),...
        'prot_');
    
    % fit melting curves
    [proteinIDs, meltingCurves.f, estimates, rmse, T_m, T_opt, T_H, fN_TH, r2adj]= ...
        fitTPPBetaFcn(meltomeFile, ec_model_prot_ids, false);
    
    T_m(r2adj<r2_t) = NaN;
    T_opt(r2adj<r2_t) = NaN;
    T_H(r2adj<r2_t) = NaN;
    fN_TH(r2adj<r2_t) = NaN;
    estimates(r2adj<r2_t, :) = NaN;
    
    tmp_meltome_ids = strtok(proteinIDs, '_');
    matchIdxMeltome = cellfun(@(x)find(ismember(tmp_meltome_ids, x)),...
        ec_model_prot_ids, 'un', 0);
    matchIdxModel = ~cellfun(@isempty, matchIdxMeltome);
    matchIdxMeltome = cell2mat(matchIdxMeltome);

    
    meltingCurves.fN_TH = nan(size(ec_model_prot_ids));
    meltingCurves.fN_TH(matchIdxModel) = fN_TH(matchIdxMeltome);
    
    meltingCurves.rmse = nan(size(ec_model_prot_ids));
    meltingCurves.rmse(matchIdxModel) = rmse(matchIdxMeltome);
    
    meltingCurves.estimates = nan(numel(ec_model_prot_ids),size(estimates,2));
    meltingCurves.estimates(matchIdxModel,:) = estimates(matchIdxMeltome,:);
    
    model.meltingCurves = meltingCurves;
    
    model.T_opt = nan(size(ec_model_prot_ids));
    model.T_opt(matchIdxModel) = T_opt(matchIdxMeltome);
    
    model.T_m = nan(size(ec_model_prot_ids));
    model.T_m(matchIdxModel) = T_m(matchIdxMeltome);
    
    model.T_H = nan(size(ec_model_prot_ids));
    model.T_H(matchIdxModel) = T_H(matchIdxMeltome);
    
else

    uniprotFile = config('uniprotFile');
    
    % match UniProt IDs with gene IDs from model
    upTable = readtable(uniprotFile,'ReadVariableNames', true,...
        'FileType', 'text');
    
    % find UniProt IDs for model gene IDs
    uniProtModel = cellfun(...
        @(x)upTable.Entry(contains(upTable.GeneNames, x, 'IgnoreCase', 1)),...
        model.genes, 'un', 0);
    model.geneUniprotID = cellfun(@char, uniProtModel, 'un', 0);
    
    % fit melting curves
    [proteinIDs, meltingCurves.f, estimates, rmse, T_m, T_opt, T_H, fN_TH, r2adj]= ...
        fitTPPBetaFcn(meltomeFile, model.geneUniprotID, false);
    proteinIDs = strtok(proteinIDs, '_');
    
    T_m(r2adj<r2_t) = NaN;
    T_opt(r2adj<r2_t) = NaN;
    T_H(r2adj<r2_t) = NaN;
    fN_TH(r2adj<r2_t) = NaN;
    estimates(r2adj<r2_t, :) = NaN;
    
    % assign melting curves to matching IDs
    matchIdxMeltome = cellfun(@(x)find(ismember(proteinIDs,x)),model.geneUniprotID,'un',0);
    matchIdxModel = ~cellfun(@isempty,matchIdxMeltome);
    matchIdxMeltome = cell2mat(matchIdxMeltome);
    
    meltingCurves.estimates = nan(numel(model.genes),3);
    meltingCurves.estimates(matchIdxModel,:) = estimates(matchIdxMeltome,:);
    
    meltingCurves.rmse = nan(size(model.genes));
    meltingCurves.rmse(matchIdxModel) = rmse(matchIdxMeltome);
    
    meltingCurves.fN_TH = nan(size(model.genes));
    meltingCurves.fN_TH(matchIdxModel) = fN_TH(matchIdxMeltome);
    
    model.meltingCurves = meltingCurves;
    
    model.T_opt = nan(size(model.genes));
    model.T_opt(matchIdxModel) = T_opt(matchIdxMeltome);
    
    model.T_m = nan(size(model.genes));
    model.T_m(matchIdxModel) = T_m(matchIdxMeltome);
    
    model.T_H = nan(size(model.genes));
    model.T_H(matchIdxModel) = T_H(matchIdxMeltome);
    
end

fprintf('\tMelting curves available for %d proteins in the model (%.4g%%)\n',...
    sum(matchIdxModel), 100*sum(matchIdxModel)/numel(matchIdxModel))

clear meltingCurves matchIdxMeltome matchIdxModel uniProtModel uniprotFile ...
    proteinIDs estimates rmse T_m T_opt T_H fN_TH meltomeFile tmp_meltome_ids

% add predicted optimal temperatures if available and set fN_TH at 100 °C
% to 1e-6 for these proteins
topt_prediction_file = config('topt_prediction_file');
topt_pred_table = readtable(topt_prediction_file);
for i = 1:numel(model.T_opt)
    if isnan(model.T_opt(i))
        if is_ec
            model.T_opt(i) = topt_pred_table.Topt(...
                ismember(topt_pred_table.ID, ec_model_prot_ids(i)));
        else
            model.T_opt(i) = topt_pred_table.Topt(...
                ismember(topt_pred_table.ID, model.geneUniprotID(i)));
        end
        model.T_H(i) = 100;
        model.meltingCurves.fN_TH(i) = 1e-6;
    end
end

clear topt_prediction_file
fprintf('Done\n')

%% Total protein content
fprintf('Fitting protein content...')

% prepare data to fit the model
proteinContentFile = config('proteinContentFile');
protContentTable = readtable(proteinContentFile,...
    'ReadVariableNames',false,'FileType','text');
% experimental temperatures
X = table2array(protContentTable(1,3:end));
% protein content measurements
Y = table2array(protContentTable(2:end,3:end));
% if dry weight to fresh weight conversion factor is given, convert to dry
% weight
DW2FW = config('DW2FW');
if exist('DW2FW','var') || ~isempty(DW2FW) || ~isnan(DW2FW)
    Y = Y/DW2FW;
end

% convert mg to g
Y = Y/1000; % [g/gFW]

% fit to gamma-pdf-like function
modelFunction = @(F,T) (1/F(1).^F(2)) .* T.^(F(2)-F(3)) .* exp(-T/F(1));
nlm = nlinfit(X,mean(Y,'omitnan'),modelFunction,[1 1 1]);
proteinModel = @(T)modelFunction(nlm,T);

model.proteinModel = proteinModel;

clear X Y modelFunction nlm proteinModel proteinContentFile
fprintf('Done\n')

%% Add gas exchange / photosynthesis parameters
fprintf('Adding gas exchange parameters...')
model.g_m = config('g_m');
model.g_s = config('g_s');
model.r_b = config('r_b');
model.C_a = config('C_a');
model.O_a = config('O_a');
model.K_c = config('K_c');
model.K_o = config('K_o');
model.k_c = config('k_c');
model.k_o = config('k_o');
model.V_c_max = config('V_c_max');
model.S_co = config('S_co');
model.E_a = config('E_a');
model.J_max = config('J_max');
model.alpha = config('alpha');
model.absorptance = config('absorptance');
model.P = config('P');
model.sqcf = config('sqcf');
model.C_ID = config('C_ID');
model.O_ID = config('O_ID');
model.RESP_ID = config('RESP_ID');
model.CO2_IMP_ID = config('CO2_IMP_ID');
model.ABS_ID = config('ABS_ID');

fprintf('Done\n');

%% Leaf mass per area
fprintf('Fitting model for leaf mass per area...')
lmaFile = config('lmaFile');
lmaTable = readtable(lmaFile, 'ReadVariableNames',false,'FileType','text');
X = celsius2kelvin(repmat(table2array(lmaTable(1,:)),1,size(lmaTable,1)-1));
Y = reshape(table2array(lmaTable(2:end,:))',...
    1,numel(table2array(lmaTable(2:end,:))));

f_sig = @(F,T) F(1) ./ (1+exp(F(2)-F(3)./T)) + F(4);
f_fit_sig = nlinfit(X,Y,f_sig,[max(Y) 1 1 min(Y)]);
model.lma_model = @(T) f_sig(f_fit_sig,T);

clear X Y lmaTable f_sig f_fit_sig lmaFile

fprintf('Done\n');

%% Fit MMRT parameters
fprintf('Fitting parameters of MMRT model using key temperatures...\n')
if is_ec
    [model.DH,model.DS,model.DCp] = fitMMRTGECKOModel(model);
else
    warning('MMRT fitting for standard GEMs is currently not implemented')
end
fprintf('Done\n')

TGEM = model;
end