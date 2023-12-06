% Create an ecModel from the refined AraCore model using GECKO 2.0
clear;clc
% add the path of the RAVEN toolbox
addpath(genpath(fullfile('C://MATLAB/RAVEN/')))
cd C:\Users\wende\MobaXterm\home\GECKO
model = readCbModel('../ArabidopsisCoreModel/AraCore_v2_1/AraCore_v2_1.mat');

toolbox = 'COBRA';
modelName = 'ecModel of the ArabidopsisCoreModel v2.1';
modelVer = '1.0.0';

% re-stucture GPR rules
to_correct = {'PSII_h', 'PSI_h', 'cplx1_m', 'cplx2_m', 'cplx3_m', 'cplx4_m',...
    'cplx5_m', 'ATPase_h', 'RBC_h', 'RBO_h', 'AGPase_h', 'bAMY1_h', 'bAMY2_h',...
    'PyrK_h', 'iCitDHNAD_m', 'SCACoAL_m', 'ATPCitL_c', 'ACoAC_h', '3IPMDA1_h',...
    '3IPMDA2_h', 'ANTS_h', 'TrpS_h', 'ADPR_c', 'CDPR_c', 'GDPR_c', 'UDPR_c'};

for i = 1:numel(to_correct)
    rxn_idx = findRxnIDs(model, to_correct(i));
    model.rules{rxn_idx} = connectGPRByOR(model.rules{rxn_idx});
end

model.rules = strtrim(model.rules);

% Update biomass reaction to include pseudometabolites for metabolite
% classes of mian biomass components
mainComponents = struct;
mainComponents.protein = {'Ala[c]', 'Arg[c]', 'Lys[c]', 'His[c]', 'Ile[c]', 'Leu[c]',...
    'Phe[c]', 'Trp[c]', 'Tyr[c]', 'Val[c]', 'Met[c]', 'Gly[c]', 'Thr[c]',...
    'Cys[c]', 'Ser[c]', 'Asn[c]', 'Glu[c]', 'Asp[c]', 'Gln[c]', 'Pro[c]'};
mainComponents.DNA = {'dATP[c]', 'dTTP[c]', 'dGTP[c]', 'dCTP[c]'};
mainComponents.RNA = {'ATP[c]', 'UTP[c]', 'GTP[c]', 'CTP[c]'};
mainComponents.carbohydrate = {'starch2[h]', 'cellulose2[c]', 'Glc[c]',...
    'Suc[c]', 'Frc[c]', 'Tre[c]', 'Mas[c]','SCA[m]', 'Fum[m]', 'Mal[m]', 'SA[h]','Orn[h]'};
mainComponents.lipid = {'M-ACP[h]'};

fNames = fieldnames(mainComponents);
for i = 1:numel(fNames)
    
    % current main component name
    compName = fNames{i};
    
    % metabolites that belong to the current class
    currentIDs = mainComponents.(compName);
    currentIdx = findMetIDs(model,currentIDs);
    
    % respective biomass coefficients
    coeff = model.S(currentIdx,model.c==1);
    
    % add compartment to pseudometabolite
    compName = [compName '[c]'];
    
    % calculate sum formula of main component while respecting the
    % stoichiometric coefficients of the individual metabolites
    compFormula = '';
    for j=1:numel(currentIdx)
        tmpFormula = model.metFormulas(currentIdx(j));
        [ematrix, elements] = getElementalComposition(tmpFormula);
        ematrix = abs(coeff(j))*ematrix;
        compFormula = [compFormula, char(elementalMatrixToFormulae(ematrix, elements))];
    end
    % compute the sum formula for the component
    [ematrix, elements] = getElementalComposition(compFormula);
    compFormula = elementalMatrixToFormulae(ematrix, elements);
    
    % add pseudometabolite to model
    model = addMetabolite(model, compName,...
        'metName', strtok(compName,'['),...
        'metFormula', compFormula);
    
    clear ematrix elements tmpFormula compFormula
    
    % add pseudoreaction to model
    model = addReaction(model, strtok(compName,'['),...
        'reactionName', [strtok(compName,'['), ' pseudoreaction'],...
        'metaboliteList', [currentIDs {compName}],...
        'stoichCoeffList', [coeff;1],...
        'reversible', false, 'printLevel', 0);
    
    % remove metabolites of current class from biomass reaction
    model.S(currentIdx, model.c==1) = 0;
    
    % add current class/pseudometabolite to biomass reaction
    model.S(findMetIDs(model, compName),model.c==1) = -1;
end
clear coeff tmpFormula compName currentIDs currentIdx i mainComponents compNames ...
    ACP_IDX

% update some metabolite names for better matching of kcat values
old_names = {'Adenosine diphosphate','Adenosine triphosphate',...
    'Orthophosphate', '(R)-3-phosphoglycerate', 'Erythrose 4-phosphate',...
    'Xylulose 5-phosphate', 'Ribulose 5-phosphate', '2,3,4,5-Tetrahydrodipicolinate',...
    '2,3-Dihydroxy-3-methylpentanoate', '2-Deoxyadenosine 5-diphosphate',...
    '2-Deoxyadenosine 5-triphosphate', '2-Deoxycytidine 5-diphosphate',...
    '2-Deoxycytidine 5-triphosphate', '2-Deoxyguanosine 5-diphosphate',...
    '2-Deoxyguanosine 5-triphosphate', '2-Deoxythymidine 5-diphosphate',...
    '2-Deoxythymidine 5-monophosphate', '2-Deoxythymidine 5-triphosphate',...
    '2-Deoxyuridine 5-diphosphate', '2-Deoxyuridine 5-monophosphate',...
    '2-Deoxyuridine 5-triphosphate', '3-Deoxy-D-arabino-heptulosonate-7-phosphate',...
    '5,6,7,8-Tetrahydrofolate', '5-Amino-1-(5-Phospho-D-ribosyl)imidazole-4-carboxamide',...
    '5-Enolpyruvyl-shikimate-3-phosphate', '5-Phosphoribosyl 1-pyrophosphate',...
    '5-Phosphoribosyl-AMP', '5-Phosphoribosyl-AMP', '5-Phosphoribosyl-anthranilate',...
    'Acetyl-Coenzyme A', 'Acetylglutamate', 'Adenosine 5-phosphate',...
    'Adenosylhomocysteine', 'Adenosylmethionine', 'Aspartyl-4-phosphate',...
    'Citrulline', 'Cytidine diphosphate', 'Cytidine triphosphate',...
    'D-Fructose 2,6-bisphosphate', 'Diaminopimelate',...
    'Erythro-1-imidazol-4-glycerol 3-phosphate', 'Glutamate 5-semialdehyde',...
    'Guanosine 5-phosphate', 'Guanosine diphosphate', 'Guanosine triphosphate',...
    'Histidinal', 'Histidinol', 'Homoserine', 'Homoserine', 'Indole-glycerol 3-phosphate',...
    'Inosine 5''phosphate', 'Malonyl-Coenzyme A', 'Malonyl-acyl-carrier-protein',...
    'Orotidine 5-phosphate', 'Sucrose-6-phosphate', 'Uridine 5-phosphate',...
    'Uridine diphosphate', 'Uridine triphosphate', 'Xanthosine 5-phosphate'};


new_names = {'adp', 'atp', 'phosphate', '3-phosphoglycerate', 'D-Erythrose 4-phosphate',...
    'D-Xylulose 5-phosphate', 'D-Ribulose 5-phosphate', 'Tetrahydrodipicolinate',...
    '2R,3R-2,3-dihydroxy-3-methylpentanoate', 'dadp', 'datp', 'dcdp', 'dctp',...
    'dgdp', 'dgtp', 'dtdp', 'dtmp', 'dttp', 'dudp', 'dump', 'dutp',...
    '3-Deoxy-D-arabino-heptulosonate 7-phosphate', 'Tetrahydrofolate',...
    '5-formamido-1-(5-phosphoribosyl)imidazole-4-carboxamide',...
    '5-enolpyruvylshikimate 3-phosphate', '5-phosphoribosyl 1-diphosphate',...
    '1-(5-phosphoribosyl)-amp', '1-(5-phosphoribosyl)-amp',...
    'n-(5-phospho-beta-d-ribosyl)anthranilate', 'acetyl-coa', 'n-acetyl-l-glutamate',...
    'amp', 's-adenosyl-l-homocysteine', 's-adenosyl-l-methionine',...
    'aspartyl phosphate', 'L-Citrulline', 'cdp', 'ctp', 'fructose 2,6-bisphosphate',...
    'll-diaminopimelate', 'd-erythro-1-(imidazol-4-yl)glycerol 3-phosphate',...
    'l-glutamate 5-semialdehyde', 'gmp', 'gdp', 'gtp', 'l-histidinal', 'l-histidinol',...
    'l-homoserine', 'l-homoserine', '(1s,2r)-1-c-(indol-3-yl)glycerol 3-phosphate',...
    'imp', 'malonyl-coa', 'malonyl-acp', 'orotidine 5''-phosphate',...
    'd-sucrose 6-phosphate', 'ump', 'udp', 'utp', 'xmp'};

for i = 1:numel(old_names)
    model.metNames(ismember(model.metNames, old_names{i})) = new_names(i);
end

% update ProtDatabase.mat
updateDatabases

% A breakpoint was set within the GECKO function "getConstrainedModel.m"
% in line 55 (in the two previous lines, I inserted two lines to save the
% raw models).
% If the breakpoint is not set, GECKO will continue with the
% adjustment/correction of kcat values, which will not finish due to the
% low predicted growth rate.
model = rmfield(model,'grRules');
cd('geckomat')
[ecModel,ecModel_batch] = enhanceGEM(model,toolbox,modelName,modelVer);
cd('..')

% revert metabolite names
load('models\ecAraCoreModel\raw_batch_ecModel.mat')
load('models\ecAraCoreModel\raw_ecModel.mat')
for i = 1:numel(old_names)
    ecModel.metNames(ismember(ecModel.metNames, new_names{i})) = old_names(i);
    ecModel_batch.metNames(ismember(ecModel_batch.metNames, new_names{i})) = old_names(i);
end

inchikeys = repmat({''}, size(ecModel.mets));
inchikeys(1:numel(model.metisinchikeyID)) = model.metisinchikeyID;
ecModel.metisinchikeyID = inchikeys;

inchikeys = repmat({''}, size(ecModel_batch.mets));
inchikeys(1:numel(model.metisinchikeyID)) = model.metisinchikeyID;
ecModel_batch.metisinchikeyID = inchikeys;

save('models\ecAraCoreModel\raw_batch_ecModel.mat', 'ecModel_batch')
save('models\ecAraCoreModel\raw_ecModel.mat', 'ecModel')