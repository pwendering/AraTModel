% Create of an enzyme- and temperature-constrained metabolic model for Arabidopsis
% thaliana

pathToTGEM = config('pathToTGEM');
tmp = load(fullfile(pathToTGEM, 'metabolic-models', 'ecAraCore_batch_raw.mat'));
model = tmp.ecModel_batch;

% add/update relevant fields for analysis with COBRA toolbox
model.rxnECNumbers = model.eccodes;
model.osenseStr = 'max';
model.csense = repelem('E', numel(model.mets), 1);

TGEM = createTGEM(model);
save(config('tgemFile'), 'TGEM');
