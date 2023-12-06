function isValid = validateTGEM(model)
% Check if a model structure has all required fields for addition of
% temperature dependent constraints based on the FvCB model, total protein 
% content and kcat values.
% For more information see config.m and createTGEM.m .

% check if model has all required fields
requiredFields = {...
    ...         --- General model fields ---
    'S',...             % stoichiometric matrix
    'b',...             % right-hand side of steady state constraints
    'c',...             % objective vector
    'osenseStr',...     % 'min' for minimization, 'max' for maximization
    'lb',...            % upper bounds for reaction fluxes
    'ub',...            % lower bounds for reaction fluxes
    'rxns',...          % reaction identifiers
    'mets',...          % metabolite identifiers
    ...         --- FvCB model ---
    'C_a',...           % ambient CO2 partial pressure
    'O_a',...           % ambient O2 partial pressure
    'K_c',...           % Michaelis-Menten constant for CO2 (RuBisCO carboxylation)
    'K_o',...           % Michaelis-Menten constant for O2 (RuBisCO oxygenation)
    'k_c',...           % turnover number for CO2 (RuBisCO carboxylation)
    'k_o',...           % turnover number for O2 (RuBisCO oxygenation)
    'V_c_max',...       % maximum carboxylation velocity
    'J_max',...         % light saturated potential rate of electron transport
    'E_a',...           % activation energies for FvCB model parameters
    'lma_model',...     % temperature model for leaf mass per area
    'P',...             % atmospheric pressure
    'g_m',...           % mesophyll conductance
    'r_b',...           % boundary layer resistance
    'sqcf',...          % correction factor for the spectral quality of light
    'alpha',...         % number of CO2 molecules released per oxygenation
    'absorptance',...   % absorptance / quantum yield of red light
    'C_ID',...          % model ID for the RuBisCO carboxylation reaction
    'O_ID',...          % model ID for the RuBisCO oxygenation reaction
    'RESP_ID',...       % model ID for CO2 export out of mitochondria
    'ABS_ID',...        % model ID for photon import / exchange reaction
    ...         --- Total protein content ---
    'proteinModel',...  % temperature model for the total protein content
    ...         --- kcat adjustment ---
    'enzymes',...       % enzyme IDs (added by GECKO toolbox)
    'DH',...            % enthalpy term for MMRT-based kcat adjustment
    'DS',...            % entropy term for MMRT-based kcat adjustment
    'DCp'...            % heat capacity term for MMRT-based kcat adjustment
    };
    
fieldFound = cellfun(@(x)isfield(model,x),requiredFields);

if ~all(fieldFound)
    error('Model does not have all required fields\nMissing fields are: %s',...
        strjoin(requiredFields(~fieldFound), ', '))
end

% check model fields
[rowDim, colDim] = size(model.S);
if size(model.b,1) ~= rowDim
    error('RHS vector (model.b) has incorrect size.')
elseif numel(model.c) ~= colDim
    error('Objective vector (model.c) has incorrect size.')
elseif ~any(ismember({'min', 'max'}, model.osenseStr))
    error('model.osenseStr must be ''min'' or ''max''.')
elseif numel(model.lb) ~= colDim
    error('Lower bound vector (model.lb) has incorrect size.')
elseif numel(model.ub) ~= colDim
    error('Upper bound vector (model.ub) has incorrect size.')
elseif numel(model.rxns) ~= colDim
    error('Reaction ID vector (model.rxns) has incorrect size.')
elseif numel(model.mets) ~= rowDim
    error('Metabolite ID vector (model.mets) has incorrect size.')
end

% check FvCB parameter fields
tempRange = celsius2kelvin([10 25 40]);
lmaTest = model.lma_model(tempRange);
if any(lmaTest<0 | isnan(lmaTest))
    error('LMA model (model.lma_model) produces negative or NaN values between 10 °C and 40°C')
end

if ~findRxnIDs(model, model.C_ID)
    error('RuBisCO carboxylation reaction ID (model.C_ID) was not found.')
elseif ~findRxnIDs(model, model.O_ID)
    error('RuBisCO oxygenation reaction ID (model.O_ID) was not found.')
elseif ~findRxnIDs(model, model.RESP_ID)
    error('Mitochondrial CO2 export reaction ID (model.RESP_ID) was not found.')
elseif ~findRxnIDs(model, model.ABS_ID)
    error('Photon uptake reaction ID (model.ABS_ID) was not found.')
end
    
% check total protein content temperature model
tempRange = [10 25 40];
ptotTest = model.proteinModel(tempRange);
if any(ptotTest<0 | isnan(ptotTest))
    error(['Protein content model (model.proteinModel) produces negative '...
        'or NaN values between 10 °C and 40°C'])
end

% check kcat adjustment field (MMRT parameters)
nEnzymes = numel(model.enzymes);

% if the model contains TFA constraints added by matTFA toolbox
if any(startsWith(model.rxns, 'DGo_')) && any(startsWith(model.rxns, 'NF_'))
    % use number of "net flux variables" to check the correct column
    % dimension
    colDim = sum(startsWith(model.rxns, 'NF_'));
end

if size(model.DH,1) ~= nEnzymes || size(model.DH,2) ~= colDim
    error('Size of matrix containing MMRT enthalpy terms (model.DH) is incorrect')
elseif size(model.DS,1) ~= nEnzymes || size(model.DS,2) ~= colDim
    error('Size of matrix containing MMRT entropy terms (model.DS) is incorrect')
elseif size(model.DCp,1) ~= nEnzymes || size(model.DCp,2) ~= colDim
    error('Size of matrix containing MMRT heat capacity terms (model.DCp) is incorrect')
end

isValid = true;

end