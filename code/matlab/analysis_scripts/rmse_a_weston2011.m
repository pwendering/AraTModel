% Calculate the RMSE as measure of (dis)agreement between the predicted net
% CO2 assimilation rates and the data set from Weston et al. 2011

% data from Weston et al. 2011
T = [24 26 28 30 32 35 36 38 40 42];
A = [10.33 13.29 15.64 17.24 17.67 17.08 15.61 13.51 11.29 8.55];
I = 800;

% load the temperature-constrained model
model_file = config('tgemFile');
tmp = load(model_file);
model = tmp.TGEM;
clear tmp

% add RubisCO activase module
model = addRCA(model);

% block uptake of H2S and chloroplastic alanine transaminase
model.ub(findRxnIDs(model, {'Im_H2S', 'Im_H2S_REV'})) = 0;
model.ub(findRxnIDs(model, {'AlaTA_h', 'AlaTA_h_REV'})) = 0;
model.ub(findRxnIDs(model, {'Bio_AA', 'Bio_CLim', 'Bio_NLim'})) = 0;

% correct THF transport reaction
model.S(findMetIDs(model, 'H_h'), findRxnIDs(model, {'Tr_THF_h', 'Tr_THF_h_REV'})) = [-1 1];

sol = simulateTempEffects(model, I, 'tempRange', celsius2kelvin(T));

A_pred_aracore = sol.A_max;

rmse_aracore = sqrt(mean((A-A_pred_aracore).^2));
[r_aracore, p_aracore] = corr(A', A_pred_aracore');

fprintf('RMSE ecAraCore to Weston et al. 2011: %.4g\n', rmse_aracore)
fprintf('Pearson correlation ecAraCore to Weston et al. 2011: %.4g (P=%.4g)\n',...
    r_aracore, p_aracore)

% FvCB model
fvcbParams = getFvCBParams();
% initialize assimulation rate matrix (T x p(CO2))
A_far = zeros(numel(T), 1);

for i=1:numel(T)
    A_far(i) = farquhar(...
        'T', celsius2kelvin(T(i)),...
        'C', 0.75*model.C_a,...
        'O', 210000,...
        'I', I,...
        'kc', fvcbParams.k_c,...
        'ko', fvcbParams.k_o,...
        'K_c', fvcbParams.K_c,...
        'K_o', fvcbParams.K_o,...
        'E_a', fvcbParams.E_a,...
        'V_c_max', fvcbParams.V_c_max,...
        'J_max', fvcbParams.J_max);
end

rmse_fvcb = sqrt(mean((A-A_far').^2));
[r_fvcb, p_fvcb] = corr(A', A_far);

fprintf('RMSE FvCB model to Weston et al. 2011: %.4g\n', rmse_fvcb)
fprintf('Pearson correlation FvCB model to Weston et al. 2011: %.4g (P=%.4g)\n', r_fvcb, p_fvcb)

