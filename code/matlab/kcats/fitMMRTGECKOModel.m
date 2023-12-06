function [DH,DS,DCp] = fitMMRTGECKOModel(model,T0)
%% [DH,DS,DCp] = fitMMRTGECKOModel(model,T0)
% Fit parameters of the Macromolecular rate theory (MMRT) model for temperature
% dependency of reaction rates based on key temperatures of protein
% stability, i.e. optimal temperature, melting temperature, and heat
% denaturation temperature.
% By default, a reference temperature of 20 °C is used.
% Reference: Hobbs et al. (2013), https://doi.org/10.1021/cb4005029
%
% Input:    struct model:       GECKO-formatted metabolic model with
%                               additional fields: T_opt,T_m, and T_H
%           scalar T0:          (opt) reference temperature [K]
% Output:   double DH:          difference in enthalpy between the ground and
%                               transition state
%           double DS:          difference in entropy between the ground and
%                               transition state
%           double DCp:         difference in heat capacity between the ground and
%                               transition state
%
% Philipp Wendering, University of Potsdam (philipp.wendering@gmail.com) 23/03/2023

%% Constants and options
NPROT = numel(model.enzymes);
NRXNS = numel(model.rxns);
ENZ_PFX = 'prot_';

% Boltzmann constant to [J K^-1]
kB = config('kB');
% Planck's constant [J s]
h = config('h');
% universal gas constant [J K^-1 mol^-1]
R = config('R');

% set reference temperature
if nargin < 2; T0 = celsius2kelvin(20); end

% set default for kcat multiplication factor at heat denaturation
% temperature, T_H
% set default for kcat multiplication factor at heat denaturation
% temperature, T_H
if isfield(model,'meltingCurves') && isfield(model.meltingCurves,'fN_TH')
    pNativeTH = model.meltingCurves.fN_TH;
    pNativeTH(pNativeTH<=0) = 1e-30;
else
    pNativeTH = repelem(1e-30,1,NPROT);
end

%% Initialize output arrays
% activation energy contributions
DH = sparse(NPROT,NRXNS);
DS = sparse(NPROT,NRXNS);
DCp = sparse(NPROT,NRXNS);

% use key temperatures that were estimated from melting curves
% (if not possible to determine, use median over all proteins)
T_opt = model.T_opt;
T_opt = celsius2kelvin(T_opt);

T_H = model.T_H;
T_H(~isreal(T_H)) = NaN;
T_H = celsius2kelvin(T_H);

% find row indices of enzymes
enzRowIdx = find(startsWith(model.mets,ENZ_PFX));
enzRowIdx = setdiff(enzRowIdx,find(ismember(model.mets,'prot_pool')));

%% loop over all proteins and fit EAAR parameters
for i=1:NPROT
    
    % check correct order of key temperatures
    if ~any(isnan([T_opt(i) T_H(i)]))
        
        % set up system of linear equations (LHS)
        A = [-1   T_opt(i)  T0-T_opt(i)+T_opt(i)*log(T_opt(i)/T0);
             -1   T_H(i)    T0-T_H(i)+T_H(i)*log(T_H(i)/T0);
             -1    0        T0-T_opt(i)
            ];
        
        % find reactions associated with the current protein
        protRxnsIdx = find(model.S(enzRowIdx(i),:)<0);
        enzKcatValues = -1./model.S(enzRowIdx(i),protRxnsIdx)/3600; %[s^-1]
        
        for j=1:numel(protRxnsIdx)
            
            % set RHS of equations using the appropriate kcat value
            B = [R*T_opt(i)*(log(enzKcatValues(j))-log(kB*T_opt(i)/h));
                 R*T_H(i)*(log(enzKcatValues(j)*pNativeTH(i))-log(kB*T_H(i)/h));
                 T_opt(i)*R];
             
            % solve system
            solution = rref([A,B]);
            DH(i,protRxnsIdx(j)) = solution(1,end);
            DS(i,protRxnsIdx(j)) = solution(2,end);
            DCp(i,protRxnsIdx(j)) = solution(3,end);
            
        end
    end
end

end
