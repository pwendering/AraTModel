function lnK = mmrtFun(T,DH,DS,DCp,T0)
%% nK = mmrtFun(T,DH,DS,DCp)
% Calculate the logarithm of the reaction rate at a given temperature using
% Macromolecular rate theory (MMRT, Hobbs et al. (2013), DOI: 10.1021/cb4005029
% The reference temperature is 20 °C by default but can be changed by
% setting T0.
% Input:
%       double T:       scalar or array of temperatures [K]
%       double DH:      difference in enthalpy between the ground and
%                       transition state
%       double DS:     	difference in entropy between the ground and
%                       transition state
%       double E_Cp:    difference in heat capacity between the ground and
%                       transition state
%       double T0:      (opt) reference temperature [K]
% Output:
%       double lnK:     scalar or array of logarithmized rates at T
if nargin < 5; T0=celsius2kelvin(20); end
% Boltzmann constant to [J K^-1]
kB = config('kB');
% Planck's constant [J s]
h = config('h');
% universal gas constant [J K^-1 mol^-1]
R = config('R');

lnK=log(kB*T/h)+(-DH+DCp*(T0-T+T.*log(T/T0))+DS*T)./(R*T);

end