function [VPD, SVP] = vpd(varargin)
%% [VPD, SVP] = vpd(rh,T)
% Calculate the Vapour Pressure Deficit (VPD) at a given temperature and
% relative humidity.
% Reference for the calculation of the saturation vapour pressure (SVP):
%   Murray (1967), DOI: 10.1175/1520-0450(1967)006<0203:OTCOSV>2.0.CO;2
% Input
%   double rh:        relative humidity in %
%   double T:         temperature in K
% Output:
%   double VPD:       VPD in Pa
%   doubel SVP:       SVP in Pa

DEFAULT_T = 293.15;
DEFAULT_RH = 70;

validTemp = @(T) all(T>273.15 & T<373.15);
validRH = @(RH) RH>=0 & RH <= 100;
p = inputParser;

addParameter(p,'T',DEFAULT_T,validTemp)
addParameter(p,'rh',DEFAULT_RH,validRH)

parse(p,varargin{:});

T = p.Results.T;
rh = p.Results.rh;

if rh < 1
    warning('Relative humidity was given below 1, double-check if given as percentage')
end

% calculate saturation vapour pressure
T = kelvin2celsius(T);
SVP = 610*10.^(7.5*T./(237.3+T)); % Pa

% calculate VPD
VPD = (1-rh/100)*SVP;

end