function j_max_adj = adj_j_max(varargin)
%% j_max_adj = adj_j_max(varargin)
% Adjust the light saturated potential rate of electron transport to a given temperature.
% Function and default parameters were taken from the average J_max(T)
% model over multiple plant species described in
%       Leuning et al. (2002), DOI: 10.1046/j.1365-3040.2002.00898.x
% 
% Input:
%   double T:                   temperature(s) in Kelvin
%                               default: 298.15 K (25 °C)
%   scalar double J_max_ref:    J_max at the reference temperature [umol m^-2 s^-1]
%                               default: 138 umol m^-2 s^-1 (measured for
%                               A. thaliana by Gandin et al. (2012),
%                               DOI: 10.1093/pcp/pcs107)
%   scalar double H_a:          activation energy [J mol^-1]
%                               default: 50300 J mol^-1
%   scalar double H_d:          deactivation energy [J mol^-1]
%                               default: 152044 J mol^-1
%   scalar double S_v:          entropy term [J mol^-1]
%                               default: 495 J mol^-1 K^-1
%   scalar double T_ref:        reference temperature (298.15 K)
%                               default: 298.15 K (25 °C)
% Output:
%   double j_max_adj:           adjusted J_max value

T_DEFAULT = 298.15;     % [K]
JMAX_DEFAULT = 138.5;     % [umol m^-2 s^-1]
HA_DEFAULT = 50300;     % [J mol^-1]
HD_DEFAULT = 152044;    % [J mol^-1]
SV_DEFAULT = 495;       % [J mol^-1 K^-1]
Tref_DEFAULT = 298.15;  % [J mol^-1]
R = 8.3145;             % [J mol^-1 K^-1]
p = inputParser;

addOptional(p, 'T', T_DEFAULT)
addOptional(p, 'J_max_ref', JMAX_DEFAULT)
addOptional(p, 'H_a', HA_DEFAULT)
addOptional(p, 'H_d', HD_DEFAULT)
addOptional(p, 'S_v', SV_DEFAULT)
addOptional(p, 'T_ref', Tref_DEFAULT)

p.parse(varargin{:})

T = p.Results.T;
J_max_ref = p.Results.J_max_ref;
H_a = p.Results.H_a;
H_d = p.Results.H_d;
S_v = p.Results.S_v;
T_ref = p.Results.T_ref;

j_max_adj = J_max_ref .* ...
    (1+exp((S_v*T_ref-H_d)/R./T_ref)) .* ...
    exp((H_a/R/T_ref)*(1-T_ref./T)) ./ ...
    (1+exp((S_v*T-H_d)/R./T));

end