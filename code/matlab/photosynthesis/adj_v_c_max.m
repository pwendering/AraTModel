function v_c_max_adj = adj_v_c_max(varargin)
%% v_c_max_adj = adj_v_c_max(varargin)
% Adjust the maximum carboxylation velocity of RuBisCO to a given temperature.
% Function and was taken from
%       Walker et al. (2013), Plant, Cell and Environment,
%       https://doi.org/10.1111/pce.12166
%       (Table 1)
% Default Parameter was taken from Heckwolf et al. (2011), Plant Journal,
%       https://doi.org/10.1111/j.1365-313x.2011.04634.x
%       (Table 2)
% 
% Input:
%   double T:               temperature(s) in Kelvin
%   double k25:             V_cmax value at 25 °C
%                           default: 96.30 umol m^-2 s^-1
%   double c:               scaling constant
%                           default: 16.66
%   double DHa:             activation energy [J mol^-1]
%                           default: 41400 J mol^-1
% Output:
%   double v_c_max_adj:     adjusted V_cmax value

T_DEFAULT = 298.15;     % [K]
K25_DEFAULT = 96.3; 	% [umol m^-2 s^-1]
C_DEFAULT = 16.66;      % unitless
DHA_DEFAULT = 41400;    % [J mol^-1]
R = 8.3145;             % [J mol^-1 K^-1]
p = inputParser;

addOptional(p,'T',T_DEFAULT)
addOptional(p,'k25',K25_DEFAULT)
addOptional(p,'c',C_DEFAULT)
addOptional(p,'DHa',DHA_DEFAULT)

p.parse(varargin{:})

T = p.Results.T;
k25 = p.Results.k25;
c = p.Results.c;
DHa = p.Results.DHa;

% adjust V_cmax
v_c_max_adj = k25 * exp(c - DHa/R./T);

end