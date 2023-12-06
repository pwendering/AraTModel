function K_c_adj = adj_Kc(varargin)
%% K_c_adj = adj_Kc(varargin)
% Adjust the Michaelis-Menton constant of RuBisCO carboxylase activity to temperature.
% Function and default parameters were taken from
%       Walker et al. (2013), Plant, Cell & Environment
%       https://doi.org/10.1111/pce.12166
% Default values are taken from Table 1.
% 
% Input:
%   double T:           temperature(s) in Kelvin
%   double c:           scaling constant
%                       default: 23.32
%   double DHa:         activation energy [J mol^-1]
%                       default: 49700 J mol^-1
% Output:
%   double K_c_adj:     adjusted K_c value

T_DEFAULT = 298.15;     % [K]
R = 8.3145;             % [J mol^-1 K^-1]
C_DEFAULT = 23.32;      % unitless
DHA_DEFAULT = 49700;    % [J mol^-1]

p = inputParser;

addOptional(p,'T',T_DEFAULT)
addOptional(p,'c',C_DEFAULT)
addOptional(p,'DHa',DHA_DEFAULT)

p.parse(varargin{:})

T = p.Results.T;
c = p.Results.c;
DHa = p.Results.DHa;

% adjust K_c
K_c_adj = 10 * exp(c - DHa/R./T); % [ubar]

end