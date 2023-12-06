function K_o_adj = adj_Ko(varargin)
%% K_o_adj = adj_Ko(varargin)
% Adjust the Michaelis-Menton constant of RuBisCO oxygenase activity to temperature.
% Function and default parameters were taken from
%       Walker et al. (2013), Plant, Cell & Environment
%       https://doi.org/10.1111/pce.12166
% Default values are taken from Table 1.
% 
% Input:
%   double T:           temperature(s) in Kelvin
%   double c:           scaling constant
%                       default: 14.72
%   double DHa:         activation energy [J mol^-1]
%                       default: 29100 J mol^-1
% Output:
%   double K_o_adj:     adjusted K_o value

T_DEFAULT = 298.15;     % [K]
R = 8.3145;             % [J mol^-1 K^-1]
C_DEFAULT = 14.72;      % unitless
DHA_DEFAULT = 29100;    % [J mol^-1]

p = inputParser;

addOptional(p,'T',T_DEFAULT)
addOptional(p,'c',C_DEFAULT)
addOptional(p,'DHa',DHA_DEFAULT)

p.parse(varargin{:})

T = p.Results.T;
c = p.Results.c;
DHa = p.Results.DHa;

% adjust K_o
K_o_adj = 10000*exp(c - DHa/R./T); % [ubar]

end