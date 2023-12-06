function s_co_adj = adj_s_co(varargin)
%% s_co_adj = adj_s_co(varargin)
% Adjust the RuBisCO specificity for CO2 over O2 (S_c/o) to a given temperature.
% Function and default parameters were taken from
%       Boyd et al. (2019), Journal of Experimental Botany,
%       https://doi.org/10.1093/jxb/ery355
% Default values are taken from Table 3 (radiolabel).
% 
% Input:
%   double T:           temperature(s) in Kelvin
%   scalar k25:         S_c/o value at 25 °C [bar bar^-1]
%                       default: 2003 bar bar^-1
%   scalar Ea:          activation energy [J mol^-1]
%                       default: -28660 J mol^-1
% Output:
%   double s_co_adj:    adjusted S_c/o value

T_DEFAULT = 298.15;     % [K]
K25_DEFAULT = 2003; 	% [bar bar^-1]
EA_DEFAULT = -28660;    % [bar bar^-1]
R = 8.3145;             % [J mol^-1 K^-1]
p = inputParser;

addOptional(p,'T',T_DEFAULT)
addOptional(p,'k25',K25_DEFAULT)
addOptional(p,'Ea',EA_DEFAULT)

p.parse(varargin{:})

T = p.Results.T;
k25 = p.Results.k25;
Ea = p.Results.Ea;

% check if all temperatures are between 5 and 40 °C
if any(T<278.15) || any(T>313.15)
    fprintf(['%s: The given temperature is outside the experimental ',...
        'range of 5-40 °C. Results may be inaccurate.\n'], mfilename)
end

% adjust S_c/o
s_co_adj = k25*exp(-Ea*(298.15-T)./(R*T*298.15));

end