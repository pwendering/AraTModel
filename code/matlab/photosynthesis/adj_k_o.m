function k_o_adj = adj_k_o(varargin)
%% k_o_adj = adj_k_o(varargin)
% Adjust the turnover number of RuBisCO oxygenase activity to temperature.
% Function and default parameters were taken from
%       Boyd et al. (2019), Journal of Experimental Botany,
%       https://doi.org/10.1093/jxb/ery355
% Default values are taken from Table 2 (MIMS).
% The function consideres the thermal breakpoint at 25 °C, by using the two
% different activation energies given in Boyd et al. (2019).
% 
% Input:
%   double T:           temperature(s) in Kelvin
%   scalar k25:         k_o value at 25 °C [s^-1]
%                       default: 1.38 s^-1
% Output:
%   double k_o_adj:     adjusted k_o value

T_DEFAULT = 298.15;     % [K]
K25_DEFAULT = 1.38; 	% [s^-1]
R = 8.3145;             % [J mol^-1 K^-1]
Ea_10_25 = 92950;       % [J mol^-1]
Ea_25_40 = 47110;       % [J mol^-1]

p = inputParser;

addOptional(p,'T',T_DEFAULT)
addOptional(p,'k25',K25_DEFAULT)

p.parse(varargin{:})

T = p.Results.T;
k25 = p.Results.k25;

% check if all temperatures are between 10 and 40 °C
if any(T<283.15) || any(T>313.15)
    fprintf(['%s: The given temperature is outside the experimental ',...
        'range of 10-40 °C. Results may be inaccurate.\n'], mfilename)
end

% adjust k_o
if T < celsius2kelvin(25)
    Ea = Ea_10_25;
else
    Ea = Ea_25_40;
end

k_o_adj = k25*exp(-Ea*(298.15-T)./(R*T*298.15));

end