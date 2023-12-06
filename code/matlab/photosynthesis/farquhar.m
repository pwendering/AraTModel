function [A,R_d,phi,gamma] = farquhar(varargin)
% Implementation of the C3 photosynthesis model originally published by
% Farquhar et al. (1980). If no input is given, runs with parameters from
% the original publication at 25 °C. Provide input using key-value pairs.
% Input (optional):
%   double T:           temperature [K]
%   double O:           O2 partial pressure [ubar]
%   double C:           CO2 partial pressure [ubar]
%   double I:           irradiance [uE m-2 s-1]
%   double V_c_max:     maximum carboxylation velocity at 25 °C [umol m-2 s-1]
%   double K_c:         Michaelis constant for CO2 at 25 °C [ubar]
%   double K_o:         Michaelis constant for O2 at 25 °C [ubar]
%   double kc:          turnover number of RuP2 carboxylase at 25 °C [s-1]
%   double ko:          turnover number of RuP2 oxygenase at 25 °C [s-1]
%   double J_max:       light saturated potential rate of electron 
%                       transport at 25 °C [uEq m-2 s-1]
%   double R_d:         respiration in the absence of light at 25 °C [umol m-2 s-1]
%   struct E_a:         activation energies for photosynthetic 
%                       parameters [J mol-1]
%                       possible fields: K_c, K_o, k_c, V_c_max, R_d, J_max
%   double TPU:         triose phosphate supply [umol m-2 s-1]
% 
% Output:
%   double A:           net CO2 assimilation rate [umol m-2 s-1]
%   double R_d:         respiration in the absence of light [umol m-2 s-1]
%   double phi:         oxygenation to carboxylation ratio
%   double gamma:       CO2 compensation point [ubar]
% 
% references:
%   Farquhar, G., von Caemmerer, S., Berry, J., A biochemical model of 
%       photosynthetic CO2 assimilation in leaves of C3 species
%       Planta 149(1), 78-90 (1980)
%   von Caemmerer S., Farquhar G., Berry J. Biochemical Model of C3 
%       Photosynthesis. In: Laisk A., Nedbal L., Govindjee (eds) 
%       Photosynthesis in silico. Advances in Photosynthesis and 
%       Respiration, vol 29. Springer, Dordrecht, 209-230 (2009)
    
% temperature [K]
defaultTemperature = 298.15;
% partial pressure of O2 [ubar]
defaultpO2 = 210000; 
% intercellular partial pressure of CO2 [ubar]
defaultpCO2 = 230;
% Irradiance [uE m-2 s-1]
defaultIrradiance = 1000;
% maximum carboxylation velocity [umol m-2 s-1]
defaultVcmax = 98;
% Michaelis constant for CO2 [ubar]
defaultKc = 460;
% Michaelis constant for O2 [ubar]
defaultKo = 330000;
% turnover number of RuP2 carboxylase [s-1]
defaultkc = 2.5;
% turnover number of RuP2 oxygenase [s-1]
defaultko = 0;
% light saturated potential rate of electron transport [uEq m-2 s-1]
defaultJmax = 210;
% "dark respiration" rate [umol m-2 s-1]
defaultRd = 1.1;
% triose phosphate supply []
defaultTPU = NaN;

% fraction of light not absorbed by chloroplasts
f = 0.15;
% empirical curvature factor for maximum electron transport rate
theta = 0.7;
% absorptance of leaves
absorptance = 0.85;

% activation energies for photosynthetic parameters [J mol-1]
defaultEa = struct;
defaultEa.K_c = 59356;
defaultEa.K_o = 35948;
defaultEa.k_c = 58520;
defaultEa.k_o = 58520;
defaultEa.V_c_max = 58520;
defaultEa.R_d = 66405;
defaultEa.J_max = 37000;

% free enthalpy [J mol-1]
H = 220000;
% entropy [J mol-1 K-1]
S = 710;
% universal gas constant [J mol-1 K-1]
R = 8.314;

% parse input
p = inputParser;

validTemp = @(x) isnumeric(x) && all(x >= 273.15) && all(x <= 333.15);

addParameter(p,'T',defaultTemperature,validTemp)
addParameter(p,'O',defaultpO2,@isnumeric)
addParameter(p,'C',defaultpCO2,@isnumeric)
addParameter(p,'I',defaultIrradiance,@isnumeric)
addParameter(p,'V_c_max',defaultVcmax,@isnumeric)
addParameter(p,'K_c',defaultKc,@isnumeric)
addParameter(p,'K_o',defaultKo,@isnumeric)
addParameter(p,'kc',defaultkc,@isnumeric)
addParameter(p,'ko',defaultko,@isnumeric)
addParameter(p,'J_max',defaultJmax,@isnumeric)
addParameter(p,'R_d',defaultRd,@isnumeric)
addParameter(p,'E_a',defaultEa,@isstruct)
addParameter(p,'TPU',defaultTPU,@isnumeric)

parse(p,varargin{:});

T = p.Results.T;
O = p.Results.O;
C = p.Results.C;
I = p.Results.I;
V_c_max = p.Results.V_c_max;
K_c = p.Results.K_c;
K_o = p.Results.K_o;
k_c = p.Results.kc;
k_o = p.Results.ko;
J_max = p.Results.J_max;
R_d = p.Results.R_d;
E_a = p.Results.E_a;
TPU = p.Results.TPU;

% fill missing activation energies with default values
fields = fieldnames(defaultEa);
for i=1:numel(fields)
    % check if field is present and if a value different from zero is assigned
    if ~isfield(E_a,fields{i}) || E_a.(fields{i})==0 || isempty(E_a.(fields{i}))
        E_a.(fields{i}) = defaultEa.(fields{i});
    end
end
        
% adjust Michaelis-Menten constants to given temperature
% K_c = adjustParameterToTemperature('p_ref',K_c,'Ea',E_a.K_c,'T_dest',T);
K_c = adj_Kc('T',T);
% K_o = adjustParameterToTemperature('p_ref',K_o,'Ea',E_a.K_o,'T_dest',T);
K_o = adj_Ko('T',T);
% adjust turnover number of RuP2 carboxylase to given temperature
k_c = adj_k_c('k25',k_c,'T',T);
% k_c = adjustParameterToTemperature('p_ref',k_c,'Ea',E_a.k_c,'T_dest',T);
% turnover number of RuP2 oxygenase [s-1]
if k_o==0
    k_o = 0.21*k_c;
else
%     k_o = adjustParameterToTemperature('p_ref',k_o,'Ea',E_a.k_o,'T_dest',T);
    k_o = adj_k_o('k25',k_o,'T',T);
end

% adjust maximum carboxylation velocity to given temperature
% V_c_max = adjustParameterToTemperature('p_ref',V_c_max,'Ea',E_a.V_c_max,'T_dest',T);
V_c_max = adj_v_c_max('k25', V_c_max, 'T', T);
% adjust "dark respiration" rate to given temperature
R_d = adjustParameterToTemperature('p_ref',R_d,'Ea',E_a.R_d,'T_dest',T);
% adjust light saturated potential rate of electron transport [uEq m-2 s-1]
% J_max = J_max*exp(((T-298.15)./298.15)*(E_a.J_max)./(R*T)) ...
%     .* ((1+exp((298.15*S-H)/298.15*R)) ...
%     ./ (1+exp((S.*T-H)./(R.*T))));
J_max = adj_j_max('J_max_ref', J_max, 'T', T);
% ratio of oxygenation to carboxylation
phi = (k_o./k_c).*(O./K_o)./(C./K_c);

% calculate CO2 compensation point, without dark respiration [ubar]
gamma_star = (K_c*O.*k_o) ./ (2*K_o.*k_c);
gamma = (gamma_star + K_c.*(1+O./K_o).*(R_d./V_c_max)) ./ (1 - R_d./V_c_max);

% determine potential rate of electron transport
I_2 = I*absorptance*(1-f)/2;
J = (I_2 + J_max - sqrt((I_2+J_max).^2 - 4*theta*I_2.*J_max)) / (2*theta);

% RuP2 saturated rate (RubBisCO-limited)
W_c = (V_c_max * C) ./ (C + K_c.*(1 + O./K_o));

% electron transport limited rate
W_j = J ./ (4 + 8*gamma_star/C);

% triose phosphate utilization limited rate
W_t = 3*TPU ./ (1 - gamma_star/C);

A = (1 - gamma_star/C).*min([W_c;W_j;W_t],[],1) - R_d;
end