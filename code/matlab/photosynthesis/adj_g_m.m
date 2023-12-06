function g_m = adj_g_m(varargin)
%% g_m = adj_g_m(varargin)
% Estimate the mesophyll conductance using either an empirical model by
% fitting experimental data to a gaussian function:
%               g_m(T) = a * exp(-((T-b)/c)^2),
% a mechanistical model by von Caemmerer and Evans 2015 (DOI:10.1111/pce.12449):
%               g_m(T) = (1/g_liq + 1/g_mem)^-1,
% a Q10-based model by Niinemets et al. 2009 (DOI:10.1093/jxb/erp063):
%               g_m(T) = g_m_ref * Q10^(T/10 - T_ref/10),
% or an empirical model described in Walker et al. 2013 (DOI:10.1111/pce.12166)
% (There is a mistake in the formula, but the correct version can be found
% in Bernacchi et al. 2002, Plant Physiology, DOI:10.1104/pp.008250
%               g_m(T) = exp(c - DHa/RT) / (1 + exp( (DS*T - DHd)/RT ))
%
% The default values for E, S_c, l, g_m_25 (ref) are specific for Arabidopsis thaliana.
%
% Input:
%       double T:       (optional) temperature [K]
%       double T_ref:   (optional) reference temperature for given g_m [K];
%                           default 293.15 K
%       double S_c:     (optional) chloroplast surface area appressing
%                           intercellular airspace per unit leaf area
%       double l:       (opional) effective path length of CO2 diffusion [m]
%       double E:       (optional) activation energy [J mol^-1]
%       double method:  (optional) either 'walker' (default),
%                           'caemmerer', 'niinemets', and 'empirical'
%       double Q10:     (optional) Q10 value for niinemets method
%       double g_m_ref: (optional) mesophyll conductance at T_ref
%       double c:       (optional) scaling constant (unitless)
%       double DHa:     (optional) activation energy [J mol^-1]
%       double DHd:     (optional) deactivation energy [J mol^-1]
%       double DS:      (optional) entropy differnce [J mol^-1 K^-1]
%
% Output:
%       double g_m:     mesophyll conductance [mol m^-2 s^-1 bar^-1]

% temperature [K]
default_T = celsius2kelvin(25);

% universal gas constant [J K^-1 mol^-1]
R = 8.31446261815324;

% chloroplast surface area appressing intercellular airspace per unit leaf
% area
default_S_c = 7;

% effective path length [m]
default_l = 0.62 * 1E-6;

% activation energy [J mol^-1]
default_E = 68.4 * 1000; % fitted over range 15–25 °C

% default method
default_method = 'walker';

% default Q10 value for niinemets method
default_q10 = 2;

% mesophyll conductance at T_ref for Niinemets method
default_g_m_ref = 0.2;

% scaling constant (unitless)
default_c = 3;

% activation energy [J mol^-1]
default_DHa = 7400;

% deactivation energy [J mol^-1]
default_DHd = 434000;

% entropy differnce [J mol^-1 K^-1]
default_DS = 1400;

% parse input
validTemp = @(x) isnumeric(x) && all(x >= 273.15) && all(x <= 353.15);
validMethod = @(x) ismember(x,...
    {'empirical', 'caemmerer', 'niinemets', 'walker'});

p = inputParser;

addParameter(p,'T', default_T, validTemp)
addParameter(p,'T_ref',default_T, validTemp)
addParameter(p,'S_c', default_S_c, @isnumeric)
addParameter(p,'l', default_l, @isnumeric)
addParameter(p,'E', default_E, @isnumeric)
addParameter(p,'method', default_method, validMethod)
addParameter(p,'Q10', default_q10, @isnumeric)
addParameter(p,'g_m_ref', default_g_m_ref, @isnumeric)
addParameter(p,'c', default_c, @isnumeric)
addParameter(p,'DHa', default_DHa, @isnumeric)
addParameter(p,'DHd', default_DHd, @isnumeric)
addParameter(p,'DS', default_DS, @isnumeric)

parse(p,varargin{:})

T = p.Results.T;
T_ref = p.Results.T_ref;
S_c = p.Results.S_c;
l = p.Results.l;
E = p.Results.E;
method = p.Results.method;
Q10 = p.Results.Q10;
g_m_ref = p.Results.g_m_ref;
c = p.Results.c;
DHa = p.Results.DHa;
DHd = p.Results.DHd;
DS = p.Results.DS;

switch method
    case 'caemmerer'
        % combined membrane permeability to CO2 at 25 °C [m s^-1]
        P_mem_25 = 1.4 * 1E-3;
        
        % diffusivity of CO2 in water [m^2 s^-1]
        D = 1.81 * 1E-6 * exp(-16900 ./ (R * T));
        
        % molar density of water [mol m^-3] times Henry coefficient for CO2 [bar^-1]
        rho_x_H = 33.06 * exp(2400*(1./T - 1/298.15));
        
        % conductance through liquid phase [mol m^-2 s^-1 bar^-1]
        g_liq = rho_x_H .* D / l;
        
        % conductance through membrane [mol m^-2 s^-1 bar^-1]
        g_mem = rho_x_H * P_mem_25 .* exp((T-298.15)*E ./ (R * 298.15 * T));
        
        % mesophyll conductance per chloroplast surface area [mol m^-2 s^-1 bar^-1]
        g_m = 1./(1./g_liq + 1./g_mem);
        
        % mesophyll conductance per leaf area [mol m^-2 s^-1 bar^-1]
        g_m = g_m * S_c;
        
    case 'niinemets'
        g_m = g_m_ref * Q10.^(T./10 - T_ref./10);
    case 'empirical'
        % experimental data from Caemmerer and Evans (2015), DOI:10.1111/pce.12449
        T_data = celsius2kelvin([15;20;25;30;35;40]);
        g_m_data = [0.133;0.175;0.221;0.220;0.194;0.124];     
        
        % fit to gaussian curve
        gauss_fit = fit(T_data,g_m_data, 'gauss1');
        % function that corresponds to 'gauss1'
        f = @(b,T) b(1)*exp(-((T-b(2))/b(3)).^2);
        beta = [gauss_fit.a1 gauss_fit.b1, gauss_fit.c1];
        % estimate conductance at given temperature
        g_m = f(beta,T);
        
    case 'walker'
        g_m = g_m_ref*exp(c-DHa/R./T)./(1+exp((DS.*T-DHd)/R./T));
end

end