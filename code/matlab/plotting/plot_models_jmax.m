% compare temperature models for Jmax
clear;clc;close all
%% Define simulation parameters and calculation of A
T_range = linspace(10,40,100);

C = 300;

I = 300;

f = 0.15;

theta = 0.7;

Rd = 1.1;

absorptance = 0.85;

I2 = I*absorptance*(1-f)/2;

gamma_star = 0.5*210000./adj_s_co(celsius2kelvin(linspace(10,40,100)));
% electron-transport-limited
A_fun = @(J_max)((I2 + J_max - sqrt((I2+J_max).^2 - 4*theta*I2*J_max)) / ...
    (2*theta).*(C-gamma_star))./(4*C+8*gamma_star)-Rd;

%% Bunce 2008
% Gaussian function
% a     value of Jmax at the optimum temperature
% b     steepness of the curve
% x0    optimum temperature
% 
% Jmax = a x exp(-0.5 ((x-x(0))/b)^2)
%
% values for plants grown at 21 °C:
% a  = 288  +/- 18
% b  = 18.1 +/- 2
% xo = 39.6 +/- 3.3

bunce_jmax_f = @(T,a)a*exp(-0.5*((T-39.6)/18.1).^2);
jmax_bunce_288 = bunce_jmax_f(T_range, 288);
jmax_bunce_138 = bunce_jmax_f(T_range, 138.5);

plot(T_range, A_fun(jmax_bunce_288), 'LineWidth', 2)
hold on
plot(T_range, A_fun(jmax_bunce_138), 'LineWidth', 2)

%% Farquhar 1980
% entropy [J mol-1 K-1]
S = 710;
% universal gas constant [J mol-1 K-1]
R = 8.314;
% free enthalpy [J mol-1]
H = 220000;

% J_max = 210;
J_max_25 = 138.5;

E_a.J_max = 37000;

farquhar_jmax_f = @(T)J_max_25*exp(((T-298.15)./298.15).*(E_a.J_max)./(R*T)) ...
    .* ((1+exp((298.15*S-H)./298.15*R)) ...
    ./ (1+exp((S*T-H)./(R*T))));

jmax_farquhar_138 = farquhar_jmax_f(celsius2kelvin(T_range));

plot(T_range, A_fun(jmax_farquhar_138), 'LineWidth', 2)

%% Leuning 2002
% modified Arrhenius function
% H_a       activation energy
% H_d       deactivation energy
% S_v       entropy term
% T_0       reference temperature (298.15 K)
% 
% Jmax = C*exp((H_a/(R*T_0))(1-T_0/T_1)) / (1 + exp((S_v*T_1-H_d)/(R*T_1)))
% C = 1 + exp((S_v*T_0 - H_d)/(R*T_0))
% 
% Average values
% H_a = 50300 J/mol
% H_d = 152044 J/mol
% S_v = 495 J/mol/K

H_a = 50300;
H_d = 152044;
S_v = 495;
T_0 = 298.2;
R = 8.314;

leuning_jmax_f = @(T)(1+exp((S_v*T_0-H_d)/R./T_0)) .* ...
    exp((H_a/R/T_0)*(1-T_0./T)) ./ (1+exp((S_v*T-H_d)/R./T));
jmax_leuning_av = leuning_jmax_f(celsius2kelvin(T_range));

% values are scaled to the maximum Jmax, so all values must be multiplied
% by the maxium Jmax
J_max = 138.5;
jmax_leuning_138 = jmax_leuning_av*J_max;

plot(T_range, A_fun(jmax_leuning_138), 'LineWidth', 2)

%% Ali 2015
% Four different models
% Model 1
% Jmax = Jmax(25) * ((0.8 * 2.4)^(0.1*(T_1-T_0))) * (1 + exp((S_v*T_1-H_d)/(R*T_1)))
% 
% They used the same values as in Farquhar et al. (1980).
% S_v       710 J/mol/K
% H_d       220000 J/mol
% T_0       not explicitly given but it will be 298.15 K

S_v = 710;
R = 8.314;
H_d = 220000;
T_0 = 298.15;

ali_jmax_model1_138_f = @(T)138.5*((0.8*2.4).^(0.1.*(T-T_0))) .* ...
    (1 + exp((S_v*T-H_d)/R/T));
jmax_ali_model1_138 = ali_jmax_model1_138_f(celsius2kelvin(T_range));

plot(T_range, A_fun(jmax_ali_model1_138), 'LineWidth', 2)

%% improve figure
legend(...
    {'Bunce^*', 'Bunce', 'Farquhar', 'Leuning', 'Ali'},...
    'box', 'off',...
    'NumColumns', 2,...
    'Location', 'n')
xlabel('Temperature (°C)')
ylabel('A_j (\mumol m^{-2} s^{-1})')
y_limits = get(gca, 'YLim');
set(gca,...
    'FontSize', 12,...
    'FontName', 'Arial',...
    'box', 'off',...
    'LineWidth', 1.3,...
    'YLim', [y_limits(1) 1.3*y_limits(2)])
set(gcf, 'OuterPosition', [-1201 247 348 364])

exportgraphics(gcf, 'j_max_models.png')