function g_s = adj_g_s(T)
%% g_s = adj_g_s(T)
% Estimate the stomatal conductance to CO2 at a given temperature by
% fitting experimental data from von Caemmerer and Evans 2015 (DOI:10.1111/pce.12449)
% to a quadratic function:
%               g_s(T) = a*T^2 + b*T + c.
% Input:
%       double T:               temperature in Kelvin
% Output:
%       double g_s:             adjusted stomatal conductance
warning off
% experimental data
T_data = celsius2kelvin([15 20 25 30 35 40]);
g_s_data = [0.218 0.251 0.273 0.283 0.274 0.303];

% fit to quadratic function
f = @(b,T) b(1)*T.^2 + b(2)*T + b(3);

% non-linear regression (min SSE)
beta = nlinfit(T_data, g_s_data, f, [1 1 1]);

% estimate stomatal conductance at the given temperature
g_s = f(beta,T);
warning on
end

