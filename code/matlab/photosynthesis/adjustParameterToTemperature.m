function p_adj = adjustParameterToTemperature(varargin)
%% p_adj = adjustParameterToTemperature(varargin)
% Input:
%   double p_ref:       parameter measured at reference tempterature T_ref
%   double Ea:          activation energy for temperature adjustment [J mol-1]
%   double T_dest:      destination temperature to which p should be
%                       adjusted [K]
%   double T_ref:       (optional) temperature at which p was measured;
%                       default: 298.15 K
% Output:
%   double p_adj:       adjusted parameter

% default reference temperature
T_DEFAULT = 298.15;
% universal gas constant
R = config('R');

% validation function for (biologically relevant) temperatures in Kelvin
validTemp = @(x) isnumeric(x) && all(x >= 273.15) && all(x <= 333.15);

%% Parse input
p = inputParser;

addParameter(p,'p_ref', @isnumeric)
addParameter(p,'Ea', @isnumeric)
addParameter(p,'T_dest', T_DEFAULT, validTemp)
addParameter(p,'T_ref', T_DEFAULT, validTemp)

parse(p,varargin{:});

p_ref = p.Results.p_ref;
Ea = p.Results.Ea;
T_dest = p.Results.T_dest;
T_ref = p.Results.T_ref;

%% Temperature adjustment
p_adj = p_ref*exp((T_dest-T_ref).*(Ea./(T_ref*R*T_dest)));

end
