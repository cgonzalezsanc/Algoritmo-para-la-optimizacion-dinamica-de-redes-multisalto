function p_wind = get_wind_energy(n_nodes, t, wind_param)
% Calculates the power extracted from wind. Wind speed follows a Weibull
% distribution.

k = wind_param.k;   % Shape parameter
l = wind_param.l;   % Scale parameter
p = wind_param.p;   % Air density
d = wind_param.d;   % Rotor blade diameter
Cp = wind_param.Cp; % Coefficient of performance
N = wind_param.N;   % Machinery efficiency

% Generate the distribution of wind speed
v_wind = wblrnd(l,k,n_nodes,t);

% Calculate the energy harvested from wind
p_wind = 0.5* p * pi*(d/2)^2 * Cp * v_wind.^3 * N; % -> habria que dividirlo por el tiempo de cada instante

end

