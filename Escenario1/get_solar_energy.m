function p_solar = get_solar_energy(n_nodes, t, solar_param)
% Calculates solar power. Solar irradiance follows a Beta distribution.

H = solar_param.H;      % Mean solar irradiation (kWh/m^2)
S = solar_param.S;      % Total solar panel area
r = solar_param.r;      % Solar panel efficiency
Pr = solar_param.Pr;    % Performance ratio

% Beta parameters
B = 2;
A = 4;
beta_mean = A/(A+B);    % mean of beta distribution

% Calculate maximum solar radiation
max_sol_rad = H/beta_mean;

% Calculates solar radiation
sol_rad = betarnd(A, B, n_nodes, t)*max_sol_rad*1000/3600; %transform into Ws

% Calculates power output
p_solar = S * r * sol_rad * Pr;

end

