function p_harv = get_energy_harvested(n_nodes, t, solar_param, wind_param, it)
% Creates a matrix which contains the energy harvested by every node over
% time. There are two energy sources: solar and wind.
% it is time interval

p_solar = get_solar_energy(n_nodes, t, solar_param);
p_wind = get_wind_energy(n_nodes, t, wind_param);

p_harv = (p_solar + p_wind)*it; % Divided by time interval


end

