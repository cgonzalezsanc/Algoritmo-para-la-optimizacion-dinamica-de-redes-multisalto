function nrg_cost_n = get_energy_cost(total_duration, nrg_cost)
% Creates a matrix which contains the cost of the energy supplied by the
% operator over time.
% Cost units are cents/kW
% Values can fluctuate between 2.299 and 3.576 cents/kW for example.

noise_err = 0.05;

noise = noise_err*nrg_cost*(rand(1,total_duration)-0.5);

nrg_cost_n = nrg_cost*ones(1,total_duration)+noise;

end