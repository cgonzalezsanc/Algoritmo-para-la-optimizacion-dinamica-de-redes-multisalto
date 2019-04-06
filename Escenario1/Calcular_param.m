function [phi, pij_h, cij_h] = Calcular_param(h, ro, pi_tot, Limit, Noise)

[pij_h] = W_filling(ro, pi_tot, h, Limit);
[cij_h] = get_capacity(pij_h, h, Noise);

phi = -1 * ro*cij_h  +  pi_tot*pij_h;  