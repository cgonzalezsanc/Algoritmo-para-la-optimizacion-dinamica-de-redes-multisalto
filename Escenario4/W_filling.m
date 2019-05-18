function [potencia] = W_filling(ro, pi, h, P_max)

potencia = ro/pi - 1/(abs(h)^2);

if (potencia < 0)
    potencia = 0;
end;
if (potencia > P_max)
    potencia = P_max;
end;