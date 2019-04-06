function [capacidad] = get_capacity(potencia, h, P_Noise)

capacidad = log2(1+(potencia * abs(h)^2)/P_Noise);