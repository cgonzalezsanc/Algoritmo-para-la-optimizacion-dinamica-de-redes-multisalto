%This function generate a Rayleigh Channel
%==========================================
%t -> Duration of the channel estimation
%it -> Sample time
%Sigma -> Channel Variation
%d -> Distance


function [H] = create_wless_chan (t, it, Sigma, d)

if (mod(t, it) ~= 0)
    H = [];
    return;
end;

G = 10;
    
Points = t/it;
H = (randn(1, Points)+ 1i.* randn(1, Points)) * sqrt(G*Sigma/2) * sqrt(1/d^2);

end