function [ R, s ] = crea_reports_Markov( n, m, alpha, Pmal, eps, rho, s1)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% n = numero di osservazioni temporali
% m = numero di nodi
% alpha = percentuale di bizantini
% Pmal = Probabilita' di flipping
% eps = errore di misura
% rho = probabilita' del modello
% s1 = stato uniziale
R = zeros(m,n);
s = zeros(1,n);
if s1 == -1
    s1 = randi(2)-1;
end
s(1) = s1;
for k = 2:n
    cs = rand;
    if cs < rho
        s(k) = ~s(k-1);
    else
        s(k) = s(k-1);
    end;    
end
Noise_dec = rand(m,n);
indx_flip = find( Noise_dec < eps );
U = repmat(s,m,1);
U(indx_flip) = ~U(indx_flip);
R = U;
%Per comodita' i bizantini sono i primi

Num_B = round(alpha*m);
R1 = R(1:Num_B,:);
Noise_flip = rand(Num_B,n);
indx_flip = find( Noise_flip < Pmal );
R1(indx_flip) = ~R1(indx_flip);
R(1:Num_B,:) = R1;




