function [ pri_si ] = calcola_pr_s_nu( rj, si, prob_nodi, eps, Pmal )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% rj = vettore dei report 
% si = stato al tempo i
% prob_nodi = probabilita' che lo stato dei nodi sia byzantino
% eps = errore di misura
% Pmal = Probabilita' di flipping
% alpha = percentuale di bizantini

m = length(rj);
eta = eps*(1-Pmal)+(1-eps)*Pmal;
deltas = zeros(m,1);
deltas(rj == si) = 1;
pri_si = ((1-eps)*deltas+eps*(~deltas)).*(1-prob_nodi)+((1-eta)*deltas+eta*(~deltas)).*prob_nodi;
end

