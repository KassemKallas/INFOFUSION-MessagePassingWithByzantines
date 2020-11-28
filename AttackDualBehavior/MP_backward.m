function [ BM, ER ] = MP_backward( R, s, eps, Pmal, alpha, rho)
% R = tutti i report
% s = tutti gli stati
% eps = errore di misura
% Pmal = Probabilita' di flipping
% alpha = percentuale di bizantini
% rho = probabilita' del modello
m = size(R,1);
n = length(s);
BM = zeros(2,n);

for it = n:-1:1
    if it == n
        MIN_R = [1 1].';
    else
        MIN_P = BM(:,it+1);
        MIN_R = zeros(2,1);
        MIN_R(1) = MIN_P(1)*(1-rho)+MIN_P(2)*rho;
        MIN_R(2) = MIN_P(2)*(1-rho)+MIN_P(1)*rho;        
    end;
    MIN_U = zeros(2,1);
    [ pri_si ] = calcola_pr_s( R(:,it), 0, eps, Pmal, alpha );
    MIN_U(1) = prod(pri_si)*2^m;
    [ pri_si ] = calcola_pr_s( R(:,it), 1, eps, Pmal, alpha );
    MIN_U(2) = prod(pri_si)*2^m;
    BM(:,it) = MIN_U.*MIN_R;
end;
BM = BM./repmat(sum(BM),2,1);
[maxv indxm] = max(BM);
EST = indxm-1;
ER = sum(xor(EST,s))/n;
end
