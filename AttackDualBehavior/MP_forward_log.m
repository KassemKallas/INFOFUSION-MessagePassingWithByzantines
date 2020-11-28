function [ FM, ER ] = MP_forward_log( R, s, eps, Pmal, alpha, rho, s1 )
% R = tutti i report
% s = tutti gli stati
% eps = errore di misura
% Pmal = Probabilita' di flipping
% alpha = percentuale di bizantini
% rho = probabilita' del modello
% s1 = stato uniziale
m = size(R,1);
n = length(s);
FM = zeros(2,n);

for it = 1:n
    if it == 1
        MIN_L = log([~s1 s1].'+1e-10);
    else
        MIN_P = exp(FM(:,it-1));
        MIN_L = zeros(2,1);
        MIN_L(1) = MIN_P(1)*(1-rho)+MIN_P(2)*rho;
        MIN_L(2) = MIN_P(2)*(1-rho)+MIN_P(1)*rho; 
        MIN_L = log(MIN_L+1e-10);
    end;
    MIN_U = zeros(2,1);
    [ pri_si ] = calcola_pr_s( R(:,it), 0, eps, Pmal, alpha );
    MIN_U(1) = log(prod(pri_si)+1e-10);
    [ pri_si ] = calcola_pr_s( R(:,it), 1, eps, Pmal, alpha );
    MIN_U(2) = log(prod(pri_si)+1e-10);
    FM(:,it) = MIN_U+MIN_L;
end;
%FM = FM./repmat(sum(FM),2,1);
[maxv indxm] = max(FM);
EST = indxm-1;
ER = sum(xor(EST,s))/n;
end
