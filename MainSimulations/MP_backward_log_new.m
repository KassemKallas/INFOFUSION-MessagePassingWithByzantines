function [ BM, ER ] = MP_backward_log_new( R, s, eps, Pmal, q, rho)
% R = tutti i report
% s = tutti gli stati
% eps = errore di misura
% Pmal = Probabilita' di flipping
% q = messaggi relativi allo stato dei nodi (q(1) = bizantini, q(2) =
% onesti)
% rho = probabilita' del modello
m = size(R,1);
n = length(s);
BM = zeros(2,n);

for it = n:-1:1
    if it == n
        MIN_R = [0 0].';
    else
        MIN_P = exp(BM(:,it+1));
        MIN_R = zeros(2,1);
        MIN_R(1) = MIN_P(1)*(1-rho)+MIN_P(2)*rho;
        MIN_R(2) = MIN_P(2)*(1-rho)+MIN_P(1)*rho;  
        MIN_R = log(MIN_R+1e-10);
    end;
    MIN_U = zeros(2,1);
    
    [ pri_si ] = calcola_pr_s_nu( R(:,it), 0, (q(1,:)./(q(1,:)+q(2,:))).', eps, Pmal );
    MIN_U(1) = sum(log(pri_si+1e-10));
    [ pri_si ] = calcola_pr_s_nu( R(:,it), 1, (q(1,:)./(q(1,:)+q(2,:))).', eps, Pmal );
    MIN_U(2) = sum(log(pri_si+1e-10));
    BM(:,it) = MIN_U+MIN_R;
end;
%BM = BM./repmat(sum(BM),2,1);
[maxv indxm] = max(BM);
EST = indxm-1;
ER = sum(xor(EST,s))/n;
end
