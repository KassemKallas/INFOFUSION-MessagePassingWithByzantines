function [ BM, ER, MESS_R0, MESS_R1 ] = MP_backward_log_def( MESS_U0, MESS_U1, s, eps, Pmal, rho)
% R = tutti i report
% s = tutti gli stati
% eps = errore di misura
% Pmal = Probabilita' di flipping
% MESS_U0, MESS_U1 = messaggi relativi allo stato dei nodi (0 = bizantini, 1 =
% onesti) provenienti dai nodi funzione 
% rho = probabilita' del modello
m = size(MESS_U0,1);
n = length(s);
BM = zeros(2,n);

for it = n:-1:1
    if it == n
        MIN_R = log([0.5 0.5]).';
    else
        MIN_P = exp(BM(:,it+1));
        MIN_R = zeros(2,1);
        MIN_R(1) = MIN_P(1)*(1-rho)+MIN_P(2)*rho;
        MIN_R(2) = MIN_P(2)*(1-rho)+MIN_P(1)*rho;  
        MIN_R_pre = log(MIN_R+1e-10);
        MIN_R(1) = log(exp(MIN_R_pre(1))./(exp(MIN_R_pre(1))+exp(MIN_R_pre(2)))+1e-10);
        MIN_R(2) = log(exp(MIN_R_pre(2))./(exp(MIN_R_pre(1))+exp(MIN_R_pre(2)))+1e-10);
    end;
    MESS_R0(it) = MIN_R(1);
    MESS_R1(it) = MIN_R(2);

    MIN_U = zeros(2,1);
    MIN_U_pre(1) = sum(MESS_U0(:,it));
    MIN_U_pre(2) = sum(MESS_U1(:,it));
    MIN_U(1) = log(exp(MIN_U_pre(1))./(exp(MIN_U_pre(1))+exp(MIN_U_pre(2)))+1e-10);
    MIN_U(2) = log(exp(MIN_U_pre(2))./(exp(MIN_U_pre(1))+exp(MIN_U_pre(2)))+1e-10);
    BM_pre = MIN_U+MIN_R;
    BM(1,it) = log(exp(BM_pre(1))./(exp(BM_pre(1))+exp(BM_pre(2)))+1e-10);
    BM(2,it) = log(exp(BM_pre(2))./(exp(BM_pre(1))+exp(BM_pre(2)))+1e-10);
end;
%BM = BM./repmat(sum(BM),2,1);
[maxv indxm] = max(BM);
EST = indxm-1;
ER = sum(xor(EST,s))/n;
end
