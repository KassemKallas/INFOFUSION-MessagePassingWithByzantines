function [ FM, ER, MESS_L0, MESS_L1 ] = MP_forward_log_def( MESS_U0, MESS_U1, s, eps, Pmal, rho, s1 )

% R = tutti i report
% s = tutti gli stati
% eps = errore di misura
% Pmal = Probabilita' di flipping
% alpha = percentuale di bizantini
% q = messaggi relativi allo stato dei nodi (q(1) = bizantini, q(2) =
% onesti)
% rho = probabilita' del modello
% s1 = stato uniziale
m = size(MESS_U0,1);
n = length(s);
FM = zeros(2,n);

for it = 1:n
    if it == 1
        if s1 == -1
            MIN_L = log([0.5 0.5].');
        else
            MIN_L = log([~s1 s1].'+1e-10);
        end;
    else
        MIN_P = exp(FM(:,it-1));
        MIN_L = zeros(2,1);
        MIN_L(1) = MIN_P(1)*(1-rho)+MIN_P(2)*rho;
        MIN_L(2) = MIN_P(2)*(1-rho)+MIN_P(1)*rho; 
        MIN_L_pre = log(MIN_L+1e-10);
        MIN_L(1) = log(exp(MIN_L_pre(1))/(exp(MIN_L_pre(1))+exp(MIN_L_pre(2)))+1e-10);
        MIN_L(2) = log(exp(MIN_L_pre(2))/(exp(MIN_L_pre(1))+exp(MIN_L_pre(2)))+1e-10);
    end;
    MESS_L0(it) = MIN_L(1);
    MESS_L1(it) = MIN_L(2);

    MIN_U = zeros(2,1);
    MIN_U_pre(1) = sum(MESS_U0(:,it));
    MIN_U_pre(2) = sum(MESS_U1(:,it));
    MIN_U(1) = log(exp(MIN_U_pre(1))/(exp(MIN_U_pre(1))+exp(MIN_U_pre(2)))+1e-10);
    MIN_U(2) = log(exp(MIN_U_pre(2))/(exp(MIN_U_pre(1))+exp(MIN_U_pre(2)))+1e-10);
    FM_pre = MIN_U+MIN_L;
    FM(1,it) = log(exp(FM_pre(1))./(exp(FM_pre(1))+exp(FM_pre(2)))+1e-10);
    FM(2,it) = log(exp(FM_pre(2))./(exp(FM_pre(1))+exp(FM_pre(2)))+1e-10);
end;
%FM = FM./repmat(sum(FM),2,1);
[maxv indxm] = max(FM);
EST = indxm-1;
ER = sum(xor(EST,s))/n;
end
