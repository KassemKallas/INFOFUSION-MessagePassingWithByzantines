function [ nu_d ] = calc_nu_d( phi_r, phi_l, nu_u )

nu_d = 0*nu_u;
m = size(nu_u,1); %Numero osservazioni
n = size(nu_u,2); %Numero nodi
A = log(phi_r+1e-20);
B = log(phi_l+1e-20);
A1 = log(1-phi_r+1e-20);
B1 = log(1-phi_l+1e-20);
for j = 1:n
    nu_log_pre1 = log(nu_u+1e-20);
    nu_log_pre1(:,j) = 0;
    nu_log_pre1 = sum(nu_log_pre1,2);
    nu_log_pre2 = log(1-nu_u+1e-20);
    nu_log_pre2(:,j) = 0;
    nu_log_pre2 = sum(nu_log_pre2,2);
    C1 = exp(A + B + nu_log_pre1);
    C2 = exp(A1 + B1 + nu_log_pre2);
    nu_d(:,j) = (C1./(C1+C2));
end;
