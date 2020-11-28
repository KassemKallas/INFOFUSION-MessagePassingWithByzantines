function [ lambda_u ] = calc_lambda_u( omega_u, lambda_d, lambda_u_old, indx_ch )
lambda_u = 0*lambda_d;
m = size(lambda_d,2); %Numero osservazioni
n = size(lambda_d,1); %Numero nodi
A = log(omega_u+1e-20);
A1 = log(1-omega_u+1e-20);
if indx_ch < 0
    indx_ch = 1:m;
end;
for i = indx_ch
    lambda_log_pre1 = log(lambda_d+1e-20);
    lambda_log_pre1(:,i) = 0;
    lambda_log_pre1 = sum(lambda_log_pre1,2);
    lambda_log_pre2 = log(1-lambda_d+1e-20);
    lambda_log_pre2(:,i) = 0;
    lambda_log_pre2 = sum(lambda_log_pre2,2);
    C1 = exp(A + lambda_log_pre1);
    C2 = exp(A1 + lambda_log_pre2);
    lambda_u(:,i) = (C1./(C1+C2));
end;
