function [ nu_u ] = calc_nu_u( lambda_u, R, eps, Pmal, nu_u_old, indx_ch )

m = size(lambda_u,2); %Numero osservazioni
n = size(lambda_u,1); %Numero nodi
PR00 = zeros(m,n);
PR01 = zeros(m,n);
PR10 = zeros(m,n);
PR11 = zeros(m,n);
if indx_ch < 0
    indx_ch = 1:m;
end;
nu_u = nu_u_old;
for i = indx_ch
    alpha = 1; %Bizantino
    PR00(i,:) = calcola_pr_s( R(:,i), 0, eps, Pmal, alpha );
    PR10(i,:) = calcola_pr_s( R(:,i), 1, eps, Pmal, alpha );
    alpha = 0;
    PR01(i,:) = calcola_pr_s( R(:,i), 0, eps, Pmal, alpha );
    PR11(i,:) = calcola_pr_s( R(:,i), 1, eps, Pmal, alpha );
end
nu_u(indx_ch,:) = (PR00(indx_ch,:).*(lambda_u(:,indx_ch).')+PR01(indx_ch,:).*(1-lambda_u(:,indx_ch).'))./(PR00(indx_ch,:).*(lambda_u(:,indx_ch).')+PR01(indx_ch,:).*(1-lambda_u(:,indx_ch).')+PR10(indx_ch,:).*(lambda_u(:,indx_ch).')+PR11(indx_ch,:).*(1-lambda_u(:,indx_ch).'));

