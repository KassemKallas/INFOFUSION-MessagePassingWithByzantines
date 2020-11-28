function [ lambda_d ] = calc_lambda_d( nu_d, R, eps, Pmal, lambda_d_old, indx_ch )

lambda_d = lambda_d_old;
m = size(nu_d,1); %Numero osservazioni
n = size(nu_d,2); %Numero nodi
PR00 = zeros(m,n);
PR01 = zeros(m,n);
PR10 = zeros(m,n);
PR11 = zeros(m,n);
if indx_ch < 0
    indx_ch = 1:m;
end;
for i = indx_ch
    alpha = 1; %Bizantino
    PR00(i,:) = calcola_pr_s( R(:,i), 0, eps, Pmal, alpha );
    PR10(i,:) = calcola_pr_s( R(:,i), 1, eps, Pmal, alpha );
    alpha = 0;
    PR01(i,:) = calcola_pr_s( R(:,i), 0, eps, Pmal, alpha );
    PR11(i,:) = calcola_pr_s( R(:,i), 1, eps, Pmal, alpha );
end
lambda_d(:,indx_ch) = ((PR00(indx_ch,:).*nu_d(indx_ch,:)+PR10(indx_ch,:).*(1-nu_d(indx_ch,:)))./(PR00(indx_ch,:).*nu_d(indx_ch,:)+PR10(indx_ch,:).*(1-nu_d(indx_ch,:))+PR01(indx_ch,:).*nu_d(indx_ch,:)+PR11(indx_ch,:).*(1-nu_d(indx_ch,:)))).';


