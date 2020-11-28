function [ omega_d ] = calc_omega_d( lambda_d )

m = size(lambda_d,2); %Numero osservazioni
n = size(lambda_d,1); %Numero nodi
omega_d = zeros(n,1);
A = log(lambda_d+1e-20);
A1 = log(1-lambda_d+1e-20);
omega_d = exp(sum(A,2))./(exp(sum(A,2))+exp(sum(A1,2)));