function [ tau_r ] = calc_tau_r( phi_r, nu_u )
tau_r = 0*phi_r;
A = exp(sum(log(nu_u+1e-20),2));
B = exp(sum(log(1-nu_u+1e-20),2));
tau_r = phi_r.*A./(phi_r.*A+(1-phi_r).*B);
end
