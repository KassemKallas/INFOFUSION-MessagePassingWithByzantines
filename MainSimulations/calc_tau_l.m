function [ tau_l ] = calc_tau_l( phi_l, nu_u )
tau_l = 0*phi_l;
A = exp(sum(log(nu_u+1e-20),2));
B = exp(sum(log(1-nu_u+1e-20),2));
tau_l = phi_l.*A./(phi_l.*A+(1-phi_l).*B);
end

