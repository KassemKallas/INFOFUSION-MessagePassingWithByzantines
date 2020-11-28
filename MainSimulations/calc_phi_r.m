function [ phi_r ] = calc_phi_r( rho, tau_r, pS10 )
phi_r = 0*tau_r;
phi_r(1) = pS10;
phi_r(2:end) = rho*tau_r(1:end-1)+(1-rho)*(1-tau_r(1:end-1));
end
