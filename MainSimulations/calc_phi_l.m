function [ phi_l ] = calc_phi_l( rho, tau_l )
phi_l = 0*tau_l;
phi_l(end) = 0.5;
phi_l(1:end-1) = rho*tau_l(2:end)+(1-rho)*(1-tau_l(2:end));
end
