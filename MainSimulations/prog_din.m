function [ MSG_OUT0, MSG_OUT1 ] = prog_din(MSGS_IN0,MSGS_IN1,Nbyz)
nstep = length(MSGS_IN0);
Metriche = zeros(nstep+1,1)-1e10;
Metriche(1) = 0;
Metriche_old = Metriche;
for g = 1:nstep
    Metriche_lin = exp(Metriche_old);
    Metriche_lin = Metriche_lin*exp(MSGS_IN1(g));
    Metriche_lin(2:g+1) = Metriche_lin(2:g+1) + exp(Metriche_old(1:g))*exp(MSGS_IN0(g));
    Metriche = log(Metriche_lin);
    Metriche_old = Metriche;
end
Metriche_lin = exp(Metriche);
MSG_OUT0_pre = Metriche_lin(Nbyz);
MSG_OUT1_pre = Metriche_lin(Nbyz+1);
MSG_OUT0 = log(MSG_OUT0_pre/(MSG_OUT0_pre+MSG_OUT1_pre));
MSG_OUT1 = log(MSG_OUT1_pre/(MSG_OUT0_pre+MSG_OUT1_pre));