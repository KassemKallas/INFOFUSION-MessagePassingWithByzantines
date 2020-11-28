clear all;
close all;

rng('default');
rng(12);
Niteraz = 10000;




Pmal = 1.0;
Pmal_eff = Pmal;
rho = 0.5; %Modello correlazione: 1-rho = probabilita' che la nuova osservazione sia uguale dalla precedente, when 0.5 no markov model
eps = 0.15; %Rumore di osservazione
% alpha = 0.45;
% alpha_eff = alpha;
s1 = -1;%-1 = unknown : stato iniziale
m = 30; %numero osservazioni
n = 20; %numero nodi
Num_iteraz_turbo = 5; %iterations of MP algorithm
ER = zeros(Num_iteraz_turbo,Niteraz);
ER_MAJ = zeros(1,Niteraz);
ER_MAJ_r = zeros(1,Niteraz);
ER_OTT_FIXED = zeros(1,Niteraz);
ER_OTT_STAT = zeros(1,Niteraz);
LLR_ER = zeros(1,Niteraz);
VARSH_ER = zeros(1,Niteraz);
%Nbyz = round(alpha_eff*n);
alpha_rw = 1; % Per il reweighted MP

%LLR param
% L = n/2;
% Nsoglie_LLR =  6;
% PH1=0.5; %stationary prob pi = pi*P
% gammas = 0:m;
% %end LLR param
%-------------------------------------------------------------------------
sim_param.PH1 = 0.5;
sim_param.T = m;
sim_param.N = n;
% sim_param.K1 = 9;
% sim_param.alfa = (sim_param.N-sim_param.K1)/sim_param.N;
sim_param.epsilon = eps;
sim_param.Nprove = 1;
sim_param.delta= 1 - sim_param.epsilon;
%sim_param.possible_system_states =  dec2bin(0:2^sim_param.T -1,sim_param.T)-'0';
sim_param.L = sim_param.N/2;
sim_param.gammas = 0:sim_param.T;
sim_param.Pd_Hp = 1-sim_param.epsilon;%detection at honest
sim_param.Pfa_Hp = sim_param.epsilon;%fa at honest
sim_param.Pd_Bp = 1-sim_param.epsilon;%detection at Byzantine
sim_param.Pfa_Bp = sim_param.epsilon;%false alarm at byzantines
sim_param.Nsoglie_LLR = 10;
sim_param.Pmal = Pmal;


%=====================Initialize Varshney and LLR=========================



% for epss = 1:2

i=0;
for ii= 6:(n/2)-1
i=i+1;
alpha = ii/n;
alpha_eff = alpha;
Nbyz = round(alpha_eff*n);
alf(i) = alpha;

sim_param.alfa = alpha;
%CDC param
sim_param.K1 = n-Nbyz;
Nhr = zeros(1,length(sim_param.gammas))';
Nbr = Nhr;
Nhr_llr = zeros(1,sim_param.Nsoglie_LLR)';
Nbr_llr = Nhr_llr;

fprintf('Alpha = %f\n',alpha);
for IT = 1:Niteraz %nb of trials
    tau_l = 0.5*ones(m,1);
    tau_r = 0.5*ones(m,1);
    phi_l = 0.5*ones(m,1);
    phi_r = 0.5*ones(m,1);
    nu_u = 0.5*ones(m,n);
    nu_d = 0.5*ones(m,n);
    lambda_d = 0.5*ones(n,m);
    lambda_u = 0.5*ones(n,m);
    omega_d = 0.5*ones(n,1);
    omega_u = 0.5*ones(n,1);
    
    % Previous one
    
    tau_l_p = tau_l;
    tau_r_p = tau_r;
    phi_l_p = phi_l;
    phi_r_p = phi_r;
    nu_u_p = nu_u;
    nu_d_p = nu_d;
    lambda_d_p = lambda_d;
    lambda_u_p = lambda_u;
    omega_d_p = omega_d;
    omega_u_p = omega_u;
    
    MESS = zeros(2,m);
    if rem(IT,50) == 0
        fprintf('Iterazione %d su %d\n',IT,Niteraz);
    end;
%     if epss == 1
[ R, s ] = crea_reports_Markov( m, n, alpha_eff, Pmal_eff, eps, rho, s1); % this generate the seq if the nb of nodes is fixed according to alpha
    %R = ALL_GEN(IT).R;
    %s = ALL_GEN(IT).s;
%     if IT==1
%     fprintf('Attack Type = %f\n',epss);
%     end
%     else
%[ R, s, Nbyz ] = crea_reports_Markov_stat( m, n, alpha, Pmal, eps, rho, s1);
 sim_param.K1 = n - Nbyz;
 sim_param.R = R;
 sim_param.P = s;
%  if IT==1
%    fprintf('Attack Type = %f\n',epss);
%  end
    %for statistical attack
    % Iterazioni turbo
    % Inizializzazione omega 
%     end
    omega_u = alpha;
    for tit = 1:Num_iteraz_turbo
        indx_ch = 1:m;
        [ lambda_u ] = calc_lambda_u( omega_u, lambda_d, lambda_u_p, indx_ch );
        %lambda_u = alpha_rw*lambda_u + (1-alpha_rw)*lambda_d_p;
        if tit > 1
%             lambda_u0 = exp(alpha_rw*log(lambda_u+1e-20) + (1-alpha_rw)*log(lambda_u_p+1e-20));
%             lambda_u1 = exp(alpha_rw*log(1-lambda_u+1e-20) + (1-alpha_rw)*log(1-lambda_u_p+1e-20));
%             lambda_u = lambda_u0./(lambda_u0+lambda_u1);
            
            lambda_u0 = alpha_rw*lambda_u + (1-alpha_rw)*lambda_d_p;
            lambda_u1 = alpha_rw*(1-lambda_u) + (1-alpha_rw)*(1-lambda_d_p);
            lambda_u = lambda_u0./(lambda_u0+lambda_u1);

%             lambda_u0 = exp(alpha_rw*log(lambda_u+1e-20) + (1-alpha_rw)*log(lambda_d_p+1e-20));
%             lambda_u1 = exp(alpha_rw*log(1-lambda_u+1e-20) + (1-alpha_rw)*log(1-lambda_d_p+1e-20));
%             lambda_u = lambda_u0./(lambda_u0+lambda_u1);

        end;
        [ nu_u ] = calc_nu_u( lambda_u, R, eps, Pmal, nu_u_p, indx_ch );
        %nu_u = alpha_rw*nu_u + (1-alpha_rw)*nu_d_p;
%         nu_u0 = exp(alpha_rw*log(nu_u+1e-20) + (1-alpha_rw)*log(nu_d_p+1e-20));
%         nu_u1 = exp(alpha_rw*log(1-nu_u+1e-20) + (1-alpha_rw)*log(1-nu_d_p+1e-20));
%         nu_u = nu_u0./(nu_u0+nu_u1);
        
        nu_u0 = alpha_rw*nu_u + (1-alpha_rw)*nu_d_p;
        nu_u1 = alpha_rw*(1-nu_u) + (1-alpha_rw)*(1-nu_d_p);
        nu_u = nu_u0./(nu_u0+nu_u1);

%         nu_u0 = exp(alpha_rw*log(nu_u) + (1-alpha_rw)*log(nu_u_p));
%         nu_u1 = exp(alpha_rw*log(1-nu_u) + (1-alpha_rw)*log(1-nu_u_p));
%         nu_u = nu_u0./(nu_u0+nu_u1);
        [ tau_r ] = calc_tau_r( phi_r, nu_u );
        %tau_r = alpha_rw*tau_r + (1-alpha_rw)*phi_l_p;
        %tau_r = exp(alpha_rw*log(tau_r) + (1-alpha_rw)*log(tau_r_p));
        [ tau_l ] = calc_tau_l( phi_l, nu_u );
        %tau_l = alpha_rw*tau_l + (1-alpha_rw)*phi_r_p;
        %tau_l = exp(alpha_rw*log(tau_l) + (1-alpha_rw)*log(tau_l_p));
        [ phi_r ] = calc_phi_r( 1-rho, tau_r, 0.5 );
        %phi_r = alpha_rw*phi_r + (1-alpha_rw)*tau_l_p;
        %phi_r = exp(alpha_rw*log(phi_r) + (1-alpha_rw)*log(phi_r_p));
        [ phi_l ] = calc_phi_l( 1-rho, tau_l );
        %phi_l = alpha_rw*phi_l + (1-alpha_rw)*tau_l_p;
        %phi_l = exp(alpha_rw*log(phi_l) + (1-alpha_rw)*log(phi_l_p));
        [ nu_d ] = calc_nu_d( phi_r, phi_l, nu_u );
        %nu_d = alpha_rw*nu_d + (1-alpha_rw)*nu_u_p;
%         nu_d0 = exp(alpha_rw*log(nu_d) + (1-alpha_rw)*log(nu_d_p));
%         nu_d1 = exp(alpha_rw*log(1-nu_d) + (1-alpha_rw)*log(1-nu_d_p));
%         nu_d = nu_d0./(nu_d0+nu_d1);
        
        nu_d0 = alpha_rw*nu_d + (1-alpha_rw)*nu_u_p;
        nu_d1 = alpha_rw*(1-nu_d) + (1-alpha_rw)*(1-nu_u_p);
        nu_d = nu_d0./(nu_d0+nu_d1);
        
        [ lambda_d ] = calc_lambda_d( nu_d, R, eps, Pmal, lambda_d_p, indx_ch );
%         nu_d0 = exp(alpha_rw*log(nu_d+1e-20) + (1-alpha_rw)*log(nu_u_p+1e-20));
%         nu_d1 = exp(alpha_rw*log(1-nu_d+1e-20) + (1-alpha_rw)*log(1-nu_u_p+1e-20));
%         nu_d = nu_d0./(nu_d0+nu_d1);
%         [ lambda_d ] = calc_lambda_d( nu_d, R, eps, Pmal );
        %lambda_d = alpha_rw*lambda_d + (1-alpha_rw)*lambda_u_p;
%         lambda_d0 = exp(alpha_rw*log(lambda_d) + (1-alpha_rw)*log(lambda_d_p));
%         lambda_d1 = exp(alpha_rw*log(1-lambda_d) + (1-alpha_rw)*log(1-lambda_d_p));
%         lambda_d = lambda_d0./(lambda_d0+lambda_d1);
        lambda_d0 = alpha_rw*lambda_d + (1-alpha_rw)*lambda_u_p;
        lambda_d1 = alpha_rw*(1-lambda_d) + (1-alpha_rw)*(1-lambda_u_p);
        lambda_d = lambda_d0./(lambda_d0+lambda_d1);

        
%         lambda_d0 = exp(alpha_rw*log(lambda_d+1e-20) + (1-alpha_rw)*log(lambda_u_p+1e-20));
%         lambda_d1 = exp(alpha_rw*log(1-lambda_d+1e-20) + (1-alpha_rw)*log(1-lambda_u_p+1e-20));
%         lambda_d = lambda_d0./(lambda_d0+lambda_d1);

        [ omega_d ] = calc_omega_d( lambda_d );
        %omega_d = exp(alpha_rw*log(omega_d) + (1-alpha_rw)*log(omega_d_p));
        
        MESS(1,:) = exp(sum(log(nu_u+1e-20),2)+log(phi_r+1e-20)+log(phi_l+1e-20));
        MESS(2,:) = exp(sum(log(1-nu_u+1e-20),2)+log(1-phi_r+1e-20)+log(1-phi_l+1e-20));
        
        PS = MESS./repmat(sum(MESS),2,1);
        PS_ALL(tit,IT).vals = PS;
        %PS_ALL2(tit,IT).val = PS;
        [maxv indxm] = max(MESS);
        EST1 = indxm-1;
        indx = find(abs(PS(1,:)-0.5) < 1e-8);
        EST1(indx) = 0;

        ER(tit,IT) = sum(xor(EST1,s))/m;

        tau_l_p = tau_l;
        tau_r_p = tau_r;
        phi_l_p = phi_l;
        phi_r_p = phi_r;
        nu_u_p = nu_u;
        nu_d_p = nu_d;
        lambda_d_p = lambda_d;
        lambda_u_p = lambda_u;
        omega_d_p = omega_d;
        omega_u_p = omega_u;
    end;
    
    % Verifica se e' necessaria inversione, since in some cases the program
    % was giving high error probability and since it's suboptimal over a
    % tanner graph (not tree) so it's convergence can't be computed
    % precisely and then this to account for the cases where measures
    % become so near and may cause decision inversion.
    EST2 = ~EST1;
    delta = eps*(1-Pmal)+(1-eps)*Pmal;
    MAT_DIFF1 = xor(R,repmat(EST1,n,1));
    metrica1 = sum(log((1-alpha)*prod(((1-eps).^(~MAT_DIFF1)).*(eps.^(MAT_DIFF1)),2)+alpha*prod(((1-delta).^(~MAT_DIFF1)).*(delta.^(MAT_DIFF1)),2)));
    MAT_DIFF2 = xor(R,repmat(EST2,n,1));
    metrica2 = sum(log((1-alpha)*prod(((1-eps).^(~MAT_DIFF2)).*(eps.^(MAT_DIFF2)),2)+alpha*prod(((1-delta).^(~MAT_DIFF2)).*(delta.^(MAT_DIFF2)),2)));
    if metrica1 > metrica2
        EST = EST1;
    else
        EST = EST2;
    end;
    ER(Num_iteraz_turbo,IT) = sum(xor(EST,s))/m;

    %Calcolo con decisione a maggioranza
    MAJ_C = sum(R);
    EST_MAJ = zeros(1,m);
    EST_MAJ(MAJ_C >= n/2) = 1;
    ER_MAJ(IT) = sum(xor(EST_MAJ,s))/m;
    
    %Calcolo con decisione a maggioranza con rimozione %VARSHNEY SCHEME
    REP = zeros(n,1);
    for g = 1:n
        REP(g) = sum(xor(R(g,:),EST_MAJ));
    end;
    [sortv indxs] = sort(REP);
    indx_honest = indxs(1:n-Nbyz);
    %indx_honest = Nbyz+1:m;
    Rh = R(indx_honest,:);
    MAJ_C = sum(Rh);
    EST_MAJ = zeros(1,m);
   EST_MAJ(MAJ_C > length(indx_honest)/2) = 1;
   %EST_MAJ(MAJ_C > n/2) = 1;
    ER_MAJ_r(IT) = sum(xor(EST_MAJ,s))/m;
    
    %Calcolo con decisioni ottime %THIS IS ONLY FOR INDEPENDENT STATE ;
    %DOES NOT WORK WHEN N is FIXED
    ER_OTT_STAT(IT) = -1; 
    if m <= 10
        s_est_ott_fixed = 0*s;
        s_est_ott_stat = 0*s;
        max_metrica_ott_fixed = -1e10;
        max_metrica_ott_stat = -1e10;
        metriche_stat = zeros(1,2^m);
        for mc = 0:2^m-1
            s_hat = dec2bin(mc,m)-48;
            diff_s = xor(s_hat(2:end),s_hat(1:end-1));
            pS = prod((1-rho).^(~diff_s).*rho.^(diff_s));
            % Ottimo statistico
            delta = eps*(1-Pmal)+(1-eps)*Pmal;
            MAT_DIFF = xor(R,repmat(s_hat,n,1));
           metrica_ott_stat = log(pS)+sum(log((1-alpha)*prod(((1-eps).^(~MAT_DIFF)).*(eps.^(MAT_DIFF)),2)+alpha*prod(((1-delta).^(~MAT_DIFF)).*(delta.^(MAT_DIFF)),2)));
        %This is right
           %  metrica_ott_stat = sum(log((1-alpha)*prod(((1-eps).^(~MAT_DIFF)).*(eps.^(MAT_DIFF)),2)+alpha*prod(((1-delta).^(~MAT_DIFF)).*(delta.^(MAT_DIFF)),2))); 
%This is wrong           
           metriche_stat(mc+1) = metrica_ott_stat;
            if metrica_ott_stat > max_metrica_ott_stat
                s_est_ott_stat = s_hat;
                max_metrica_ott_stat = metrica_ott_stat;
            end;
        end;
        ER_OTT_STAT(IT) = sum(xor(s_est_ott_stat,s))/m;
    end;
    
    if ER_OTT_STAT(IT) == 0 && ER(Num_iteraz_turbo,IT) > 0
        check = 1;
    end;
    
 [results] = CDC(sim_param); 
%============= Results Varshney and LLR====================================
Nhr = Nhr + results.Nerr_hr;
Nbr = Nbr + results.Nerr_br;

Nhr_llr = Nhr_llr + results.Nerr_hr_LLR;
Nbr_llr = Nbr_llr + results.Nerr_br_LLR;   
%====================================END RES==============================    
%     if rem(IT,100) == 0
%         fprintf('Error rate con Turbo decoding = %f\n Error rate con rimozione hard = %f\nError rate rx ottimo modello statistico = %f\n',mean(ER(Num_iteraz_turbo,1:IT)),mean(ER_MAJ_r(1:IT)),mean(ER_OTT_STAT(1:IT)));
%     end;    
end;
for ig = 1:length(sim_param.gammas)
    PDr(ig) = 1-Nhr(ig)/Niteraz/m;
    PFAr(ig) = Nbr(ig)/Niteraz/m;
end;
for is = 1:sim_param.Nsoglie_LLR
    PDr_LLR(is) = 1-Nhr_llr(is)/Niteraz/m;
    PFAr_LLR(is) = Nbr_llr(is)/Niteraz/m;
end;
Varsh(i) = min(PFAr + 1-PDr);
LLR(i) = min(PFAr_LLR+1-PDr_LLR);



fprintf('Errr rate con decisione a maggioranza = %f\n',mean(ER_MAJ));
fprintf('Errr rate con decisione a maggioranza dopo la rimozione = %f\n',mean(ER_MAJ_r));
fprintf('Error rate rx ottimo modello statistico = %f\n',mean(ER_OTT_STAT));
fprintf('Error rate Varshney = %f\n',Varsh(i));
fprintf('Error rate LLR = %f\n',LLR(i));

for l = 1:Num_iteraz_turbo
    fprintf('Error rate con Turbo decoding dopo %d iterazioni = %f\n',l,mean(ER(l,:)));
end;

err_maj(i) = mean(ER_MAJ);
err_maj_r(i) = mean(ER_MAJ_r);
err_opt(i) = mean(ER_OTT_STAT);
err_mp(i) = mean(ER(end,:));
err_varsh(i) = Varsh(i);
err_llr(i) = LLR(i);
end


% %fprintf('Errr rate con algoritmo BCJR prima della rimozione = %f\n',mean(ER));
% %fprintf('Errr rate con algoritmo BCJR dopo la rimozione = %f\n',mean(ER_r));
figure(1);
set(gca,'fontsize',20);
hold on;
plot(alf,err_maj,'blue-+','LineWidth',2,'markers',8);
set(gca,'yscale','log');
hold on;
plot(alf,err_varsh,'green-o','LineWidth',2,'markers',8);
set(gca,'yscale','log');
hold on;
plot(alf,err_llr,'cyan-s','LineWidth',2,'markers',8);
set(gca,'yscale','log');
hold on;
plot(alf,err_mp,'black-*','LineWidth',2,'markers',8);
set(gca,'yscale','log');
if m<=10
hold on;
plot(alf,err_opt,'red-x','LineWidth',2,'markers',8);
set(gca,'yscale','log');
end
grid on;
if m<=10
legend('Majority','Hard Isolation Scheme','Soft Isolation Scheme','Message Passing','Optimal','Location','northwest')
else
legend('Majority','Hard Isolation Scheme','Soft Isolation Scheme','Message Passing','Location','northwest')
end
xlabel('\alpha');
ylabel('log(P_e)');
xlabel('\alpha','FontSize',20);
ylabel('log(P_e)','FontSize',20);
%saveas(gcf,'Pe.fig');
