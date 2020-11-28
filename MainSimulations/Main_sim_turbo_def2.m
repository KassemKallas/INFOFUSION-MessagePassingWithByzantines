rng('default');
rng(12);
Niteraz = 10000;

Pmal = 1;
rho = 0.5;
eps = 0.1;
alpha = 0.45;
s1 = -1;%-1 = unknown
m = 4; %numero osservazioni
n = 20; %numero nodi
Num_iteraz_turbo = 5;
ER = zeros(Num_iteraz_turbo,Niteraz);
Nbyz = round(alpha*n);
alpha_rw = 0.9; % Per il reweighted MP

for IT = 1:Niteraz
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
    
    [ R, s ] = crea_reports_Markov( m, n, alpha, Pmal, eps, rho, s1);
    %R = ALL_GEN(IT).R;
    %s = ALL_GEN(IT).s;

    %[ R, s, Nbyz ] = crea_reports_Markov_stat( m, n, alpha, Pmal, eps, rho, s1);
    % Iterazioni turbo
    % Inizializzazione omega 
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
(m);
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
        EST = indxm-1;
        indx = find(abs(PS(1,:)-0.5) < 1e-8);
        EST(indx) = 0;

        ER(tit,IT) = sum(xor(EST,s))/m;

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
    
    %Calcolo con decisione a maggioranza
    MAJ_C = sum(R);
    EST_MAJ = zeros(1,m);
    EST_MAJ(MAJ_C > n/2) = 1;
    ER_MAJ(IT) = sum(xor(EST_MAJ,s))/m;
    
    %Calcolo con decisione a maggioranza con rimozione
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
    ER_MAJ_r(IT) = sum(xor(EST_MAJ,s))/m;
    if rem(IT,100) == 0
        fprintf('Errr rate con Turbo decoding = %f\n',mean(ER(Num_iteraz_turbo,1:IT)));
    end;    
end;
fprintf('Errr rate con decisione a maggioranza = %f\n',mean(ER_MAJ));
fprintf('Errr rate con decisione a maggioranza dopo la rimozione = %f\n',mean(ER_MAJ_r));
for l = 1:Num_iteraz_turbo
    fprintf('Errr rate con Turbo decoding dopo %d iterazioni = %f\n',l,mean(ER(l,:)));
end;

%fprintf('Errr rate con algoritmo BCJR prima della rimozione = %f\n',mean(ER));
%fprintf('Errr rate con algoritmo BCJR dopo la rimozione = %f\n',mean(ER_r));


