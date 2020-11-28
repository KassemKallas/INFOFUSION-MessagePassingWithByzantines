rng('default');
rng(55);
Niteraz = 1000;

Pmal = 1;
rho = 0.1;
eps = 0.1;
alpha = 0.45;
STAT_NODE = 1; %0 = FIXED
s1 = -1;%-1 = unknown
n = 30;
m = 20;
Num_iteraz_turbo = 10;
ER = zeros(Num_iteraz_turbo,Niteraz);
Nbyz = round(alpha*m);
alpha_rw = 0.25; % Per il reweighted MP

for IT = 1:Niteraz
    if rem(IT,50) == 0
        fprintf('Iterazione %d su %d\n',IT,Niteraz);
    end;
    %[ R, s ] = crea_reports_Markov( n, m, alpha, Pmal, eps, rho, s1);
    ALL_GEN(IT).R = R;
    ALL_GEN(IT).s = s;
    [ R, s, Nbyz ] = crea_reports_Markov_stat( n, m, alpha, Pmal, eps, rho, s1);
    % Inizializzazione dei messaggi relativi ai valori degli stati
    MESS_L0 = zeros(1,n);
    MESS_L1 = zeros(1,n);
    MESS_R0 = zeros(1,n);
    MESS_R1 = zeros(1,n);
    MESS_U0 = zeros(m,n);
    MESS_U1 = zeros(m,n);
    MESS_D0 = zeros(m,n);
    MESS_D1 = zeros(m,n);
    % Inizializzazione dei messaggi relativi allo stato dei nodi
    MESS_NU_D0 = zeros(n,m);
    MESS_NU_D1 = zeros(n,m);
    MESS_NU_U0 = zeros(n,m);
    MESS_NU_U1 = zeros(n,m);
    
    MESS_L0_prev = zeros(1,n);
    MESS_L1_prev = zeros(1,n);
    MESS_R0_prev = zeros(1,n);
    MESS_R1_prev = zeros(1,n);
    MESS_U0_prev = zeros(m,n);
    MESS_U1_prev = zeros(m,n);
    MESS_D0_prev = zeros(m,n);
    MESS_D1_prev = zeros(m,n);
    % Inizializzazione dei messaggi relativi allo stato dei nodi
    MESS_NU_D0_prev = zeros(n,m);
    MESS_NU_D1_prev = zeros(n,m);
    MESS_NU_U0_prev = zeros(n,m);
    MESS_NU_U1_prev = zeros(n,m);
 
    MESS_OM_D0 = zeros(m,1);
    MESS_OM_D1 = zeros(m,1);
    MESS_OM_U0 = zeros(m,1);
    MESS_OM_U1 = zeros(m,1);
    
    % Iterazioni turbo
    
    for tit = 1:Num_iteraz_turbo
        % Primo step: calcolo dei messaggi che arrivano ai nodi stato per
        % effetto delle osservazioni e dello stato dei bizantini
        Pbiz = exp(sum(MESS_NU_D0))./(exp(sum(MESS_NU_D0))+exp(sum(MESS_NU_D1)));
 
        if tit == 1 || STAT_NODE == 1
            MESS_OM_U0 = log(alpha)*ones(m,1);
            MESS_OM_U1 = log(1-alpha)*ones(m,1);
        else
            MESS_OM_D0_pre = sum(MESS_NU_D0);
            MESS_OM_D1_pre = sum(MESS_NU_D1);
            MESS_OM_D0 = log(exp(MESS_OM_D0_pre)./(exp(MESS_OM_D0_pre)+exp(MESS_OM_D1_pre))+1e-10);
            MESS_OM_D1 = log(exp(MESS_OM_D1_pre)./(exp(MESS_OM_D0_pre)+exp(MESS_OM_D1_pre))+1e-10);
            for q = 1:m
                [MESS_OM_U0(q),MESS_OM_U1(q)] = prog_din(MESS_OM_D0([1:q-1 q+1:end]),MESS_OM_D1([1:q-1 q+1:end]),Nbyz);
            end;
        end;    

        for g = 1:n
            if tit == 1
                MESS_NU_U0(g,:) = log(alpha);
                MESS_NU_U1(g,:) = log(1-alpha);
            else
                MESS_NU_D0_til_n = MESS_NU_D0;
                MESS_NU_D0_til_n(g,:) = 0; % Escludo il ramo verso cui devo mandare il messaggio
                MESS_NU_D1_til_n = MESS_NU_D1;
                MESS_NU_D1_til_n(g,:) = 0; % Escludo il ramo verso cui devo mandare il messaggio
                MESS_NU_U0_pre = sum(MESS_NU_D0_til_n)+MESS_OM_U0.';
                MESS_NU_U1_pre = sum(MESS_NU_D1_til_n)+MESS_OM_U1.';
                MESS_NU_U0(g,:) = alpha_rw*log(exp(MESS_NU_U0_pre)./(exp(MESS_NU_U0_pre)+exp(MESS_NU_U1_pre))+1e-10)+(1-alpha_rw)*MESS_NU_U0_prev(g,:);
                MESS_NU_U1(g,:) = alpha_rw*log(exp(MESS_NU_U1_pre)./(exp(MESS_NU_U0_pre)+exp(MESS_NU_U1_pre))+1e-10)+(1-alpha_rw)*MESS_NU_U1_prev(g,:);
            end;
            PS_all_n = exp(MESS_NU_U0(g,:))./(exp(MESS_NU_U0(g,:))+exp(MESS_NU_U1(g,:))); %Probabilita' che i nodi siano bizantini sulla base dei messaggi
            rj = R(:,g);
            prs0 = calcola_pr_s_nu( rj, 0, PS_all_n.', eps, Pmal );
            prs1 = calcola_pr_s_nu( rj, 1, PS_all_n.', eps, Pmal );
            MESS_U0_pre = log(prs0);
            MESS_U1_pre = log(prs1); 
            
            MESS_U0(:,g) = alpha_rw*log(exp(MESS_U0_pre)./(exp(MESS_U0_pre)+exp(MESS_U1_pre))+1e-10)+(1-alpha_rw)*MESS_U0_prev(:,g);
            MESS_U1(:,g) = alpha_rw*log(exp(MESS_U1_pre)./(exp(MESS_U0_pre)+exp(MESS_U1_pre))+1e-10)+(1-alpha_rw)*MESS_U1_prev(:,g);
        end;
        % A questo punto possiamo far girare le forward and backward
        % propagation
        [ BM, ER_f, MESS_R0, MESS_R1 ] = MP_backward_log_def( MESS_U0, MESS_U1, s, eps, Pmal, rho);
        [ FM, ER_b, MESS_L0, MESS_L1 ] = MP_forward_log_def( MESS_U0, MESS_U1, s, eps, Pmal, rho, s1 );
        
        for g = 1:m
            MESS_D0_til_n = MESS_U0;
            MESS_D0_til_n(g,:) = 0; % Escludo il ramo verso cui devo mandare il messaggio
            MESS_D0_pre = sum(MESS_D0_til_n) + MESS_L0 + MESS_R0;
            MESS_D1_til_n = MESS_U1;
            MESS_D1_til_n(g,:) = 0; % Escludo il ramo verso cui devo mandare il messaggio
            MESS_D1_pre = sum(MESS_D1_til_n) + MESS_L1 + MESS_R1;
            MESS_D0(g,:) = alpha_rw*log(exp(MESS_D0_pre)./(exp(MESS_D0_pre)+exp(MESS_D1_pre))+1e-10)+(1-alpha_rw)*MESS_D0_prev(g,:);
            MESS_D1(g,:) = alpha_rw*log(exp(MESS_D1_pre)./(exp(MESS_D0_pre)+exp(MESS_D1_pre))+1e-10)+(1-alpha_rw)*MESS_D1_prev(g,:);
        end;
        % Infine si calcolano i messaggi di ritorno ai nodi sul loro stato 
        for g = 1:n
            rj = R(:,g);
            prs0 = calcola_pr_s_nu( rj, 0, ones(m,1), eps, Pmal ); % stato 0, nodo bizantino
            prs1 = calcola_pr_s_nu( rj, 1, ones(m,1), eps, Pmal ); % stato 1, nodo bizantino
            MESS_NU_D0_pre = log(prs0.*exp(MESS_D0(:,g)) + prs1.*exp(MESS_D1(:,g))+1e-10);
            prs0 = calcola_pr_s_nu( rj, 0, zeros(m,1), eps, Pmal ); % stato 0, nodo onesto
            prs1 = calcola_pr_s_nu( rj, 1, zeros(m,1), eps, Pmal ); % stato 1, nodo onesto
            MESS_NU_D1_pre = log(prs0.*exp(MESS_D0(:,g)) + prs1.*exp(MESS_D1(:,g))+1e-10);
            MESS_NU_D0(g,:) = (alpha_rw*log(exp(MESS_NU_D0_pre)./(exp(MESS_NU_D0_pre)+exp(MESS_NU_D1_pre))+1e-10)).'+(1-alpha_rw)*MESS_NU_D0_prev(g,:);
            MESS_NU_D1(g,:) = (alpha_rw*log(exp(MESS_NU_D1_pre)./(exp(MESS_NU_D0_pre)+exp(MESS_NU_D1_pre))+1e-10)).'+(1-alpha_rw)*MESS_NU_D1_prev(g,:);
        end;
        % Calcolo della probabilita' che i nodi siano bizantini
        Pbiz_pre = exp(sum(MESS_NU_D0))./(exp(sum(MESS_NU_D0))+exp(sum(MESS_NU_D1)));
        %PBIZ_ALL(IT).val = Pbiz_pre;
        MESS = FM+[MESS_R0;MESS_R1];
        PS = exp(MESS)./repmat(sum(exp(MESS)),2,1);
        %PS_ALL(tit,IT).val = PS;
        [maxv indxm] = max(MESS);
        EST = indxm-1;
        indx = find(abs(PS(1,:)-0.5) < 1e-8);
        EST(indx) = 0;
        ER(tit,IT) = sum(xor(EST,s))/n;
        
        MESS_L0_prev = MESS_L0;
        MESS_L1_prev = MESS_L1;
        MESS_R0_prev = MESS_R0;
        MESS_R1_prev = MESS_R1;
        MESS_U0_prev = MESS_U0;
        MESS_U1_prev = MESS_U1;
        MESS_D0_prev = MESS_D0;
        MESS_D1_prev = MESS_D1;
        % Inizializzazione dei messaggi relativi allo stato dei nodi
        MESS_NU_D0_prev = MESS_NU_D0;
        MESS_NU_D1_prev = MESS_NU_D1;
        MESS_NU_U0_prev = MESS_NU_U0;
        MESS_NU_U1_prev = MESS_NU_U1;

    end;
    
    %Calcolo con decisione a maggioranza
    MAJ_C = sum(R);
    EST_MAJ = zeros(1,n);
    EST_MAJ(MAJ_C > m/2) = 1;
    ER_MAJ(IT) = sum(xor(EST_MAJ,s))/n;
    
    %Calcolo con decisione a maggioranza con rimozione
    REP = zeros(m,1);
    for g = 1:m
        REP(g) = sum(xor(R(g,:),EST_MAJ));
    end;
    [sortv indxs] = sort(REP);
    indx_honest = indxs(1:m-Nbyz);
    %indx_honest = Nbyz+1:m;
    Rh = R(indx_honest,:);
    MAJ_C = sum(Rh);
    EST_MAJ = zeros(1,n);
    EST_MAJ(MAJ_C > length(indx_honest)/2) = 1;
    ER_MAJ_r(IT) = sum(xor(EST_MAJ,s))/n;
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


