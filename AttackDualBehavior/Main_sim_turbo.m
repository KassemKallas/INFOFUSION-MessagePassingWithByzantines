Niteraz = 3000;

Pmal = 0.5;
rho = 0.5;
eps = 0.1;
alpha = 0.4;
s1 = -1;%-1 = unknown
n = 20;
m = 10;
Num_iteraz_turbo = 5;
ER = zeros(Num_iteraz_turbo,Niteraz);
Nbyz = round(alpha*m);
for IT = 1:Niteraz
    if rem(IT,50) == 0
        fprintf('Iterazione %d su %d\n',IT,Niteraz);
    end;
    %[ R, s ] = crea_reports_Markov( n, m, alpha, Pmal, eps, rho, s1);
    [ R, s, Nbyz ] = crea_reports_Markov_stat( n, m, alpha, Pmal, eps, rho, s1);
    % Inizializzazione dei messaggi di stato dei nodi
    q = zeros(2,m);
    q(1,:) = alpha;
    q(2,:) = 1-alpha;
    % Iterazioni turbo
    
    for tit = 1:Num_iteraz_turbo
        [ BM, ER_f ] = MP_backward_log_new( R, s, eps, Pmal, q, rho);
        [ FM, ER_b ] = MP_forward_log_new( R, s, eps, Pmal, q, rho, s1 );
        
        %MESS = BM.*FM;
        %MESS = MESS./repmat(sum(MESS),2,1);
        MESS = BM+FM;
        PS = exp(MESS)./repmat(sum(exp(MESS)),2,1);
        %Adesso occorre calcolare q per il prossimo giro
        Q = zeros(2,m);
        for g = 1:m
            for nu = 0:1 %0 = Byzantino, 1 = onesto
                if nu == 1
                    iota = eps;
                else
                    iota = eps*(1-Pmal)+(1-eps)*Pmal;
                end;
                M = zeros(2,n);
                for stato = 0:1 %0 = idle, 1 = busy
                    Prob = PS(stato+1,:);
                    deltas = zeros(n,1);
                    deltas(R(g,:) == stato) = 1;
                    M(stato+1,:) = ((1-iota)*deltas+iota*(~deltas)).'.*Prob;
                end;
                Q(nu+1,g) = sum(log(sum(M)));
            end;
        end;
        q = exp(Q)./repmat(sum(exp(Q)),2,1);
        [maxv indxm] = max(MESS);
        EST = indxm-1;
        ER(tit,IT) = sum(xor(EST,s))/n;
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
    
end;
fprintf('Errr rate con decisione a maggioranza = %f\n',mean(ER_MAJ));
fprintf('Errr rate con decisione a maggioranza dopo la rimozione = %f\n',mean(ER_MAJ_r));
for l = 1:Num_iteraz_turbo
    fprintf('Errr rate con Turbo decoding dopo %d iterazioni = %f\n',l,mean(ER(l,:)));
end;

%fprintf('Errr rate con algoritmo BCJR prima della rimozione = %f\n',mean(ER));
%fprintf('Errr rate con algoritmo BCJR dopo la rimozione = %f\n',mean(ER_r));


