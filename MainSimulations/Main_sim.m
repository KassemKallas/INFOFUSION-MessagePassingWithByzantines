Niteraz = 5000;

Pmal = 1;
rho = 0.1;
eps = 0.2;
alpha = 0.45;
s1 = 0;
n = 30;
m = 20;

for IT = 1:Niteraz
    if rem(IT,50) == 0
        fprintf('Iterazione %d su %d\n',IT,Niteraz);
    end;
    [ R, s ] = crea_reports_Markov( n, m, alpha, Pmal, eps, rho, s1);
    %rho = 0.5;
    
    %[ BM, ER_f ] = MP_backward( R, s, eps, Pmal, alpha, rho);
    %[ FM, ER_b ] = MP_forward( R, s, eps, Pmal, alpha, rho, s1 );
    
    [ BM, ER_f ] = MP_backward_log( R, s, eps, Pmal, alpha, rho);
    [ FM, ER_b ] = MP_forward_log( R, s, eps, Pmal, alpha, rho, s1 );
    
    %MESS = BM.*FM;
    %MESS = MESS./repmat(sum(MESS),2,1);
    MESS = BM+FM;
    [maxv indxm] = max(MESS);
    EST = indxm-1;
    ER(IT) = sum(xor(EST,s))/n;
    
    %Calcolo con decisione a maggioranza
    MAJ_C = sum(R);
    EST_MAJ = zeros(1,n);
    EST_MAJ(MAJ_C > m/2) = 1;
    ER_MAJ(IT) = sum(xor(EST_MAJ,s))/n;
    
    PS = exp(MESS)./repmat(sum(exp(MESS)),2,1);
    %Individuazione bizantini
    % Calcolo di PU, probabilita' delle osservazioni
    PU = 0*PS;
    PU(1,:) = PS(1,:)*(1-eps)+PS(2,:)*eps;
    PU(2,:) = PS(2,:)*(1-eps)+PS(1,:)*eps;
    %Creazione delle cattive reputazioni per ciascun nodo
    nu_all = zeros(m,n);
    for j = 1:m
        nu_all(j,:) = PU(R(j,:)+1+[0:n-1]*2)*(1-Pmal)+PU(~R(j,:)+1+[0:n-1]*2)*Pmal;
    end;
    nu = mean(nu_all,2);%sum(log(nu_all+1e-10),2);
    %Creazione delle buone reputazioni per ciascun nodo
    mu_all = zeros(m,n);
    for j = 1:m
        %Probabilita' che l'osservazione u sia uguale al report
        mu_all(j,:) = PU(R(j,:)+1+[0:n-1]*2);
    end;

    mu = mean(mu_all,2);%sum(log(nu_all+1e-10),2);
    %
    Nbyz = round(alpha*m);
    [sortv indxs] = sort(-mu);
    indx_honest = indxs(1:m-Nbyz);
    
    Rh = R(indx_honest,:);
    Pmal_r = 0;%Pmal;
    %[ BM, ER_f ] = MP_backward( Rh, s, eps, Pmal, alpha, rho);
    %[ FM, ER_b ] = MP_forward( Rh, s, eps, Pmal, alpha, rho, s1 );
    [ BM, ER_f ] = MP_backward_log( Rh, s, eps, Pmal_r, alpha, rho);
    [ FM, ER_b ] = MP_forward_log( Rh, s, eps, Pmal_r, alpha, rho, s1 );
    %MESS = BM.*FM;
    %MESS = MESS./repmat(sum(MESS),2,1);
    MESS = BM+FM;
    
    [maxv indxm] = max(MESS);
    EST = indxm-1;
    ER_r(IT) = sum(xor(EST,s))/n;
    
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
fprintf('Errr rate con decisione a maggioranza dopo la rimozione ideale = %f\n',mean(ER_MAJ_r));
fprintf('Errr rate con algoritmo BCJR prima della rimozione = %f\n',mean(ER));
fprintf('Errr rate con algoritmo BCJR dopo la rimozione = %f\n',mean(ER_r));


