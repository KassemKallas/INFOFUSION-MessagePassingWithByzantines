function[results] = CDC(sim_param)

Nprove = sim_param.Nprove;
PH1 = sim_param.PH1;
% PH0 = 1-PH1;
N = sim_param.N;
T = sim_param.T;
L=sim_param.L;
Pmal=sim_param.Pmal;
R = sim_param.R;
P = sim_param.P;
alfa = sim_param.alfa;
% epsilon = sim_param.epsilon;

% delta = sim_param.delta;

% delta_B = (1-delta)*(1-Pmal) +(delta)*Pmal;


%========Varsh and LLR parameters=================================
Pd_Hp = sim_param.Pd_Hp;
Pfa_Hp = sim_param.Pfa_Hp;
Pd_Bp = sim_param.Pd_Bp;
Pfa_Bp = sim_param.Pfa_Bp;

gammas = sim_param.gammas;

Pd_H = Pd_Hp;
Pfa_H = Pfa_Hp;
Pd_B = Pmal*(1-Pd_Bp)+(1-Pmal)*Pd_Bp;
Pfa_B = Pmal*(1-Pfa_Bp)+(1-Pmal)*Pfa_Bp;

Nsoglie_LLR = sim_param.Nsoglie_LLR;
Nerr_h = 0;
Nerr_b = 0;
N0=0;
N1=0;
Nerr_hr = zeros(length(gammas),1);
Nerr_br = zeros(length(gammas),1);
Nerr_hr_LLR = zeros(Nsoglie_LLR,1);
Nerr_br_LLR = zeros(Nsoglie_LLR,1);
Nerr_H = zeros(length(gammas),1);
Nerr_B = zeros(length(gammas),1);
Nerr_H_LLR = zeros(Nsoglie_LLR,1);
Nerr_B_LLR = zeros(Nsoglie_LLR,1);

%========End Varsh and LLR parameters=================================
% possible_states = sim_param.possible_system_states;

for np = 1:Nprove
%     if rem(np,10) == 0
%         fprintf('Simulation %d su %d\n',np,Nprove);
%     end;
%     
      %  ---------generate nodes randomly according to probability to be byzantine----
%     gen_prob = rand(N,1);
%     HO = find(gen_prob < 1-alfa);
%     BY = find(gen_prob >= 1-alfa);                                             % the probabilities are pseudorandom values drawn
                                                                               %from the standard uniform distribution on the open interval(0,1).
    
%     K = length(HO);
%     M = length(BY);
K = sim_param.K1;
M = N - K;
   % --------------FINISH GENERATING NODES---------------------------------
%   rd = rand(1,T);
%     P = zeros(1,T);
%     P(rd < PH1) = 1; %system state
    UH = zeros(K,T); % decision at honest
    UB = zeros(M,T); % decision at Byzantines
    D = zeros(1,T);
    LLRs_OUT = zeros(N,T);
    for t = 1:T
%         if P(t) == 1 % if the system state was 1
%             UH(:,t) = 1; % the reports of honest were first equal to one
%             GH = rand(K,1);
%             UH(GH < epsilon,t) = 0; % switch there report for honest
%             UB(:,t) = 1;
%             GB = rand(M,1);
%             UB(GB < delta_B,t) = 0;        
%         else % if the system state was zero
%             GH = rand(K,1);
%             UH(GH < epsilon,t) = 1;
%             GB = rand(M,1);
%             UB(GB < delta_B,t) = 1;
%         end;
        col_tmp = R(:,t);
        UB(:,t) = col_tmp(1:M);
        UH(:,t) = col_tmp(M+1:N);

   %================================ Varshney and LLR initialization =============    
        Prob_err = Pmal*M/N;
        
        U_ALL = [UH(:,t);UB(:,t)];
        Num_ones = length(find(U_ALL == 1));
        Num_zeros = length(find(U_ALL == 0));
        
        
        
       % alpha = M/N;
       alpha = alfa;
        P1 = (1-alpha)*Pfa_H+alpha*Pfa_B;
        P2 = (1-alpha)*Pd_H+alpha*Pd_B;
        PUH0 = ((1-Prob_err)*P1+Prob_err*(1-P1))^Num_ones * ((1-Prob_err)*(1-P1)+Prob_err*P1)^Num_zeros;
        PUH1 = ((1-Prob_err)*P2+Prob_err*(1-P2))^Num_ones * ((1-Prob_err)*(1-P2)+Prob_err*P2)^Num_zeros;
        for dec = 1:N
            if U_ALL(dec) == 0
                PUH0d = PUH0/(((1-Prob_err)*(1-P1)+Prob_err*P1));
                PUH1d = PUH1/((1-Prob_err)*(1-P2)+Prob_err*P2);
                Px0U = (1-Prob_err)*((1-P1)*(1-PH1)*PUH0d+(1-P2)*PH1*PUH1d);
                Px1U = Prob_err*(P1*(1-PH1)*PUH0d+P2*PH1*PUH1d);
            else
                PUH0d = PUH0/((1-Prob_err)*P1+Prob_err*(1-P1));
                PUH1d = PUH1/((1-Prob_err)*P2+Prob_err*(1-P2));
                Px0U = Prob_err*((1-P1)*(1-PH1)*PUH0d+(1-P2)*PH1*PUH1d);
                Px1U = (1-Prob_err)*(P1*(1-PH1)*PUH0d+P2*PH1*PUH1d);
            end;
            LLRs_OUT(dec,t) = abs(log(Px0U/Px1U));
        end;
        
%================================ Varshney and LLR initialization =============    


        if sum([UH(:,t);UB(:,t)]) >= L % obtain the majority rule result
            D(t) = 1;
        else
            D(t) = 0;
        end;
        
        R_matrix(:,t) = U_ALL;
        
    end
    
%============================Decoding Varshney and LLR ==========================================
REL = sum(LLRs_OUT,2);
    if  (max(REL)- min(REL)) > 0 
        
        for i=1:Nsoglie_LLR
         SOGLIA_LLRs(i) = min(REL) + ((max(REL)-min(REL))/Nsoglie_LLR)*i;
        end
    else
        SOGLIA_LLRs = mean(REL)*ones(1,Nsoglie_LLR);
    end;
    
    for is = 1:Nsoglie_LLR
        SOGLIA_LLR = SOGLIA_LLRs(is);
        Nerr_H_LLR(is) = Nerr_H_LLR(is) + length(find(REL(1:K) < SOGLIA_LLR));
        Nerr_B_LLR(is) = Nerr_B_LLR(is) + length(find(REL(K+1:N) >= SOGLIA_LLR));
    end;
    
    
    Dall_H = repmat(D,K,1);
    Errs_H = xor(Dall_H,UH);
    Dall_B = repmat(D,M,1);
    Errs_B = xor(Dall_B,UB);
    eta_H = sum(Errs_H,2);
    eta_B = sum(Errs_B,2);
    for ig = 1:length(gammas)
        gamma = gammas(ig);
        Nerr_H(ig) = Nerr_H(ig) + length(find(eta_H > gamma));
        Nerr_B(ig) = Nerr_B(ig) + length(find(eta_B <= gamma));
    end;
    
   
    
    indx = find(P == 1);
    Nerr_h = Nerr_h + length(find(D(indx) == 0));
    indx = find(P == 0);
    Nerr_b = Nerr_b + length(find(D(indx) == 1));
    % variabili che servono nel caso non simmetrico
    N0 = N0 + numel(find(P==0));
    N1 = N1 + numel(find(P==1));
   
    %Valutazione delle prestazioni dopo rimozione
    for ig = 1:length(gammas)
        Dr = 0*D;
        gamma = gammas(ig);
        indxH = find(eta_H <= gamma);
        indxB = find(eta_B <= gamma);
        if length(indxH)+length(indxB) == 0
            indxH = 1:length(eta_H);
            indxB = 1:length(eta_B);
        end;        
        for t = 1:T
            %Decisione
            if sum([UH(indxH,t);UB(indxB,t)]) >= (length(indxH)+length(indxB))/2
                Dr(t) = 1;
            else
                Dr(t) = 0;
            end;
        end
        indx = find(P == 1);
        Nerr_hr(ig) = Nerr_hr(ig) + length(find(Dr(indx) == 0));
        indx = find(P == 0);
        Nerr_br(ig) = Nerr_br(ig) + length(find(Dr(indx) == 1));
    end;
    %Valutazione delle prestazioni dopo rimozione nel caso LLR
    for is = 1:Nsoglie_LLR
        Dr = 0*D;
        SOGLIA_LLR = SOGLIA_LLRs(is);
        indxH = find(REL(1:K) >= SOGLIA_LLR);
        indxB = find(REL(K+1:N) >= SOGLIA_LLR);
        if length(indxH)+length(indxB) == 0
            indxH = 1:length(eta_H);
            indxB = 1:length(eta_B);
        end;
        for t = 1:T
            %%%%Decisione
            if sum([UH(indxH,t);UB(indxB,t)]) >= (length(indxH)+length(indxB))/2
                Dr(t) = 1;
            else
                Dr(t) = 0;
            end;            
        end
        indx = find(P == 1);
        Nerr_hr_LLR(is) = Nerr_hr_LLR(is) + length(find(Dr(indx) == 0));
        indx = find(P == 0);
        Nerr_br_LLR(is) = Nerr_br_LLR(is) + length(find(Dr(indx) == 1));
    end;
    
%================== END DECODING VARSHNEY AND LLR ==============================





%===================DECODING Independent states====================

% for i=1:2^T % means over all possibile states
%     select_state = possible_states(i,:); % take a vector from truth table
%     
%     nonzeros = bsxfun(@minus,R_matrix,select_state);
%     
%     neq_state = T - sum(nonzeros~=0,2);
%     
%     Tminusneq = bsxfun(@minus,T,neq_state); 
%     
%     for pmal_dec_guess = 1:1:6
%     Pmal_guess = (pmal_dec_guess+4)/10;
%     delta_B_guess = (1-delta)*(1-Pmal_guess) +(delta)*Pmal_guess;    
%     
%     f1 = (1-alfa)*((1-epsilon).^neq_state).*((epsilon.^Tminusneq));
%     f2 = (alfa)*((1-delta_B_guess).^neq_state).*((delta_B_guess.^Tminusneq));
%     f3 = f1 + f2;
%     prod_PrRsT(i,pmal_dec_guess) = log(prod(f3));
% 
%     end
% end
% %===================End DECODING Independent states====================
% 
% 
% %===================Decision Independent states====================
% [val idx]=max(prod_PrRsT,[],1);
% 
% for iii=1:length(val)
%     
% decision_fusion(iii,:) = possible_states(idx(iii),:);
% err_nb_eq4 (np,iii) = numel(find(P~=decision_fusion(iii,:)));
% 
% end
% results.decision_fusion= decision_fusion;
%===================End decision Independent states====================

diff_majority = P-D;

err_nb_majority (np) = nnz(diff_majority); % number of error decision for majority rule

total_nb_trials(np) = length(P);

end;

% error_eq4 = sum(err_nb_eq4)/sum(total_nb_trials);
% 
% 
% 
% results.error_eq4 = error_eq4;

results.Nerr_hr = Nerr_hr;
results.Nerr_br = Nerr_br;
results.Nerr_hr_LLR = Nerr_hr_LLR;
results.Nerr_br_LLR = Nerr_br_LLR;
results.Nerr_h = Nerr_h;
results.Nerr_b = Nerr_b;
%================RESULTS VARSHNEY and LLR=================================
% results.PD = 1-Nerr_h/Nprove/T;
% results.PFA = Nerr_b/Nprove/T;
% for ig = 1:length(gammas)
%     results.PDr(ig) = 1-Nerr_hr(ig)/Nprove/T;
%     results.PFAr(ig) = Nerr_br(ig)/Nprove/T;
%     results.PD_IDB(ig) = 1-Nerr_H(ig)/K/Nprove;
%     results.PFA_IDB(ig) = Nerr_B(ig)/M/Nprove;
%     results.P_ISO_H(ig) = Nerr_H(ig)/K/Nprove;
%     results.P_ISO_B(ig) = 1 - Nerr_B(ig)/M/Nprove;
% end;
% for is = 1:Nsoglie_LLR
%     results.PD_IDB_LLR(is) = 1-Nerr_H_LLR(is)/K/Nprove; % 1 - P_ISO^H 
%     results.PFA_IDB_LLR(is) = Nerr_B_LLR(is)/M/Nprove;  % 1 - P_ISO^B = P_NONISO^B 
%     results.P_ISO_H_LLR(is) = Nerr_H_LLR(is)/K/Nprove; % P_ISO^H 
%     results.P_ISO_B_LLR(is) = 1 - Nerr_B_LLR(is)/M/Nprove;  % P_ISO^B 
%     results.PDr_LLR(is) = 1-Nerr_hr_LLR(is)/Nprove/T;
%     results.PFAr_LLR(is) = Nerr_br_LLR(is)/Nprove/T;
% end;


%================END RESULTS VARSHNEY and LLR================================
% error_majority = sum(err_nb_majority)/sum(total_nb_trials);
% results.error_majority = error_majority;
end

















