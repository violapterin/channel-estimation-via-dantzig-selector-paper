%% Simulation: Proposed Hybrid Beamforming System
% Reference: "Hybrid Beamforming Design for Multiuser Transmission in Massive MIMO Millimeter Wave System"
% Author:Chen-Chieh Hong
% Date: 2017/07/21
% System: Multi-user scenario (MU-MISO)
% Idea:
% 1. Number of RF chains less than users: Eigen Beam
% 2. Number of RF chains equal to users: GSINR-FB in analog part 
% 3. number of RF chains larger than users: Decomposition Design,
%                                           Orthogonal Space Design
% Brief description: 
% 1. Scope (from outer to inner): SNR -> channel realization -> number of bits in phase shifter -> number of RF chains 
% 2. Consider homogenerous case, i.e. the distance between each user to
%    BS is equal. Thus the definition of SNR is received SNR (pathloss)
% 3. The region of number of RF chains: Num_UE-1 ~ 2*Num_UE
% 4. Result: averge sum rate
%% =================== Initialization ==================== %%
tic;clc;clear;close all;
% ------------------------System Parameters---------------------------------
Num_BS_Antennas = 64;                        % number of BS antennas
BSAntennas_Index = 0:1:Num_BS_Antennas-1;    % Indices of the BS Antennas
Num_UE = 8;                                 % number of users
Num_BS_RFchains_t = Num_UE-1:1:2*Num_UE;     % number of BS RF chains
%Num_BS_RFchains_t = Num_UE-1:1:2*Num_UE;     % number of BS RF chains
AvgRate_com = zeros(length(Num_BS_RFchains_t),11);     % average rate when Q=inf
AvgRate_com_q_1 = zeros(length(Num_BS_RFchains_t),11); % average rate when Q=1
AvgRate_com_q_2 = zeros(length(Num_BS_RFchains_t),11); % average rate when Q=2
AvgRate_com_q_3 = zeros(length(Num_BS_RFchains_t),11); % average rate when Q=3
P = Num_UE;                                  % Total power constrain (equal to number of users)
Num_UE_Antennas = 1;                         % number of UE antennas
UEAntennas_Index = 0:1:Num_UE_Antennas-1;    % Indices of the UE Antennas
Num_Qbits = [0 1 2];                         % Number of phase shifters quantization bits, 0: infinite resolution
% ---------------------Channel Parameters ---------------------------------
Num_paths = 15; % Number of channel paths
%% =================== Simulation ============================================= %%
SNR_dBa=-10:2:10;
% nITER= 1000;
nITER= 10;
% ------------------------Different SNR-------------------------  
for SNR_ind=1:length(SNR_dBa)
    SNR_ind
    SNR=10^(.1*SNR_dBa(SNR_ind));
    No=P/SNR;
    % ----------------------Iteration of Independent Realizations-------------------------
    for ITER= 1:nITER
        ITER;
        % --------------------------Channel Generation---------------------------------                
        % Channel for different user
        for K = 1:Num_UE          
            AoD=pi*rand(1,Num_paths)-pi/2;
            AoD_t(ITER,K,:)=AoD;
            AoA=pi*rand(1,Num_paths)-pi/2;
            AoA_t(ITER,K,:)=AoA;
            alpha=(sqrt(1/2)*sqrt(1/Num_paths)*(randn(1,Num_paths)+1j*randn(1,Num_paths)));
            Alpha_t_nor(ITER,K,:)=abs(alpha);
            Channel=zeros(Num_UE_Antennas,Num_BS_Antennas);
            % different path
            for path=1:1:Num_paths
                Abh(:,path)=sqrt(1/Num_BS_Antennas)*exp(1j*pi*BSAntennas_Index*sin(AoD(path)));
                Amh(:,path)=sqrt(1/Num_UE_Antennas)*exp(1j*pi*UEAntennas_Index*sin(AoA(path)));
                Channel=Channel+sqrt(Num_BS_Antennas*Num_UE_Antennas)*alpha(path)*Amh(:,path)*Abh(:,path)';
            end
            H(K,:)= Channel;
        end        
        % --------------------------Different Quantized Version---------------------------------
        for num_quant_ind = 1:length(Num_Qbits)           
            if Num_Qbits(num_quant_ind) == 0
                Q_set = -1;
            else               
                Q_set = 2*pi*(0:2^Num_Qbits(num_quant_ind)-1)/2^Num_Qbits(num_quant_ind);
            end
            %---------------------- Different RF chains-------------------------------               
            for rf_ind=1:length(Num_BS_RFchains_t)
                Num_BS_RFchains = Num_BS_RFchains_t(rf_ind);
                %---------------------- Analog: Eigen beam (RF chain < UE)-------------------------------               
                if Num_BS_RFchains < Num_UE                              
                    H_temp_UL=H';                       
                    S_Cov = zeros(Num_BS_Antennas,Num_BS_Antennas);
                    for k=1:Num_UE
                        S_Cov =  S_Cov + H_temp_UL(:,k)*H_temp_UL(:,k)';                                                               
                    end
                    N_Cov =  No*eye(Num_BS_Antennas);                                
                    [V_eig,D_eig] = eig(S_Cov,N_Cov);
                    [y_eig,ind_eig] = sort(max(D_eig),'descend');
                    V_RF_less_f = zeros(Num_BS_Antennas,Num_BS_RFchains);
                    V_RF_less = zeros(Num_BS_Antennas,Num_BS_RFchains);
                    for k=1:Num_BS_RFchains
                        V_RF_less_f(:,k) = V_eig(:,ind_eig(k));
                        if Q_set == -1
                            V_RF_less(:,k)=exp(1j*angle(V_RF_less_f(:,k)));
                        else
                            V_RF_less(:,k)=exp(1j*angle(V_RF_less_f(:,k)));
                            for F_ind = 1:Num_BS_Antennas
                                [y,Q_ind] = min(abs(V_RF_less(F_ind,k)-exp(1j*Q_set)));
                                V_RF_less(F_ind,k) = exp(1j*Q_set(Q_ind));
                            end
                        end               
                    end
                    H_eq_less=H*V_RF_less;
                %---------------------- Analog: GMSINR-FB (RF chain = UE)-------------------------------                   
                else
                    H_temp_UL=H';            
                    % RF chain = UE
                    V_RF_f = zeros(Num_BS_Antennas,Num_UE);
                    V_RF = zeros(Num_BS_Antennas,Num_UE);
                    for k=1:Num_UE
                        S_Cov = zeros(Num_BS_Antennas,Num_BS_Antennas);
                        S_Cov = H_temp_UL(:,k)*H_temp_UL(:,k)';
                        [x,inter_ind]=find(k~=[1:Num_UE]);
                        N_Cov = zeros(Num_BS_Antennas,Num_BS_Antennas);
                        for u=1:length(inter_ind)
                            N_Cov =  N_Cov + H_temp_UL(:,inter_ind(u))*H_temp_UL(:,inter_ind(u))';
                        end
                        N_Cov =  N_Cov + No*eye(Num_BS_Antennas); %+ (((pi)^2)/6)*eye(Num_BS_Antennas);                      
                        [V_f,D_f] = eig(S_Cov,N_Cov);
                        [y,ind] = max(max(D_f));
                        V_RF_f(:,k) = V_f(:,ind);              
                        V_RF_f(:,k)=sqrt(1/(V_RF_f(:,k)'*V_RF_f(:,k)))*V_RF_f(:,k);
                        % Infinite
                        if  Q_set == -1
                            V_RF(:,k)=sqrt(1/Num_BS_Antennas)*exp(1j*angle(V_f(:,ind)));
                        % Quantization
                        else
                            V_RF(:,k)=exp(1j*angle(V_f(:,ind)));
                            for F_ind = 1:Num_BS_Antennas
                                [y,Q_ind] = min(abs(V_RF(F_ind,k)-exp(1j*Q_set)));
                                V_RF(F_ind,k) = sqrt(1/Num_BS_Antennas)*exp(1j*Q_set(Q_ind));
                            end
                        end               
                    end
                    H_eq=H*V_RF;
                    %---------------------- Analog: Decom / Ortho (UE < RF chain <= 2UE)-------------------------------                   
                    if Num_BS_RFchains > Num_UE
                        % Decomposition Design
                        % Calculate SINR_UL for RF chain = UE
                        for k=1:Num_UE
                            desired_PW_UL(k) = abs((V_RF(:,k)'*H_temp_UL(:,k)))^2;
                            [x,Non_zero_ind]=find(k~=[1:Num_UE]);               
                            N_Cov = zeros(Num_BS_Antennas,Num_BS_Antennas);
                            for j=1:length(Non_zero_ind)
                                N_Cov =  N_Cov + H_temp_UL(:,Non_zero_ind(j))*H_temp_UL(:,Non_zero_ind(j))';             
                            end
                            N_Cov =  N_Cov + sqrt(1/Num_BS_Antennas)*No*eye(Num_BS_Antennas);
                            interference_PW_UL(k) = real(V_RF(:,k)'*N_Cov*V_RF(:,k));
                            SINR_UL(k) = desired_PW_UL(k)/interference_PW_UL(k);
                        end
                        [y_order,ind_order] = sort(SINR_UL,'descend');                       
                        % more RF cahins allocate to the best user, combination
                        V_RF_decom = zeros(Num_BS_Antennas,Num_BS_RFchains);
                        V_RF_decom(:,1:Num_UE) = V_RF;            
                        for j=1:Num_BS_RFchains-Num_UE
                            nor_f = 2/max(abs(V_RF_f(:,ind_order(j))));                
                            V_RF_f_nor = nor_f*V_RF_f(:,ind_order(j));              
                            % Structure 1               
                            Theta = angle(V_RF_f_nor) + acos(abs(V_RF_f_nor)/2);
                            Phi = angle(V_RF_f_nor) - acos(abs(V_RF_f_nor)/2);
                            if  Q_set == -1                    
                                % Structure 1
                                V_RF_decom(:,ind_order(j)) = sqrt(1/Num_BS_Antennas)*exp(1j*Theta);
                                V_RF_decom(:,Num_UE+j) = sqrt(1/Num_BS_Antennas)*exp(1j*Phi);                    
                            else                                      
                                % Structure 1
                                V_RF_decom(:,ind_order(j)) = exp(1j*Theta) ;                    
                                V_RF_decom(:,Num_UE+j) = exp(1j*Phi);
                                for F_ind = 1:Num_BS_Antennas
                                    [y_t,Q_ind_t] = min(abs(V_RF_decom(F_ind,ind_order(j))-exp(1j*Q_set)));
                                    V_RF_decom(F_ind,ind_order(j)) = sqrt(1/Num_BS_Antennas)*exp(1j*Q_set(Q_ind_t));
                                    [y_p,Q_ind_p] = min(abs(V_RF_decom(F_ind,Num_UE+j)-exp(1j*Q_set)));
                                    V_RF_decom(F_ind,Num_UE+j) = sqrt(1/Num_BS_Antennas)*exp(1j*Q_set(Q_ind_p));                        
                                end
                            end
                        end                       
                        V_RF_decom_set = zeros(Num_BS_Antennas,2);
                        Coeff = zeros(Num_BS_RFchains-Num_UE,2);
                        for tt=1:Num_BS_RFchains-Num_UE
                            %nor_f = 2/max(abs(V_RF_5_f(:,ind_order(tt))));
                            nor_f = (2/sqrt(Num_BS_Antennas))/max(abs(V_RF_f(:,ind_order(j))));
                            V_RF_f_nor = nor_f*V_RF_f(:,ind_order(tt));
                            V_RF_decom_set = [V_RF_decom(:,ind_order(tt)) V_RF_decom(:,Num_UE+tt)];
                            % LS
                            Coeff(tt,:)=(V_RF_decom_set'*V_RF_decom_set)\(V_RF_decom_set'*V_RF_f_nor); % non-singular in quantized version ?
%                             % pseudo inverse
%                             [U_com,S_com,V_com] = svd(V_RF_5_com_set);
%                             Sub_S = S_com(1:2,1:2);
%                             S_com_i = zeros(2,Num_BS_Antennas);
%                             S_com_i(1:2,1:2) = inv(Sub_S);
%                             Coeff(tt,:) = V_com*S_com_i*U_com'*V_RF_5_f_nor;
                        end           
                        Combination = eye(Num_BS_RFchains,Num_UE);
                        for tt=1:Num_BS_RFchains-Num_UE
                            Combination(ind_order(tt),ind_order(tt)) = Coeff(tt,1);
                            Combination(Num_UE+tt,ind_order(tt)) = Coeff(tt,2);
                        end
                        V_RF_decom_total=V_RF_decom*Combination;
                        H_eq_decom=H*V_RF_decom_total;
                        
                        % Orthogonal Space Design                                     
                        V_RF_ortho_com = eye(Num_BS_Antennas,Num_UE);
                        Cov_V_RF_ortho_com = zeros(Num_BS_Antennas,Num_BS_Antennas);
                        for j=1:Num_UE
                            % orthogonal projection
                            q_1 = V_RF(:,j);                
                            V_RF_ortho_com(:,j) = V_RF_f(:,j) - ((q_1'*V_RF_f(:,j))/(q_1'*q_1))*q_1;
                            Cov_V_RF_ortho_com = Cov_V_RF_ortho_com + V_RF_ortho_com(:,j)*V_RF_ortho_com(:,j)';
                        end            
                        [V_com,D_com] = eig(Cov_V_RF_ortho_com);
                        %[V_com,D_com] = eigs(Cov_V_RF_5_comple,Num_UE,'lm');
                        %[V_com,D_com] = eigs(Cov_V_RF_5_comple,Num_BS_RFchains-Num_UE,'lm');
                        %[V_com,D_com] = eigs(Cov_V_RF_5_comple);
                        V_RF_ortho = zeros(Num_BS_Antennas,Num_BS_RFchains);
                        V_RF_ortho(:,1:Num_UE) = V_RF;
                        [y_com,ind_com] = sort(max(D_com),'descend');
                        for j=1:Num_BS_RFchains-Num_UE
                            if Q_set == -1
                                V_RF_ortho(:,Num_UE+j) = sqrt(1/Num_BS_Antennas)*exp(1j*angle(V_com(:,ind_com(j))));
                            else
                                V_RF_ortho(:,Num_UE+j) = exp(1j*angle(V_com(:,ind_com(j))));
                                for F_ind = 1:Num_BS_Antennas
                                    [y,Q_ind] = min(abs(V_RF_ortho(F_ind,Num_UE+j)-exp(1j*Q_set)));
                                    V_RF_ortho(F_ind,Num_UE+j) = sqrt(1/Num_BS_Antennas)*exp(1j*Q_set(Q_ind));
                                end
                            end
                        end
                        H_eq_ortho=H*V_RF_ortho;                      
                    end
                end
            end
            %----------------------Digital: ZF / MMSE-------------------------------                    
            % Hybrid, less, digital MMSE
            V_D_less_hat = H_eq_less'/(H_eq_less*H_eq_less'+No*eye(Num_UE));                                      
            Q_less = (V_D_less_hat)'*(V_RF_less)'*(V_RF_less)*V_D_less_hat;                                      
            % Hybrid, equal, digital ZF
            V_D_hat = H_eq'/(H_eq*H_eq');                                                      
            Q = (V_D_hat)'*(V_RF)'*(V_RF)*V_D_hat;                   
            % Hybrid, large_decom, digital ZF
            V_D_decom_hat = H_eq_decom'/(H_eq_decom*H_eq_decom');                                                      
            Q_decom = (V_D_decom_hat)'*(V_RF_decom_total)'*(V_RF_decom_total)*V_D_decom_hat;
            % Hybrid, large_ortho, digital ZF
            V_D_ortho_hat = H_eq_ortho'/(H_eq_ortho*H_eq_ortho');                                                      
            Q_ortho = (V_D_ortho_hat)'*(V_RF_ortho)'*(V_RF_ortho)*V_D_ortho_hat;
            % Fully, digital ZF
            V_D_fully_zf_hat = H'/(H*H'); %H'*inv(H*H');
            Q_fully_zf = (V_D_fully_zf_hat)'*V_D_fully_zf_hat;                     

            % equal power allocation         
            for k=1:Num_UE
                p_a_e(k) = P/Num_UE;
            end                       
            for k=1:Num_UE                
                Pr_less(k,k) = real((1/Q_less(k,k))*p_a_e(k));               
                Pr(k,k) = real((1/Q(k,k))*p_a_e(k));
                Pr_decom(k,k) = real((1/Q_decom(k,k))*p_a_e(k));
                Pr_ortho(k,k) = real((1/Q_ortho(k,k))*p_a_e(k));
                Pr_fully_zf(k,k) = real((1/Q_fully_zf(k,k))*p_a_e(k));               
            end
            V_D_less = V_D_less_hat*sqrt(Pr_less);
            V_D = V_D_hat*sqrt(Pr);        
            V_D_decom = V_D_decom_hat*sqrt(Pr_decom);
            V_D_ortho = V_D_ortho_hat*sqrt(Pr_ortho);
            V_D_fully_zf = V_D_fully_zf_hat*sqrt(Pr_fully_zf);                            
            % ------------------------Sum Rate Calculation (per SNR per channel)---------------------------------                      
            % Hybrid, RF chain < Num_UE, A: Eigen-beam, D: MMSE
            for ind_D = 1:Num_UE
                desired_PW_p = abs((H(ind_D,:)*V_RF_less*V_D_less(:,ind_D)))^2;
                [x,Non_zero_ind]=find(ind_D~=[1:Num_UE]);
                interference_PW_p = sum(abs((H(ind_D,:)*V_RF_less*V_D_less(:,Non_zero_ind))).^2);
                SINR_less(num_quant_ind,ind_D) = desired_PW_p/(interference_PW_p+No);
            end
            R_less(num_quant_ind,ITER) = sum(log2(1+SINR_less(num_quant_ind,:)));            
            % Hybrid, RF chain = Num_UE, A: GMSINR-FB, D: ZF
            for ind_D = 1:Num_UE
                desired_PW_p = abs((H(ind_D,:)*V_RF*V_D(:,ind_D)))^2;
                [x,Non_zero_ind]=find(ind_D~=[1:Num_UE]);
                interference_PW_p = sum(abs((H(ind_D,:)*V_RF*V_D(:,Non_zero_ind))).^2);
                SINR(num_quant_ind,ind_D) = desired_PW_p/(interference_PW_p+No);
            end
            R(num_quant_ind,ITER) = sum(log2(1+SINR(num_quant_ind,:)));                        
            % Hybrid, RF chain > Num_UE, A: GMSINR-FB + Decom, D: ZF
            for ind_D = 1:Num_UE
                desired_PW_p = abs((H(ind_D,:)*V_RF_decom_total*V_D_decom(:,ind_D)))^2;
                [x,Non_zero_ind]=find(ind_D~=[1:Num_UE]);
                interference_PW_p = sum(abs((H(ind_D,:)*V_RF_decom_total*V_D_decom(:,Non_zero_ind))).^2);
                SINR_decom(num_quant_ind,ind_D) = desired_PW_p/(interference_PW_p+No);
            end
            R_decom(num_quant_ind,ITER) = sum(log2(1+SINR_decom(num_quant_ind,:)));
            % Hybrid, RF chain > Num_UE, A: GMSINR-FB + Ortho, D: ZF
            for ind_D = 1:Num_UE
                desired_PW_p = abs((H(ind_D,:)*V_RF_ortho*V_D_ortho(:,ind_D)))^2;
                [x,Non_zero_ind]=find(ind_D~=[1:Num_UE]);
                interference_PW_p = sum(abs((H(ind_D,:)*V_RF_ortho*V_D_ortho(:,Non_zero_ind))).^2);
                SINR_ortho(num_quant_ind,ind_D) = desired_PW_p/(interference_PW_p+No);
            end
            R_ortho(num_quant_ind,ITER) = sum(log2(1+SINR_ortho(num_quant_ind,:)));
            % Fully digital, ZF
            for ind_D = 1:Num_UE
                desired_PW_p_fully_zf = abs((H(ind_D,:)*V_D_fully_zf(:,ind_D)))^2;
                [x,Non_zero_ind]=find(ind_D~=[1:Num_UE]);
                interference_PW_p_fully_zf = sum(abs((H(ind_D,:)*V_D_fully_zf(:,Non_zero_ind))).^2);
                SINR_fully_zf(num_quant_ind,ind_D) = desired_PW_p_fully_zf/(interference_PW_p_fully_zf+No);
            end
            R_fully_zf(num_quant_ind,ITER) = sum(log2(1+SINR_fully_zf(num_quant_ind,:)));                                 
        end
    end
    % ------------------------Sum Rate Calculation (per SNR)---------------------------------
    for num_quant_ind = 1:length(Num_Qbits)       
        AvgRate_less(num_quant_ind,SNR_ind) = sum(R_less(num_quant_ind,:))/nITER;
        AvgRate(num_quant_ind,SNR_ind) = sum(R(num_quant_ind,:))/nITER;
        AvgRate_decom(num_quant_ind,SNR_ind) = sum(R_decom(num_quant_ind,:))/nITER;
        AvgRate_ortho(num_quant_ind,SNR_ind) = sum(R_ortho(num_quant_ind,:))/nITER;
        AvgRate_fully_zf(num_quant_ind,SNR_ind) = sum(R_fully_zf(num_quant_ind,:))/nITER;    
    end
end
% % ------------------------Sum Rate Calculation (per SNR)---------------------------------
% for num_quant_ind = 1:length(Num_Qbits)       
%     AvgRate_less(num_quant_ind,SNR_ind) = sum(R_less(num_quant_ind,:))/nITER;
%     AvgRate(num_quant_ind,SNR_ind) = sum(R(num_quant_ind,:))/nITER;
%     AvgRate_decom(num_quant_ind,SNR_ind) = sum(R_decom(num_quant_ind,:))/nITER;
%     AvgRate_ortho(num_quant_ind,SNR_ind) = sum(R_ortho(num_quant_ind,:))/nITER;
%     AvgRate_fully_zf(num_quant_ind,SNR_ind) = sum(R_fully_zf(num_quant_ind,:))/nITER;    
% end
%% ================ Result  =================================================== %%
figure(1);
plot(SNR_dBa,AvgRate_less(1,:),'-ro');
hold on
plot(SNR_dBa,AvgRate(1,:),'-bo');
hold on
plot(SNR_dBa,AvgRate_decom(1,:),'-ko');
hold on
plot(SNR_dBa,AvgRate_ortho(1,:),'-k>');
hold on
plot(SNR_dBa,AvgRate_fully_zf(1,:),'-go');
hold on

grid on 
axis([-10 10 0 80]);
xlabel('SNR (dB)','FontSize',12);
ylabel('Spectral Efficiency (bps/ Hz)','FontSize',12);
titleStr=sprintf('UE = %d, Path = %d',Num_UE,Num_paths);
title(titleStr)
legnedStr1=sprintf('Proposed less, N^{RF} = %d', Num_UE-1);
legnedStr2=sprintf('Proposed, N^{RF} = %d', Num_UE);
legnedStr3=sprintf('Proposed decom, N^{RF} = %d', Num_BS_RFchains);
legnedStr4=sprintf('Proposed ortho, N^{RF} = %d', Num_BS_RFchains);
legnedStr5=sprintf('Fully ZF, N^{RF} = %d', Num_BS_Antennas);
legend(legnedStr1,legnedStr2,legnedStr3,legnedStr4,legnedStr5);

figure(2);
plot(SNR_dBa,AvgRate_less(2,:),'-ro');
hold on
plot(SNR_dBa,AvgRate(2,:),'-bo');
hold on
plot(SNR_dBa,AvgRate_decom(2,:),'-ko');
hold on
plot(SNR_dBa,AvgRate_ortho(2,:),'-k>');
hold on
plot(SNR_dBa,AvgRate_fully_zf(2,:),'-go');
hold on
grid on 
axis([-10 10 0 80]);
xlabel('SNR (dB)','FontSize',12);
ylabel('Spectral Efficiency (bps/ Hz)','FontSize',12);
titleStr=sprintf('UE = %d, Path = %d, Quantized',Num_UE,Num_paths);
title(titleStr)
legnedStr1=sprintf('Proposed less, N^{RF} = %d', Num_UE-1);
legnedStr2=sprintf('Proposed, N^{RF} = %d', Num_UE);
legnedStr3=sprintf('Proposed decom, N^{RF} = %d', Num_BS_RFchains);
legnedStr4=sprintf('Proposed ortho, N^{RF} = %d', Num_BS_RFchains);
legnedStr5=sprintf('Fully ZF, N^{RF} = %d', Num_BS_Antennas);
legend(legnedStr1,legnedStr2,legnedStr3,legnedStr4,legnedStr5);

toc
