clearvars
clc
%% Notation
% ---------------------------------------
% |K      |Number of potential users    |
% |epsilon|Active probability           |
% |L      |Pilot sequence length        |
% |N      |Number of antennas           |
% |S      |Pilot sequence               |
% |gamma  |Active indicator and pathloss|
% |H      |Small-scale fading           |
% |Y      |Received signals             |
% ---------------------------------------
%% Simulation setting

% Total Number of Users in Cell Network
N = 500;

% Number of Antennas at each AP
M = 40;

% Binomial Activizty percentage: Activity pattern
epsilon = 0.05;

% Pilot Sequence Length
L = 30;

% SNR
SNR = 10;

% Noise Parameter :  bandwidth = 1 MHz
sigma_sqr_dBm = -109;
sigma_sqr = 10^((sigma_sqr_dBm-30)/10); % Real noise variance
sigma_sqrN = 1; % Normalized Sigma2

% Transmit Power Constraint in W
TxPow = 200e-3; % 23dBm

% Number of Monte-Carlo simulations
monte = 100;

% Variables for evaluation
G_real      = repmat(zeros(N,M),[1 1 monte]);
Active_List = zeros(N,monte);
%% Generate Signature Sequence (multiplied by transmit power)
S = (1/sqrt(2*L)*complex(randn(L,N),randn(L,N)));
S = S./vecnorm(S,2,1);
%% Estimation Initialization
G_hat_CVAMP    = repmat(zeros(N,M),[1 1 monte]);
G_hat_CAMP     = repmat(zeros(N,M),[1 1 monte]);
P_a_CVAMP      = zeros(N,monte);
P_a_CAMP       = zeros(N,monte);
%% Generate received signals
for i = 1:1:monte
    %     Random active devices (Binomial distribution)
    %     Ac_list = binornd(1,epsilon,[K,1]);
    %     Active_List(:,i) = Ac_list;

    %     Random active devices (fixed number)
    Ac_idx = randperm(N);
    Ac_list = zeros(N,1);
    Ac_list(Ac_idx(1:N*epsilon)) = 1;

    %     Index of active devices
    idx = find(Ac_list);
    Active_List(:,i) = Ac_list;

    %     Pathloss (distance uniform distributed in [0.5,1]km)
    %     Square area
    %     AP = [0,0];
    %     UT=unifrnd(-3/2,3/2,K,2);
    %     dist = sqrt(sum((UT - AP).^2,2));
    %     dist = sqrt(sum(UT .^2,2));
    % dist = 0.05 + 0.45 * rand(1,K);
    % K_rician = 10.^(1.3-3*dist);
    % K_rician = 6*rand(1,K);
    %     Circular area
    %     dist = 1 + 2 * rand(1,K);
    %     Fixed distance
    %     dist = 1 * ones(1,K);
    % Beta = 0.2 * 10.^((-128.1 - 37.6*log10(dist))/10);
    %     Pathloss for normalized noise (variance=1)
    % Beta = Beta/sigma_sqr;%     SNR = -10 * log10(0.2) + 25;
    %     Beta = 10^(SNR/10) * ones(1,K);
    %     Calculate gamma_n=a_n*sqrt(beta_n)

    % idx_snr = find(Beta>10);
    % Beta(idx_snr) = 10;
    % idx_snr2 = find(Beta<10);
    % Beta(idx_snr2) = 10;

    Beta = 10^(SNR/10)*ones(1,N);
    gamma(idx,i) = Beta(idx);

    % Received signal
    H = (1/sqrt(2)*complex(randn(N,M),randn(N,M)));
    % H_smv = reshape(H,1,N*K);
    % H_smv = H_smv.';

    % Rician fading
    % N_rician = randperm(K);
    % Percent_rician = 0;
    % random_phase_shift = sqrt(-1)*2*pi*rand(K*Percent_rician,1);
    % random_phase_shift = zeros(K*Percent_rician,1);
    % H(N_rician(1:Percent_rician*K),:) = diag(sqrt(K_rician(N_rician(1:Percent_rician*K))./(K_rician(N_rician(1:Percent_rician*K))+1)))*exp(random_phase_shift * linspace(0,N-1,N)) + diag(sqrt(1./(K_rician(N_rician(1:Percent_rician*K))+1))) * H(N_rician(1:Percent_rician*K),:);

    W = (sqrt(sigma_sqrN))*(1/sqrt(2)*complex(randn(L,M),randn(L,M)));
    X = diag(sqrt(L*(gamma(:,i))));
    Y_real = S*X*H;
    Y = Y_real + W;
    % X_smv = repmat(diag(X),[1 N]);
    % Y_smv = S_smv*diag(reshape(X_smv,[N*K,1]))*H_smv;
    % Y_smv = Y_smv + reshape(W,L*N,1);
    % G_smv_real(:,i) = diag(reshape(X_smv,[N*K,1]))*H_smv;
    % Record G
    G_real(:,:,i) = X * H;
    %     Record pathloss
    pathloss = Beta';
    %     pathloss_real = Beta_noised';

    fprintf('Set %d\n', i);
    %% Activity detection and Channel estimation
    % AMP
    clear optEM
    optEM.SNRdB = SNR;
    optEM.lambda = 1;
    optEM.active_mean = 0; %Initialize at correct active mean
    % optEM.active_weights = 1; %Trivial weights since 1-term mixture
    % optEM.active_var = 1; %Initialize correct active variance
    optEM.noise_var = 1; %Initialize at correct noise variance

    % Turn off EM learning of all but the sparsity-rate, lambda
    optEM.learn_lambda= false;
    optEM.learn_mean = false;
    optEM.learn_weights = false;
    optEM.learn_var = false;
    optEM.learn_noisevar = false;

    % Run, time, and check performance of EM
    % BGAMP
    optEM.heavy_tailed = false; % since operating on a sparse signal
    % [G_hat_AMP(:,:,i), EMfin] = EMBGAMP(Y,S,optEM,[]);
    % [G_hat_AMP(:,:,i),conv_AMP(i)] = AMP(Y,S,1,Beta'*L,epsilon);
    % [~,~,~,~,~,conv_AMP(i)] = EMBGAMP(Y,S,optEM,[]);

    % row AMP
    % [~,G_hat_rowAMP(:,:,i),~,~,~,conv_rowAMP(i)] = noisyCAMPmmseforKLS(S,K,N,L,Y,G_real(:,:,i),100,0.5,Beta'*L,1,1);

    % CAMP
    % [G_hat_CAMP(:,:,i),P_a_CAMP(:,i),conv_CAMP(i),recov_it] = CAMP(Y,S,1,Beta'*L,epsilon);

    % GMMV-AMP
    % [G_hat_gmmvamp(:,:,i), lambda, conv_gmmvamp] = gmmv_amp(Y, S, 0.3, 30, 1e-2, 0);
    % P_a_CAMP(:,i) = mean(lambda,2);

    % EM-AMP
    % [G_hat_EMAMP(:,:,i),P_a_EMAMP(:,i),conv_EMAMP(i)] = EMAMP(Y,S,sigma_sqrN,Beta'*L,epsilon*ones(K,1));
    % [G_hat_EMAMP(:,:,i),P_a_EMAMP(:,i),conv_EMAMP(i)] = EMAMP(Y,S,sigma_sqrN,Beta'*L,Ac_list);

    % BO-VAMP
    % [G_hat_BOVAMP(:,:,i)] = ADVAMP(Y,S,0.1,epsilon*ones(K,1),L*Beta');
    % [G_hat_EMBOVAMP(:,:,i),P_a_EMBOVAMP(:,i),conv_EMVAMP(i)] = EMVAMP(Y,S,Beta'*L,epsilon*ones(K,1));
    % [G_hat_EMBOVAMP(:,:,i),P_a_EMBOVAMP(:,i),conv_EMVAMP(i)] = EMVAMP(Y,S,Beta'*L,Ac_list);
    % [G_hat_EMBOVAMP_decoupled(:,:,i),~] = EMVAMP_decoupled(Y,S);
    [G_hat_CAMP(:,:,i),P_a_CAMP(:,i),conv_CAMP(i),~,NMSECVAMP(:,:,i)] = CVAMP(Y,S,sigma_sqrN,Beta'*L,G_real(:,:,i));
    % if sum(P_a_CAMP(:,i))>50
    %     break;
    % end

    % [G_hat_CVAMP_smv(:,i),~,conv_CAMP(i),~,NMSEsmvCVAMP(:,:,i)] = smvCVAMP(Y_smv,S_smv,sigma_sqrN,L*Beta(1)'*ones(K*N,1),epsilon,G_smv_real(:,i));
    % [G_hat_CAMP(:,:,i),P_a_CAMP(:,i),conv_CAMP(i),recov_it] = CVAMP_independentS(Y,S,sigma_sqrN,Beta'*L,epsilon);
    % [G_hat_CVAMP(:,:,i),P_a_CVAMP(:,i),conv_CVAMP(i)] = OracleCVAMP(Y,S,sigma_sqrN,Beta'*L,Ac_list);
    % [G_hat_CVAMP(:,:,i),P_a_CVAMP(:,i),conv_CVAMP(i)] = BSVAMP(Y,S,sigma_sqrN,Beta'*L,epsilon);

    % cov
    % [gamma_hat(:,i)] = su_decode_activity_pattern(sigma_sqrN,Beta',Y,S);

    % GHVI
    % [G_hat_GHVI(:,:,i),z_GHVI(:,i),~] = VIAD_GH_cellfree(Y,S);
    % [G_hat_GHVI(:,:,i),~] = VIAD_GG_cellfree(Y,S);

    % CE
    % [CE_hat(:,:,i)] = LSCE(Y,S);
    % MSE_CVAMP=MSE_CVAMP+norm(G_real(:,:,i)-G_hat_CVAMP(:,:,i),'fro')^2/monte;
    % MSE_CAMP=MSE_CAMP+norm(G_real(:,:,i)-G_hat_AMP(:,:,i),'fro')^2/monte;
end
%% END
fprintf('Simulation Finished\n');
% NMSE(G_hat_CVAMP,G_real)
% NMSE(G_hat_AMP,G_real)
% NMSE(CE_hat,G_real)
% x=P_a_CAMP;
% x(Active_List==0)=0;
% y=sum(x,1);
% ac=mean(y)
% x=P_a_CAMP;
% x(Active_List==1)=0;
% x(x<=1e-8)=0;
% y=sum(x);
% inac=mean(y)
% 25-ac+inac