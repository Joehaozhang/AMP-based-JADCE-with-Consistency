clear all
clc
%% Notation
% ----------------------------------------
% |D       |Area size                    |
% |K       |Number of AP                 |
% |N       |Number of potential users    |
% |epsilon |Active ratio                 |
% |L       |Pilot sequence length        |
% |M       |Number of antennas           |
% |S       |Pilot sequence               |
% |gamma   |Active indicator and pathloss|
% |H       |Small-scale fading           |
% |Y       |Received signals             |
% ----------------------------------------
%% Simulation setting
% Area Size(km)
D = 3;

% Number of BS
K = 1;

% Total Number of Users in Cell Network
N = 500;

% Number of Antennas at each AP
M = 40;

% Binomial Activity percentage: Activity pattern
epsilon = 0.05;

% Pilot Sequence Length
L = 30;

% Noise Parameter : bandwidth = 1 MHz
sigma_sqr_dBm = -109;
sigma_sqr = 10^((sigma_sqr_dBm-30)/10); % Real noise variance
sigma_sqrN = 1; % Normalized Sigma2

% Transmit Power Constraint in W
TxPow = 200e-3; % 23dBm

% Number of Monte-Carlo simulations
monte = 10;
%% Generate Signature Sequence (multiplied by transmit power)
S = (1/sqrt(2*L)*complex(randn(L,N),randn(L,N)));
S = S./vecnorm(S,2,1);
%%    Locations of BS and Users
%     Uniformly distributed APs
AP = [-D/3,3*D/8;0,3*D/8;D/3,3*D/8;-D/3,D/8;0,D/8;D/3,D/8;-D/3,-D/8;0,-D/8;D/3,-D/8;-D/3,-3*D/8;0,-3*D/8;D/3,-3*D/8];

%     Random APs
% AP=unifrnd(-D/2,D/2,K,2);
%% Variables for evaluation initialization
Y_real      = zeros(L,M,K);
Y           = Y_real;
gamma       = zeros(N,monte,K);
G_real      = repmat(zeros(N,M),[1 1 M monte]);
Active_List = zeros(N,monte);
Beta        = zeros(N,K);
Beta_true   = zeros(N,K);
Beta2       = zeros(N,K);
%% Estimation Initialization
G_hat_ghvi     = repmat(zeros(N,M),[1 1 K monte]);
G_hat_AMP      = repmat(zeros(N,M),[1 1 K monte]);
z_hat_ghvi     = zeros(N,monte);
Pa_AMP      = zeros(N,monte);
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
    Active_List(idx,i)=1;
    %%     Uniformly allocate the UT locations in a DxD area
    UT=unifrnd(-D/2,D/2,N,2);
    %%     Uniformly allocate the UT locations in a circular area
    for k=1:K
        %     Distance-based pathloss (distance from M BS)
        dist = sqrt(sum((UT - AP(k,:)).^2,2));
        Beta(:,k) = 0.2 * 10.^((-128.1 - 37.6*log10(dist))/10);

        %     Pathloss with noise
        Beta_true(:,k) = 0.2 * 10.^((-128.1 - 37.6*log10(dist)+0*rand(N,1))/10);

        %     Pathloss for normalized noise (variance=1)
        Beta(:,k) = Beta(:,k)/sigma_sqr;
    end

    maxBeta           = max(Beta,[],2);
    maxBeta2          = maxBeta;
    idx_snr           = find(maxBeta > 10);
    maxBeta2(idx_snr) = 10;
    idx_snr2          = find(maxBeta < 10);
    maxBeta2(idx_snr2)= 10;
    PowControl        = maxBeta2./maxBeta;
    Beta2             = PowControl.*Beta_true/sigma_sqr;
    idx2              = find(maxBeta2>0);
    idx               = intersect(idx,idx2);

    for k=1:K
        gamma(idx,i,k) = Beta2(idx,k);

        % NLoS channel
        H = (1/sqrt(2)*complex(randn(N,M),randn(N,M)));

        % Noise
        W = (1/sqrt(2)*complex(randn(L,M),randn(L,M)));

        % Large-scale fading coefficients
        X = diag(sqrt(L*gamma(:,i,k)));

        %     Record pathloss
        pathloss = Beta2';

        % Record true channel
        G_real(:,:,k,i) = X * H;

        % True signal
        Y_real(:,:,k) = S*X*H;

        % Observed signal
        Y(:,:,k) = Y_real(:,:,k) + W;
    end
    fprintf('Trial %d\n', i);
    %% Activity detection and Channel estimation
    [G_hat_AMP(:,:,:,i), Pa_AMP(:,i)] = CVAMP(Y,S,1,Beta2*L);
end
%% Estimation END
fprintf('Simulation Finished\n');
%% Optional step: Dominant channel-based detection
% Dominant AP selection for devices
DominantAPSelection();