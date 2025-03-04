clearvars
clc
%% Notation
% ---------------------------------------
% |N      |Number of potential users    |
% |epsilon|Active probability           |
% |L      |Pilot sequence length        |
% |N      |Number of antennas           |
% |S      |Pilot sequence               |
% |gamma  |Active indicator and pathloss|
% |H      |Small-scale fading           |
% |Y      |Received signals             |
% ---------------------------------------
%% Simulation setting
% Area Size(km)
D = 1;

% Total Number of Users in Cell Network
N = 500;

% Number of Antennas at each AP
M = 40;

% Binomial Activizty percentage: Activity pattern
epsilon = 0.05;

% Pilot Sequence Length
L = 30;

% Noise Parameter :  bandwidth = 1 MHz
sigma_sqr_dBm = -109;
sigma_sqr = 10^((sigma_sqr_dBm-30)/10); % Real noise variance
sigma_sqrN = 1; % Normalized Sigma2

% Transmit Power Constraint in W
TxPow = 200e-3; % 23dBm

% Number of Monte-Carlo simulations
monte = 10;

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
    AP = [0,0];
    UT=unifrnd(-D/2,3/2,N,2);
    dist = sqrt(sum((UT - AP).^2,2));
    Beta = 0.2 * 10.^((-128.1 - 37.6*log10(dist))/10);

    %     Pathloss for normalized noise (variance=1)
    Beta = Beta/sigma_sqr;

    % SNR
    SNR = 10;

    % Power allocation
    Beta = 10^(SNR/10) * ones(1,N);

    % Combine activity indicator and large-scale fading
    gamma = zeros(N,1);
    gamma(idx,i) = Beta(idx);

    % Device channel
    H = (1/sqrt(2)*complex(randn(N,M),randn(N,M)));

    % Noise
    W = (sqrt(sigma_sqrN))*(1/sqrt(2)*complex(randn(L,M),randn(L,M)));

    % Large-scale fading coefficients
    X = diag(sqrt(L*(gamma(:,i))));


    % Record true channel
    G_real(:,:,i) = X * H;

    % True signal
    Y_real = S*X*H;

    % Observed signal
    Y = Y_real + W;

    fprintf('Trial %d\n', i);

    [G_hat_CAMP(:,:,i),P_a_CAMP(:,i)] = AMP_with_consistent_sparsity(Y,S,sigma_sqrN,Beta'*L,'AMP');
end
%% END
fprintf('Simulation Finished\n');