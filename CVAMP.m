function [X,Pa] = CVAMP(Y,S,gamma_w,lsfc,AMP_option)
%EMVAMP 此处显示有关此函数的摘要
%Expectation-Maximization framework for activity detection and channel
%estimation in IRS-assisted Massive MIMO system, where E-step is
%implemented with vector AMP algorihtm and M-step follow a
%block-coordinate-descent optimization.
%   此处显示详细说明
%The received signal model is: Y=S*X+W
%where S is the sensing matrix, W is the additive white Gaussian noise with
%unknown noise precision $\gamma_w$ $W_{lm}~CN(0,1/\gamma_w)$.
%For E-step, we place a spike and slab prior for each entry in X that
%$x_{nm}~(1-p_n)\delta(x_{nm}+p_n CN(0,v_nm)$ where activity probability of
%$n$th device $p_n$ and variance $v_{nm}$ are unknown.

%% System Size Extraction
[L,M] = size(Y);
[~,K] = size(S);
%% Hyper-parameters Initialization
p = 0.5*ones(K,1);
V = lsfc * ones(1,M);
%% Variable Initialization
Pi = zeros(K,M);
Mu = zeros(K,M);
Sigma = zeros(K,M);
X_hat = zeros(K,M);
X_var = zeros(K,M);
Gamma = ones(M,1)./(L+L*diag(Y'*Y)/norm(S,'fro')^2);
aclist = 1:K;
totlist = 1:K;
%% Algorithm Parameter
MAXITER = 200;
Damp = 0.03;
Threshold = 1e-4;
relative_change = zeros(1,MAXITER);
%% SVD
[Bar_U,Bar_S,Bar_V] = svd(S,'econ');
Rank_S = rank(Bar_S);
Y_tilde = Bar_S \ Bar_U' * Y;
R = S'*Y;
%% Iteration Process
for t=1:MAXITER
    %% E-step
    X_pre = X_hat; % X record for damp
    Gamma_pre = Gamma; % Gamma record for damp

    % Pre-computation for AMP
    prod_temp = ones(K,1);
    for i=1:length(totlist)
        n = totlist(i);
        for m=1:M
            prod_temp(n) = prod_temp(n) * (V(n,m) * Gamma(m) + 1)...
                * exp( - Gamma(m)^2 * V(n,m) * norm(R(n,m),2)^2/(V(n,m)...
                * Gamma(m) + 1));
        end
    end
    prod_temp(prod_temp<1e-6) = 1e-6;

    % AMP iterations across M antennas
    for m=1:M
        alpha_m = 0;
        for i=1:length(aclist)
            n = aclist(i);
            % Run several standard AMP iterations before the proposed
            % method since the mild condition is not satisfied. If M and L
            % are large, this step could be short.
            % Example: M=L=200, t<2; M=L=120, t<3; M=L=100, t<4; M=L=80,
            % t<5; M=40, L=30, t<14.
            if t<(14)
                Pi(n,m) = (1 + ((1-p(n))/p(n)) * (V(n,m) * Gamma(m) + 1)...
                    * exp( - Gamma(m)^2 * V(n,m) * norm(R(n,m),2)^2/(V(n,m)...
                    * Gamma(m) + 1)))^(-1);
            else
                Pi(n,m) = (1 + ((1-p(n))/p(n)) * prod_temp(n))^(-1);
                if Pi(n,m)<1e-8
                    Pi(n,m)=1e-8;
                end
            end
            Mu(n,m) = V(n,m) * Gamma(m)/(V(n,m) * Gamma(m) + 1) * R(n,m); % Mean of channel
            Sigma(n,m) = real(V(n,m)/(V(n,m) * Gamma(m) + 1)); % Variance of channel
            X_hat(n,m) = Damp * X_pre(n,m) + (1-Damp) * Pi(n,m) * Mu(n,m); % Posterior mean
            X_var(n,m) = Pi(n,m) * Sigma(n,m); % Posterior variance
            phi_temp = 1/Pi(n,m);
            omega_temp = 1 + (Gamma(m)^2 * V(n,m) * (phi_temp-1) * norm(R(n,m),2)^2)/((Gamma(m)*V(n,m) + 1) * phi_temp);
            alpha_temp = (Gamma(m)*V(n,m))/(Gamma(m)*V(n,m)+1) * (omega_temp/phi_temp);
            alpha_m = alpha_m + alpha_temp/K;
        end
        R_tilde = (X_hat(:,m) - alpha_m * R(:,m))/(1 - alpha_m);
        gamma_tilde = real(Gamma(m) * (1 - alpha_m)/alpha_m);
        d = gamma_w * (gamma_w * Bar_S .* Bar_S + gamma_tilde * eye(Rank_S))^(-1)...
            * diag(Bar_S .* Bar_S);
        Gamma(m) = Damp * Gamma_pre(m) + (1-Damp) * real(gamma_tilde * mean(d) / (K/Rank_S - mean(d))); % Update Gamma(m)
        R(:,m) = R_tilde + (K/Rank_S) * Bar_V * diag(d/mean(d)) * (Y_tilde(:,m) - Bar_V' * R_tilde);
    end

    %% M-step
    %Update $p_n$
    p = real(mean(Pi,2));
    p=real(p);
    p(p<1e-8) = 1e-8;
    %% Stop criteria
    relative_change(t) = norm(X_hat-X_pre,'fro')^2/norm(X_hat,'fro')^2;
    if t>5 && relative_change(t) < Threshold
        break;
    end
end
fprintf('Method: CVAMP, it %d: relative_change = %g\n', t, relative_change(t));

% Generate output
X = X_hat;
Pa = p;