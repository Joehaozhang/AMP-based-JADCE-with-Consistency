function [X,Pa] = CVAMP_cellfree(Y,S,gamma_w,lsfc)
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
[L,M,K] = size(Y);
[~,N] = size(S);
%% Hyper-parameters Initialization
p = 0.5*ones(N,1);
V = repmat( zeros(N,M), [1 1 K]);
for k=1:K
    V(:,:,k) = lsfc(:,k)*ones(1,M);
end
%% Variable Initialization
Pi = repmat(zeros(N,M), [1 1 K]);
Mu = repmat(zeros(N,M), [1 1 K]);
Sigma = repmat(zeros(N,M), [1 1 K]);
X_hat = repmat(zeros(N,M), [1 1 K]);
X_var = repmat(zeros(N,M), [1 1 K]);
Gamma = zeros(M,K);
R     = repmat(zeros(N,M), [1 1 K]);
for k=1:K
    Gamma(:,k) = ones(M,1)./(L+2*L*diag(Y(:,:,k)'*Y(:,:,k))/norm(S,'fro')^2);
    R(:,:,k)   = S'*Y(:,:,k);
end
aclist = 1:N;
%% Algorithm Parameter
MAXITER = 200;
Damp = 0;
Threshold = 1e-4;
relative_change = zeros(1,MAXITER);
%% SVD
[Bar_U,Bar_S,Bar_V] = svd(S,'econ');
Rank_S = rank(Bar_S);
Y_tilde = repmat(zeros(L,M),[1 1 K]);
for k=1:K
    Y_tilde(:,:,k) = Bar_S \ Bar_U' * Y(:,:,k);
    R(:,:,k)   = S'*Y(:,:,k);
end
%% Iteration Process
for t=1:MAXITER
    %% E-step
    X_pre = X_hat; % X record for damp
    Gamma_pre = Gamma; % Gamma record for damp
    if t>0
        Damp = 0.03; % Damp hyper-parameter
    end

    % Pre-computation for AMP
    sum_temp = ones(N,1);
    for i=1:length(aclist)
        n = aclist(i);
        for k=1:K
            for m=1:M
                sum_temp(n) = sum_temp(n) * (V(n,m,k) * Gamma(m,k) + 1)...
                    * exp( - Gamma(m,k)^2 * V(n,m,k) * norm(R(n,m,k),2)^2/(V(n,m,k)...
                    * Gamma(m,k) + 1));
            end
        end
    end
    sum_temp(sum_temp<1e-6) = 1e-6;

    % AMP iterations across M antennas
    for k=1:K
        for m=1:M
            alpha_m = 0;
            for i=1:length(aclist)
                n = aclist(i);
                if t<(2)
                    Pi(n,m,k) = (1 + ((1-p(n))/p(n)) * (V(n,m,k) * Gamma(m,k) + 1)...
                        * exp( - Gamma(m,k)^2 * V(n,m,k) * norm(R(n,m,k),2)^2/(V(n,m,k)...
                        * Gamma(m,k) + 1)))^(-1); % Posterior activity probability
                else
                    Pi(n,m,k) = (1 + ((1-p(n))/p(n)) * sum_temp(n))^(-1);
                    if Pi(n,m)<1e-8
                        Pi(n,m)=1e-8;
                    end
                end
                Mu(n,m,k) = V(n,m,k) * Gamma(m,k)/(V(n,m,k) * Gamma(m,k) + 1) * R(n,m,k); % Mean of channel
                Sigma(n,m,k) = real(V(n,m,k)/(V(n,m,k) * Gamma(m,k) + 1)); % Variance of channel
                X_hat(n,m,k) = Damp * X_pre(n,m,k) + (1-Damp) * Pi(n,m,k) * Mu(n,m,k); % Posterior mean
                X_var(n,m,k) = Pi(n,m,k) * Sigma(n,m,k); % Posterior variance
                phi_temp = 1/Pi(n,m,k);
                omega_temp = 1 + (Gamma(m,k)^2 * V(n,m,k) * (phi_temp-1) * norm(R(n,m,k),2)^2)/((Gamma(m,k)*V(n,m,k) + 1) * phi_temp);
                alpha_temp = (Gamma(m,k)*V(n,m,k))/(Gamma(m,k)*V(n,m,k)+1) * (omega_temp/phi_temp);
                alpha_m = alpha_m + alpha_temp/N;
            end
            R_tilde = (X_hat(:,m,k) - alpha_m * R(:,m,k))/(1 - alpha_m);
            gamma_tilde = real(Gamma(m,k) * (1 - alpha_m)/alpha_m);
            d = gamma_w * (gamma_w * Bar_S .* Bar_S + gamma_tilde * eye(Rank_S))^(-1)...
                * diag(Bar_S .* Bar_S);
            Gamma(m,k) = Damp * Gamma_pre(m,k) + (1-Damp) * real(gamma_tilde * mean(d) / (N/Rank_S - mean(d))); % Update Gamma(m,k)
            R(:,m,k) = R_tilde + (N/Rank_S) * Bar_V * diag(d/mean(d)) * (Y_tilde(:,m,k) - Bar_V' * R_tilde);
        end
    end

    %% M-step
    %Update $p_n$
    p = real(mean(Pi,[2 3]));
    p=real(p);
    p(p<1e-8) = 1e-8;

    if (t>0)
        aclist = find(p>1e-8);
        inaclist = p==1e-8;
        X_hat(inaclist,:)=0;
    end
    %% Stop criteria
    relative_change(t) = norm(X_hat-X_pre,'fro')^2/norm(X_hat,'fro')^2;
    if relative_change(t) < Threshold
        break;
    end
end
fprintf('Method: CVAMP, it %d: relative_change = %g\n', t, relative_change(t));
%% Generate output
X = X_hat;
Pa = p;