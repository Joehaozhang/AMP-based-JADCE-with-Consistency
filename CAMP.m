function [X,Pa] = CAMP(Y,S,gamma_w,lsfc)
%% System Size Extraction
[L,M] = size(Y);
[~,N] = size(S);
%% Hyper-parameters Initialization
p = 0.5*ones(N,1);
V = lsfc*ones(1,M);
%% Variable Initialization
Pi = zeros(N,M);
Mu = zeros(N,M);
Sigma = zeros(N,M);
X_hat = zeros(N,M);
X_var = zeros(N,M);
Gamma = ones(M,1)./(L+2*L*diag(Y'*Y)/norm(S,'fro')^2);
% Gamma = 1/(50*L) * ones(M,1);
R = S'*Y;
Vk = zeros(L,M);
aclist = 1:N;
%% Algorithm Parameter
T_max = 200;
recov_it = zeros(T_max,1);
Damp = 0.03;
Convergence_thr = 1e-4;
normalized_change = zeros(1,T_max);
state_conv = 0;
%% Iteration Process
for t=1:T_max
    X_pre = X_hat;
    Gamma_pre = Gamma;
    sum_temp = ones(N,1);
    for i=1:length(aclist)
        n = aclist(i);
        for m=1:M
            sum_temp(n) = sum_temp(n) * (V(n,m) * Gamma(m) + 1)...
                * exp( - Gamma(m)^2 * V(n,m) * norm(R(n,m),2)^2/(V(n,m)...
                * Gamma(m) + 1));
        end
    end

    sum_temp(sum_temp<1e-2) = 1e-6;
    for m=1:M
        alpha_m = 0;
        for i=1:length(aclist)
            n = aclist(i);
            if t<(14)
                Pi(n,m) = (1 + ((1-p(n))/p(n)) * (V(n,m) * Gamma(m) + 1)...
                * exp( - Gamma(m)^2 * V(n,m) * norm(R(n,m),2)^2/(V(n,m)...
                * Gamma(m) + 1)))^(-1); % Posterior activity probability
            else
                Pi(n,m) = (1 + ((1-p(n))/p(n)) * sum_temp(n))^(-1);
                if Pi(n,m)<1e-8
                    Pi(n,m)=1e-8;
                end
            end
            Mu(n,m) = V(n,m) * Gamma(m)/(V(n,m) * Gamma(m) + 1) * R(n,m);
            Sigma(n,m) = real(V(n,m)/(V(n,m) * Gamma(m) + 1));
            X_hat(n,m) = Damp * X_pre(n,m) + (1-Damp) * Pi(n,m) * Mu(n,m);
            X_var(n,m) = Pi(n,m) * Sigma(n,m);
            phi_temp = 1/Pi(n,m);
            omega_temp = 1 + (Gamma(m)^2 * V(n,m) * (phi_temp-1) * norm(R(n,m),2)^2)/((Gamma(m)*V(n,m) + 1) * phi_temp);
            alpha_temp = (Gamma(m)*V(n,m))/(Gamma(m)*V(n,m)+1) * (omega_temp/phi_temp);
            alpha_m = alpha_m + alpha_temp/N;
        end
        Vk(:,m) = Y(:,m) - S*X_hat(:,m) + (N/L)*alpha_m*Vk(:,m);
        R(:,m) = X_hat(:,m) + S'*Vk(:,m);
        Gamma_temp = (norm(Vk(:,m),2)^2/L)^(-1);
        % Gamma_temp = real(1/(1/gamma_w+(N/L)*(mean(Mu(:,m).*conj(Mu(:,m)) + X_var(:,m)))));
        % Gamma_temp = real(1/(1/gamma_w+(N/L)*(mean(X_var(:,m)))));
        Gamma(m) = Damp * Gamma_pre(m) + (1-Damp) * Gamma_temp;
    end


    % M-step

    %Update $p_n$
    p = real(mean(Pi,2));
    p=real(p);
    p(p<1e-8) = 1e-8;
    % aclist = find(p>1e-8);
    % recov_it(t) = length(find(p>(1-1e-8)));
    % 
    % aclist = find(p>1e-8);
    % inaclist = p==1e-8;
    % X_hat(inaclist,:)=0;

    % Stop criteria
    normalized_change(t) = norm(X_hat-X_pre,'fro')^2/norm(X_hat,'fro')^2;
    if normalized_change(t) < Convergence_thr
        state_conv = 1;
        break;
    end
end
fprintf('Method: CAMP, it %d: Conv_state = %d, relative_change = %g\n', t, state_conv, normalized_change(t));

X = X_hat;
Pa = p;