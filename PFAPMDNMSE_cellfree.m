function [PFA_PMD_PD_MSE] = PFAPMDNMSE_cellfree(G_hat,Active_List,N_thr,Ga_hat,Ga_real,Ganorm_real,Ganorm_hat)

%% Notation
% ----------------------------------------
% |G_hat       |Channel Estimation       |
% |Active_List |Active device list       |
% |N_thr       |Number of thresholds     |
% ----------------------------------------
%% Initialization
[N,M,K,monte] = size(G_hat);
N_total       = N * monte;
N_active      = length(find(Active_List > 0));
Thr_min       = 1/1e2;
Thr_max       = 10 * Thr_min;
threshold     = linspace(Thr_min, Thr_max, N_thr);
PFA           = zeros(N_thr,1);
PMD           = zeros(N_thr,1);
PD            = zeros(N_thr,1);
MSE           = zeros(N_thr,1);
%% PFAPMDNMSE
for i=1:N_thr
    MD      = 0;
    FA      = 0;
    MSE_hat = 0;
    Thr_i   = threshold(i) * max(Ganorm_real);
    for j=1:monte
        Idx_hat  = find(Ganorm_hat(:,1,j) > Thr_i(:,:,j));
        Idx_real = find(Active_List(:,j) > 0);
        if isempty(Idx_real)
            break;
        end
        MD       = MD + length(setdiff(Idx_real,Idx_hat));
        FA       = FA + length(setdiff(Idx_hat,Idx_real));
        G_hat    = Ga_hat(:,:,j);
        G_hat(setdiff(1:N,Idx_real),:) = 0;
        G_hat(setdiff(1:N,Idx_hat),:) = 0;
        G_real   = Ga_real(:,:,j);
        G_hat    = G_hat(Idx_real,:);
        G_real   = G_real(Idx_real,:);
        Err      = G_hat - G_real;
        MSE_new  = diag(Err * Err')./diag(G_real * G_real');
        MSE_hat  = MSE_hat + sum(MSE_new)/length(Idx_real);
    end
    Pmd = MD/N_active;
    Pfa = FA/(N_total - N_active);
    Pd  = 1 - Pmd;

    PD(i)  = Pd;
    PMD(i) = Pmd;
    PFA(i) = Pfa;
    MSE(i) = MSE_hat/monte;
end
PFA_PMD_PD_MSE = [PFA PMD PD MSE];