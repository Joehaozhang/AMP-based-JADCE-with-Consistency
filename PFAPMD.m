function [PFA_PMD] = PFAPMD(z_hat,Active_List,N_thr)
%% Notation
% ----------------------------------------
% |G_hat       |Channel Estimation       |
% |Active_List |Active device list       |
% |N_thr       |Number of thresholds     |
% ----------------------------------------
%% Initializations
[N,monte] = size(z_hat);
N_total = N * monte;
N_active      = length(find(Active_List > 0));
Thr_min       = 1/1e2;
Thr_max       = 1;
threshold     = linspace(Thr_min, Thr_max, N_thr);
PFA           = zeros(N_thr,1);
PMD           = zeros(N_thr,1);
%% PMDPFA
for i=1:N_thr
    MD    = 0;
    FA    = 0;
    Thr_i = threshold(i);
    for j=1:monte
        Idx_hat  = find(z_hat(:,j) > Thr_i);
        Idx_real = find(Active_List(:,j) > 0);
        if isempty(Idx_real)
            break;
        end
        MD       = MD + length(setdiff(Idx_real,Idx_hat));
        FA       = FA + length(setdiff(Idx_hat,Idx_real));
    end
    Pmd = MD/N_active;
    Pfa = FA/(N_total - N_active);

    PMD(i) = Pmd;
    PFA(i) = Pfa;
end

PFA_PMD = [PFA PMD];
end