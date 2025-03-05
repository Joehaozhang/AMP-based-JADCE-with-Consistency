function [X,Pa] = CAMP_JADCE(Y,S,gamma_w,lsfc,AMP_option)
% CAMP_JADCE
% AMP-based Activity Detection and Channel Estimation with Consistency
% Author: Hao Zhang
%
% Last updated: 2025/03/05
% -----------------------------------------------------------------------
%                               INPUTs                                  |
% Y: input matrix (received signals in APs)                             |
% S: pilot signals transmitted                                          |
% gamma_w: noise variance                                               |
% lsfc: large-scale fading coefficients                                 |
% AMP_option: 'AMP' or 'vector AMP' algorithm for message passing       |
% ----------------------------------------------------------------------|

% Determine single-cell or cell-free
[L,M,K] = size(Y);

% Notice that the difference between single-cell and cell-free algorithm
% is number of variables and loops
if K==1
    % single-cell algorithm
    [X,Pa] = AMP_singlecell(Y,S,gamma_w,lsfc,AMP_option);
elseif K>1
    % cell-free algorithm
    [X,Pa] = AMP_cellfree(Y,S,gamma_w,lsfc,AMP_option);
end