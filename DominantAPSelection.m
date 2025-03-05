[N,M,K,monte]       = size(G_hat_AMP);
G_hat_dominant_AMP  = zeros(N,M,monte);
G_real_dominant     = zeros(N,M,monte);
DominantAP          = zeros(N,monte);

for j=1:monte
    Gnorm2   = G_real(:,:,:,j).*conj(G_real(:,:,:,j));
    norm2sum = sum(Gnorm2,2);
    [~,DominantAP(:,j)] = max(norm2sum,[],3);
end

for j=1:monte
    for n=1:N
        G_hat_dominant_AMP(n,:,j) = G_hat_AMP(n,:,DominantAP(n,j),j);
        G_real_dominant(n,:,j)    = G_real(n,:,DominantAP(n,j),j);
    end
end

Gnorm2sum_real     = sum(G_real_dominant .* conj(G_real_dominant),2);
Gnorm2sum_hat_AMP  = sum(G_hat_dominant_AMP .* conj(G_hat_dominant_AMP),2);