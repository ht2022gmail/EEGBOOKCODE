function Xepoch = epoch(X, t0, tau)

Nchannels = size(X, 1);
Nepochs = length(t0);
Xepoch = zeros(Nchannels, tau, Nepochs);
for kk=1:Nepochs
    idx = t0(kk):(t0(kk)+tau-1);
    Xepoch(:,:,kk) = X(:,idx);
end