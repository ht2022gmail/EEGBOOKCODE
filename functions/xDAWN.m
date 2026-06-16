function [U, Lambda, A] =  xDAWN(X, tidx, tau)

[Nchannels, Ntime] = size(X); % continuous data
Ntrials = length(tidx);

D = zeros(tau, Ntime);
for tt=tidx
    D(:,tt:tt+tau-1) = D(:,tt:tt+tau-1) + eye(tau);
end
A = X*pinv(D);

CS = (A*D)*(A*D)'/Ntime/Ntrials;
CR = X*X'/Ntime/Ntrials;

[U, Lambda] = eig(CS, CR);
Lambda = diag(Lambda);
[Lambda, idx] = sort(Lambda, 'descend');
U = U(:,idx);
