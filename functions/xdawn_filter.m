function [U, A, S, Lambda] = xdawn_filter(X, onsets, tau)

[Nchannels, Ntime] = size(X);
Ntrials = sum(onsets==1);
onsets_idx = find(onsets==1); 

% Computing ERPs:
D = zeros(tau, Ntime);
for kk=1:Ntrials 
    idx = onsets_idx(kk):(onsets_idx(kk)+tau-1);
    D(:, idx) = D(:, idx) + eye(tau);
end
A = X*pinv(D);

% Computing xDAWN spatial filters:
[U, Lambda] = eig((A*D)*(A*D)'/Ntime, X*X'/Ntime);
Lambda = diag(Lambda);
[Lambda, idx] = sort(Lambda, 'descend');
U = U(:, idx);
S = U'*A;