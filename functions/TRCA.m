function [W, Lambda, Y, Xb, S, Q, U, V] = TRCA(X, t, tau)

% X: continuous data (N channels x T samplings)
% t: timing vector (K dims)
% tau: block duration (scalar)

[N, T] = size(X);
K = length(t);

% normalization:
X = bsxfun(@minus, X, mean(X,2));
X = bsxfun(@rdivide, X, std(X,[],2));

Xb = zeros(N, tau, K);
for k=1:K
    Xb(:,:,k) = bsxfun(@minus, X(:,t(k):t(k)+tau-1), ...
        mean(X(:,t(k):t(k)+tau-1),2)); 
end

% computation of S matrix:
U = sum(Xb,3)/K;
V = zeros(N,N);
for k=1:K
    V = V+Xb(:,:,k)*Xb(:,:,k)'/K;
end
S = K/(K-1)/tau*(U*U'-V/K);
% S = K/(K-1)/tau*U*U';
Q = X*X'/T;

% eigendecomposition:
[W, Lambda] = eig(S, Q); 
Lambda = diag(Lambda);
[Lambda, index] = sort(Lambda, 'descend');
W = W(:,index);

% sign 
Y = W*X;
for nn=1:N
    y = W(:,nn)'*X;
    W(:,nn) = sign(y*X(7,:)')*W(:,nn);
end

