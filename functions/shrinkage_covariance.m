function [S_shrinkage, S_sample, rho] = shrinkage_covariance(X)

% X: data matrix D*N (D: input dimensions, N: number of samples)
% Blankertz, B., Lemm, S., Treder, M., Haufe, S., & Müller, K. R. (2011). 
% Single-trial analysis and classification of ERP components—a tutorial. NeuroImage, 56(2), 814-825.
% Höhne, J., Bartz, D., Hebart, M. N., Müller, K. R., & Blankertz, B. (2016). 
% Analyzing neuroimaging data with subclasses: A shrinkage approach. NeuroImage, 124, 740-751.
%

[D, N] = size(X); 

% compute sample mean and covariance:
mu_sample = mean(X, 2);
X = X-repmat(mu_sample,1,N);
S_sample = 1/(N-1)*X*X';

% compute optimal shirinkage parameter:
nu = trace(S_sample)/D;
numer = 0;
for k=1:N
    numer = numer + norm(X(:,k)*X(:,k)'-S_sample, 'fro')^2;
end
denom = norm(S_sample-nu*eye(D), 'fro')^2;

rho= 1/N/(N-1)*numer/denom;   % Eq. (4) of Hohne et al. (2016)
% rho= 1/(N*(N-1))*numer/denom;   % Eq. (13) of Blankertz et al. (2011)

% compute the shrunk covariance matrix:
S_shrinkage = (1-rho)*S_sample + rho*nu*eye(D);
    

