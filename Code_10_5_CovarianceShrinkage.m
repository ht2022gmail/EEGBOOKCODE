clear all; close all; clc;
% Defining the random variables:
D = 50; % dimensions of random variable vector
N = 100;    % # samples
S0 = diag(10.^linspace(2,0,D));

%% Shrinkage estimator:
Niters = 100;
eig_S0 = sort(real(eig(S0)), 'descend');
eig_sample = zeros(D, Niters);
eig_shrinkage = zeros(D, Niters);

for jj=1:Niters
    % multivariate normal distribution:
    X = mvnrnd(zeros(D,1), S0, N)'; 
    % sample mean and covariance:
    [S_shrinkage, S_sample, rho] = shrinkage_covariance(X);
    eig_sample(:,jj) = sort(real(eig(S_sample)), 'descend');
    eig_shrinkage(:,jj) = sort(real(eig(S_shrinkage)), 'descend');
end
    
%% Visualzing the results:
figure(1); clf; hold on
plot(1:D, eig_S0, 'k.-', 'linewidth', 2, 'markersize', 20);
% plot(1:D, mean(eig_sample,2), 'bx-');
errorbar(1:D, mean(eig_sample,2), std(eig_sample,[],2), 'b--', 'linewidth', 2)
% plot(1:D, mean(eig_shrinkage,2), 'r.-');
errorbar(1:D, mean(eig_shrinkage,2), std(eig_shrinkage,[],2), 'r--', 'linewidth', 2)
xlim([0 25])

set(gca, 'plotboxaspectratio', [1. 1 1]);
ti = title(['$p=$', num2str(D),  ', $n=$', num2str(N), ', $n/p=$', num2str(N/D)]);
set(ti, 'interpreter', 'latex', 'fontsize', 16);
xl = xlabel('Order', 'fontsize', 14);
yl = ylabel('Eigenvalues', 'fontsize', 14);
leg = legend('True', 'Sample', 'Shrinkage', 'fontsize', 11);
set(gcf, 'color', 'w');
