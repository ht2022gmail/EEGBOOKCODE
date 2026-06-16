clear all; close all; clc;

%% Creatiging artificial data:
Nsamples = 500;
C0 = diag([3 0.01]);
R1 = [cos(pi/6) -sin(pi/6); sin(pi/6) cos(pi/6)];
R2 = [cos(pi/3) -sin(pi/3); sin(pi/3) cos(pi/3)];
C1 = R1*C0*R1';
C2 = R2*C0*R2';
X = [mvnrnd([0;0], C1, Nsamples)' mvnrnd([0;0], C2, Nsamples)'];
X = bsxfun(@minus, X, mean(X,2));
% X - repmat(mean(X,2),1,2*Nsamples);

%% Principal component analysis
[U, S, V] = svd(X, 'econ');
Y = U'*X;

figure(1); 
subplot(2,3,1); hold on;
plot(X(1,1:Nsamples), X(2,1:Nsamples), 'x', 'color', [0.4 0.4 0.4]);
plot(X(1,Nsamples+1:2*Nsamples), X(2,Nsamples+1:2*Nsamples), 'ko');
plot(3*[0 U(1,1)], 3*[0 U(2,1)], 'color', [0.6 0.6 0.6], 'linewidth', 3);
plot(3*[0 U(1,2)], 3*[0 U(2,2)], 'color', [0.6 0.6 0.6], 'linewidth', 3);
axis([-6 6 -6 6]); axis square;
title('Principal Component Analysis');
xlabel('x_1'); ylabel('x_2');

subplot(2,3,4); hold on;
plot(Y(1,1:Nsamples), Y(2,1:Nsamples), 'x', 'color', [0.4 0.4 0.4]);
plot(Y(1,Nsamples+1:2*Nsamples), Y(2,Nsamples+1:2*Nsamples), 'ko');
axis([-4 4 -4 4]); axis square;
title('Principal Component Analysis');
xlabel('PC_1'); ylabel('PC_2');

%% Signal and Reference covariance matrices
CS = X(:,1:Nsamples)*X(:,1:Nsamples)'/(Nsamples-1)-X(:,(Nsamples+1):(2*Nsamples))*X(:,(Nsamples+1):(2*Nsamples))'/(Nsamples-1);
CR = X*X'/(2*Nsamples-1);

%% Generalized eigenvalue decomposition (GED)
[W, Lambda_GED] = eig(CS, CR);
Y = W'*X;

subplot(2,3,2); hold on;
plot(X(1,1:Nsamples), X(2,1:Nsamples), 'x', 'color', [0.4 0.4 0.4]);
plot(X(1,Nsamples+1:2*Nsamples), X(2,Nsamples+1:2*Nsamples), 'ko');
plot(3*[0 W(1,1)], 3*[0 W(2,1)], 'color', [0.6 0.6 0.6], 'linewidth', 3);
plot(3*[0 W(1,2)], 3*[0 W(2,2)], 'color', [0.6 0.6 0.6], 'linewidth', 3);
axis([-6 6 -6 6]); axis square;
title('Generalized Eigendecomposition');
xlabel('x_1'); ylabel('x_2');

subplot(2,3,5); hold on;
plot(Y(1,1:Nsamples), Y(2,1:Nsamples), 'x', 'color', [0.4 0.4 0.4]);
plot(Y(1,Nsamples+1:2*Nsamples), Y(2,Nsamples+1:2*Nsamples), 'ko');
axis([-3 3 -3 3]); axis square;
title('Generalized Eigendecomposition');
xlabel('GED_1'); ylabel('GED_2');

%% Joint diagonalization (JD)
[P, DR] = eig(CR);
[Q, Lambda_JD] = eig(DR^(-1/2)*P'*CS*P*DR^(-1/2));
W = P*DR^(-1/2)*Q;
Y = W'*X;

subplot(2,3,3); hold on;
plot(X(1,1:Nsamples), X(2,1:Nsamples), 'x', 'color', [0.4 0.4 0.4]);
plot(X(1,Nsamples+1:2*Nsamples), X(2,Nsamples+1:2*Nsamples), 'ko');
plot(3*[0 W(1,1)], 3*[0 W(2,1)], 'color', [0.6 0.6 0.6], 'linewidth', 3);
plot(3*[0 W(1,2)], 3*[0 W(2,2)], 'color', [0.6 0.6 0.6], 'linewidth', 3);
axis([-6 6 -6 6]); axis square;
title('Joint Diagonalization');
xlabel('x_1'); ylabel('x_2');

subplot(2,3,6); hold on;
plot(Y(1,1:Nsamples), Y(2,1:Nsamples), 'x', 'color', [0.4 0.4 0.4]);
plot(Y(1,Nsamples+1:2*Nsamples), Y(2,Nsamples+1:2*Nsamples), 'ko');
axis([-3 3 -3 3]); axis square;
title('Joint Diagonalization');
xlabel('JD_1'); ylabel('JD_2');

set(gcf, 'color', 'w');

