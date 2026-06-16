%% Generating spiral data that cannot be separated linearly
clear all; close all; clc;
r1 = linspace(0.1,5,200);
theta1 = linspace(0, 6*pi, 200);
X = [[r1.*cos(theta1); r1.*sin(theta1)] [r1.*cos(theta1+pi); r1.*sin(theta1+pi)]]';
Y = [ones(200,1); -ones(200,1)];

N = size(X, 1); N1 = N/2; N2 = N/2;

Ngrid = 200;

figure(1); clf; hold on;
plot(X(1:200,1), X(1:200,2), 'bx');
plot(X(201:400,1), X(201:400,2), 'ro');
set(gca, 'plotboxaspectratio', [1 1 1]); % axis off;
set(gcf, 'color', 'w');

figure(2); figure(3); figure(4);
pause;

%% LDA:
mu1 = mean(X(1:N1,:))'; mu2 = mean(X((N1+1):end,:))'; 
C1 = cov(X(1:N1,:)); C2 = cov(X((N1+1):end,:));
w = inv(C1+C2)*(mu1-mu2);
b = (mu1+mu2)'*w/2;

figure(2); clf; hold on;
y0 = -(w(1)*(-5)-b)/w(2);
y1 = -(w(1)*(-5)-b)/w(2);
ff1 = fill([-5 5 5 -5 -5], [y0 y1 5 5 y0], [1.0 0.8 0.8]);
ff2 = fill([-5 5 5 -5 -5], [y0 y1 -5 -5 y0], [0.8 0.8 1.0]);
set(ff1, 'edgecolor', [1 1 1]); set(ff2, 'edgecolor', [1 1 1]);
plot(X(1:N1, 1), X(1:N2, 2), 'bx'); 
plot(X((N1+1):end, 1), X((N1+1):end, 2), 'ro');
% x0 = linspace(-5,5,100);
% plot(x0, -(w(1)*x0-b)/w(2), 'k-');
xlabel('$x$', 'Interpreter','latex', 'FontSize',12); 
ylabel('$y$', 'Interpreter','latex', 'FontSize',12);
box on;
set(gca, 'plotboxaspectratio', [1 1 1], 'xtick', [], 'ytick', []);
set(gcf, 'color', 'w');


%% Kernel LDA:
disp('Computing Kernel LDA ...'); tic;
% Computing Gram matrix
K = zeros(N, N);
% sigma = 0.2^2;
sigma = 0.3^2;
for ii=1:N
    for jj=ii:N
        K(ii, jj) = exp(-sum((X(ii,:)-X(jj,:)).^2)/sigma^2);
        K(jj, ii) = K(ii, jj);
    end
end
onevec = ones(N,1);
Kc = K - onevec*onevec'*K/N - K*onevec*onevec'/N + (onevec*onevec'/N)*K*(onevec*onevec'/N); 

% Kernel LDA eigenvalue problem:
y1 = zeros(N,1); y1(1:N1) = 1;
y2 = zeros(N,1); y2((N1+1):end)=1;

MM = N1*N2/N*Kc*(y1/N1-y2/N2)*(y1/N1-y2/N2)'*Kc;
NN = Kc*(eye(N) - y1*y1'/N1 - y2*y2'/N2)*Kc;
[V,D] = eig(MM,NN);
D = diag(D);
[D, idx] = sort(D, 'ascend');
V = V(:, idx);
alpha = V(:,end);

% Generalization
[xx, yy] = meshgrid(linspace(-5,5,Ngrid), linspace(-5,5,Ngrid));
xtest = [xx(:) yy(:)];
ytest_KLDA = zeros(size(xtest,1),1);
for ii=1:length(xtest)
    ktest = exp(-sum((xtest(ii,:)-X).^2,2)/sigma);
    ktestc = ktest - Kc*ones(N,1)/N - ones(N,1)*ones(N,1)'*ktest/N + (ones(N,1)*ones(N,1)'/N)*Kc*(ones(N,1)/N);
    ytest_KLDA(ii) = 2*(alpha'*ktestc>0)-1;
end

% Visualizing KLDA result:
figure(3); clf; hold on
for ii=1:length(xtest)
    if ytest_KLDA(ii)==1
        plot(xtest(ii,1), xtest(ii,2), '.', 'color', [0.8 0.8 1.0], 'markersize', 15);
    else
        plot(xtest(ii,1), xtest(ii,2), '.', 'color', [1.0 0.8 0.8], 'markersize', 15);
    end
end

plot(X(1:200,1), X(1:200,2), 'bx');
plot(X(201:400,1), X(201:400,2), 'ro');
xlabel('$x$', 'Interpreter','latex', 'FontSize',12); 
ylabel('$y$', 'Interpreter','latex', 'FontSize',12);
box on;
set(gca, 'plotboxaspectratio', [1 1 1], 'xtick', [], 'ytick', []);
set(gcf, 'color', 'w');

toc

%% Support Vector Machine (SVM):
disp('Computing Support Vector Machine ...'); tic;
% Gaussian kernel:
model_gaussian = fitcsvm(X,Y,'kernelfunction','gaussian');

% Generalization
[xx, yy] = meshgrid(linspace(-5,5,Ngrid), linspace(-5,5,Ngrid));
xtest = [xx(:) yy(:)];
ytest_SVM = zeros(size(xtest,1),1);
for ii=1:length(xtest)
    ytest_SVM(ii) = model_gaussian.predict([xtest(ii,1) xtest(ii,2)]);
end

% Visualizing SVM result:
figure(4); clf; hold on
for ii=1:length(xtest)
    if ytest_SVM(ii)==1
        plot(xtest(ii,1), xtest(ii,2), '.', 'color', [0.8 0.8 1.0], 'markersize', 15);
    else
        plot(xtest(ii,1), xtest(ii,2), '.', 'color', [1.0 0.8 0.8], 'markersize', 15);
    end
end

plot(X(1:200,1), X(1:200,2), 'bx');
plot(X(201:400,1), X(201:400,2), 'ro')
xlabel('$x$', 'Interpreter','latex', 'FontSize',12); 
ylabel('$y$', 'Interpreter','latex', 'FontSize',12);;
box on;
set(gca, 'plotboxaspectratio', [1 1 1], 'xtick', [], 'ytick', []);
set(gcf, 'color', 'w');

toc