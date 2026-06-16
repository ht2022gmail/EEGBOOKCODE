clear all; close all; clc;

params.N = 50;
params.K = 0.1*params.N;
params.freq = 10;    % frequency
params.omega = 2*pi*params.freq;    % angular frequency
% params.omega = 2*pi*params.freq*(1+0.1*randn(params.N,1));
params.A = ones(params.N)-eye(params.N);

theta0 = 2*pi*rand(params.N,1)-pi;
T = 2;
dt = 0.001;
t = 0:dt:T;

%% original formulation:
[t, theta] = ode45(@(t,x) kuramoto_original(t,x,params), t, theta0);
theta = angle(exp(1i*theta));
order_parameter = abs(sum(exp(1i*theta),2))/params.N;


%% analytical formulation:
x0 = exp(1i*theta0);
x = zeros(length(t), params.N);
x(1,:) = x0';

A = expm(dt*(1i*params.omega*eye(params.N) + params.K/params.N*params.A));
for nn=1:length(t)-1
    tmpx = A*x(nn,:)'; % dynamics (described in the paper)
    tmpx = exp(1i*angle(tmpx)); % norm normalization (not described in the paper)
    x(nn+1,:) = tmpx';
end
theta1 = angle(x);
order_parameter1 = abs(sum(exp(1i*theta1),2))/params.N;

figure(2); clf; 
tiledlayout(1,2);
nexttile;
plot(t, order_parameter1, 'k');
axis([0 T 0 1]);
nexttile;
imagesc(1:params.N, 0:dt:T,  theta1);

%% Visualization:
figure(1); clf; 
subplot(131); hold on;
plot(t, order_parameter, 'k', 'LineWidth', 2);
plot(t, order_parameter1, 'r--', 'LineWidth', 2);
axis([0 T 0 1]);
xlabel('Time (s)'); ylabel('Order parameter');
set(gca,'PlotBoxAspectRatio', [1 1 1]);

subplot(132); 
imagesc(1:params.N, 0:dt:T,  theta);
xlabel('Unit'); ylabel('Time');
title('Original Kuramoto model');
set(gca,'PlotBoxAspectRatio', [1 1 1]);

subplot(133); 
imagesc(1:params.N, 0:dt:T,  theta);
xlabel('Unit'); ylabel('Time');
title('Complex Kuramoto model');
set(gca,'PlotBoxAspectRatio', [1 1 1]);

set(gcf, 'Color', 'w')


%% Functions:

function dxdt = kuramoto_original(t, x, params)
    omega = params.omega;
    N = params.N;
    K = params.K;
    A = params.A;
    x_diff = repmat(x',N,1) - x*ones(1,N);
    dxdt = omega + K/N*sum(A.*sin(x_diff),2);
end

% function dxdt = kuramoto_complex(t, x, params)
%     omega = params.omega;
%     N = params.N;
%     K = params.K;
%     A = ones(N)-eye(N);
%     dxdt = (1i*omega*eye(N) + K/N*A)*x;
% end