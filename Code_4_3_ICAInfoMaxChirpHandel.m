clear all; close all; clc;
addpath('./functions');

%% Create artificial data:
N = 2;  % # source signals
T = 1.3e4;  % # signal length
Fs=10000;
t = 0:1/Fs:(T-1)/Fs;
load chirp; s1=y(1:T); s1=s1/std(s1);
load handel; s2=y(1:T); s2=s2/std(s2);

% Combine sources into vector variable s.
s=[s1, s2]';

% Generate observed data by random mixing:
A = randn(N,N);
x = A*s;

%% Independent component analyses (JADE, fastICA, InfoMax):
% JADE
Wjade = jadeR(x);
yjade = Wjade*x;

% fastICA
[yfica, Afica, Wfica] = fastica(x);

% infoMax ICA
[yimax, Wimax] = infomaxICA(x);

% Amari Index
disp(['Amari Index (Infomax): ' num2str(AmariIndex(Wimax*A))]);
disp(['Amari Index (FastICA): ' num2str(AmariIndex(Wfica*A))]);
disp(['Amari Index (JADE):    ' num2str(AmariIndex(Wjade*A))]);

%% Visualizing the results:
figure(1); clf;
subplot(2,5,1); plot(t, s(1,:), 'k');
set(gca, 'plotboxaspectratio', [1 1 1])
xlabel('Time (s)');
title('Original $s_1$', 'interpreter', 'latex')
subplot(2,5,6); plot(t, s(2,:), 'k');
set(gca, 'plotboxaspectratio', [1 1 1])
xlabel('Time (s)');
title('Original $s_2$', 'interpreter', 'latex')

subplot(2,5,2); plot(t, x(1,:), 'k');
set(gca, 'plotboxaspectratio', [1 1 1])
xlabel('Time (s)');
title('Mixed $x_1$', 'interpreter', 'latex')
subplot(2,5,7); plot(t, x(2,:), 'k');
set(gca, 'plotboxaspectratio', [1 1 1])
xlabel('Time (s)');
title('Mixed $x_2$', 'interpreter', 'latex')

subplot(2,5,3); plot(t, yimax(1,:), 'k');
set(gca, 'plotboxaspectratio', [1 1 1])
xlabel('Time (s)');
title('InfoMax ICA $y_1$', 'interpreter', 'latex')
subplot(2,5,8); plot(t, yimax(2,:), 'k');
set(gca, 'plotboxaspectratio', [1 1 1])
xlabel('Time (s)');
title('InfoMax ICA $y_2$', 'interpreter', 'latex')

subplot(2,5,4); plot(t, yfica(1,:), 'k');
set(gca, 'plotboxaspectratio', [1 1 1])
xlabel('Time (s)');
title('FastICA $y_1$', 'interpreter', 'latex')
subplot(2,5,9); plot(t, yfica(2,:), 'k');
set(gca, 'plotboxaspectratio', [1 1 1])
xlabel('Time (s)');
title('FastICA $y_2$', 'interpreter', 'latex')

subplot(2,5,5); plot(yjade(1,:), 'k');
set(gca, 'plotboxaspectratio', [1 1 1])
title('JADE $y_1$', 'interpreter', 'latex')
subplot(2,5,10); plot(yjade(2,:), 'k');
set(gca, 'plotboxaspectratio', [1 1 1])
title('JADE $y_2$', 'interpreter', 'latex')

set(gcf, 'color', 'w');


%% Sonifying the results:
disp('Playing s1 (chirp) ...');
soundsc(s1);
pause;

disp('Playing s2 (Handel) ...');
soundsc(s2);
pause;

disp('Playing x1 ...');
soundsc(x(1,:));
pause; 

disp('Playing x2 ...');
soundsc(x(2,:));
pause; 

disp('Playing y1 (InfoMax) ...');
soundsc(yimax(1,:));
pause; 

disp('Playing y2 (InfoMax) ...');
soundsc(yimax(2,:));

