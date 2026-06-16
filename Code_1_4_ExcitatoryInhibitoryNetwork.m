clear all; close all; clc
addpath('./functions');

%% Setting parameters for the excitatory-inhibitory network:
% See Chapter 7.5 of Theoretical Neuroscience for details.
% We change the WEE parameter for Hopf bifurcations.
% params.WEE = 1.25;  % this is stable.
params.WEE = 1.4;   % this is oscillatory.
params.WIE = 1;
params.WEI = -1;
params.WII = 0;
params.IE = 10;
params.II = -10;
params.tauE = 10/1000;
params.tauI = 30/1000;

r0 = [50; 30];
[t, r] = ode45(@(t,r) excitatoryinhinitorynetwork(t,r,params), [0:(1/1000):10], r0);

%% Visualizing the results:
figure(1); clf; hold on;
plot(t, r(:,1), 'k', 'linewidth', 1.5);
plot(t, r(:,2), 'r', 'linewidth', 1.5);
xlim([0 2])
xlabel('Time (s)', 'interpreter', 'latex');
ylabel('Activity (Hz)', 'interpreter', 'latex');
leg = legend('Excitatory activity $r_{\rm{E}}$', 'Inhibitory activity $r_{\rm{I}}$');
set(leg, 'interpreter', 'latex');
set(gca, 'xaxislocation', 'origin', 'yaxislocation', 'origin', 'plotboxaspectratio', [2 1 1]);
set(gcf, 'color', 'w');

figure(2); clf; hold on;
plot(r(:,1), r(:,2), 'k', 'linewidth', 1.5);
axis([-10 75 -10 60]);
xlabel('Excitatory activity $r_{\rm{E}}$ (Hz)', 'interpreter', 'latex');
ylabel('Inhibitory activity $r_{\rm{I}}$ (Hz)', 'interpreter', 'latex');
set(gca, 'xaxislocation', 'origin', 'yaxislocation', 'origin', 'plotboxaspectratio', [1 1 1]);
set(gcf, 'color', 'w');

% Native FFT:
r = bsxfun(@minus, r, mean(r));
fr = fft(r);
freq = [0:(length(t)-1)]/length(t)*1000;

figure(3); clf; hold on;
plot(freq, 10*log10(abs(fr(:,1)).^2), 'k', 'linewidth', 1.5)
xlim([1 30])
xlabel('Frequency (Hz)', 'interpreter', 'latex');
ylabel('Power (dB)', 'interpreter', 'latex');
set(gca, 'xaxislocation', 'origin', 'yaxislocation', 'origin', 'plotboxaspectratio', [2 1 1]);
set(gcf, 'color', 'w');

% Welch's method:
% r = bsxfun(@minus, r, mean(r));
% Nsamples = length(r);
% fsampling = 1000;
% win = fsampling * 2; 
% overlap = 0.5*win;  % 50% overlap 
% nfft = win;         % # samples for FFT
% [pr_welch,f_welch] = ...
%     pwelch(r(:,1), win, overlap, nfft, fsampling);
% 
% figure(3); clf; 
% plot(f_welch, pr_welch, 'k', 'linewidth', 2);
% xlim([1 20])
% xlabel('Frequency (Hz)', 'interpreter', 'latex');
% ylabel('$10 \log_{10}$(Power) (dB)', 'interpreter', 'latex');
% set(gca, 'xaxislocation', 'origin', 'yaxislocation', 'origin', 'plotboxaspectratio', [1.5 1 1]);
% set(gcf, 'color', 'w');

%%
% WEE = params.WEE;
WIE = params.WIE;
WEI = params.WEI;
WII = params.WII;
IE = params.IE;
II = params.II;
tauE = params.tauE;
tauI = params.tauI;

WEE = [1.15:0.05:1.45];

lambda1 = ((WEE-1)/tauE + (WII-1)/tauI + sqrt(((WEE-1)/tauE-(WII-1)/tauI).^2+4*WEI*WIE/tauE/tauI))/2;
relambda1 = real(lambda1);
imlambda1 = imag(lambda1);

figure(4); clf;
subplot(211);
plot(WEE, relambda1, 'k'); 
xlim([WEE(1) WEE(end)]);
set(gca, 'xaxislocation', 'origin', 'yaxislocation', 'origin', 'plotboxaspectratio', [2 1 1]);

subplot(212);
plot(WEE, imlambda1/(2*pi), 'k'); 
xlim([WEE(1) WEE(end)]);
set(gca, 'xaxislocation', 'origin', 'yaxislocation', 'origin', 'plotboxaspectratio', [2 1 1]);

set(gcf, 'color', 'w');
