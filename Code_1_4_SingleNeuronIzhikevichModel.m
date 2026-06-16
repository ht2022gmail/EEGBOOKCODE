clear all; close all; clc
addpath('./functions');

%% Single neuron (Izhikevich model):
% Resonator (RZ) neuron
a = 0.1;
b = 0.3;
c = -65;
d = 2;
v0 = -65;    % Initial values of v
u0 = b*v0;   % Initial values of u
firings=[];  % spike timings

T = 1000000;  % in ms
I = zeros(T, 1);
I(1:100) = -10; I(100:110) = 10;
v = zeros(T, 1);
u = zeros(T, 1);
v(1) = v0; u(1) = u0;

for t=1:T-1            % simulation of 1000 ms
    if v(t)>=30
        firings =[firings; t];
        v(t-1) = 50;
        v(t) = c;
        u(t) = u(t) + d;
    end
  vtmp = v(t) + 0.5*(0.04*v(t)^2+5*v(t)+140-u(t)+I(t)); % step 0.5 ms
  v(t+1) = vtmp +0.5*(0.04*vtmp^2+5*vtmp+140-u(t)+I(t)); % for numerical
  u(t+1) = u(t) + a*(b*v(t+1)-u(t));                 % stability
end


figure(1); clf; hold on;
fill([0.001 0.1 0.1 0.001 0.001], [-100 -100 60 60 -100], [0.9 0.9 0.9], 'edgecolor', 'none')
fill([0.10 0.11 0.11 0.10 0.10], [-100 -100 60 60 -100], [0.7 0.7 0.7], 'edgecolor', 'none')
plot([0:T-1]/1000, v, 'k', 'linewidth', 2);
plot((firings-2)/1000, 50, 'r.', 'markersize', 10);
xlim([0 0.3]); ylim([-100 60]);
xlabel('Time (s)', 'interpreter', 'latex');
ylabel('Potential (mV)', 'interpreter', 'latex');
set(gca, 'plotboxaspectratio', [2 1 1]);
set(gcf, 'color', 'w');

v = v(1:T>100); % v = v - mean(v);
Nsamples = length(v);
fsampling = 1000;
win = fsampling * 2; 
overlap = 0.5*win;  % 50% overlap 
nfft = win;         % # samples for FFT
[pv_welch,f_welch] = ...
    pwelch(v, win, overlap, nfft, fsampling);

figure(2); clf; 
plot(f_welch, 10*log10(pv_welch), 'k', 'linewidth', 2);
xlim([10 80])
xlabel('Frequency (Hz)', 'interpreter', 'latex');
ylabel('$10 \log_{10}$(Power) (dB)', 'interpreter', 'latex');
set(gca,  'plotboxaspectratio', [1.5 1 1]);
set(gcf, 'color', 'w');

%% Using FFT (this result is too noisy to show):
% freq = [0:Nsamples-1]/Nsamples*fsampling;
% fv = fft(v);
% figure(3); clf; 
% plot(freq, 10*log10(abs(fv).^2), 'k', 'linewidth', 2);
% xlim([10 80])
% xlabel('Frequency (Hz)', 'interpreter', 'latex');
% ylabel('$10 \log_{10}$(Power) (dB)', 'interpreter', 'latex');
% set(gca, 'xaxislocation', 'origin', 'yaxislocation', 'origin', 'plotboxaspectratio', [2 1 1]);
% set(gcf, 'color', 'w');
