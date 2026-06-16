clear all; close all; clc;
load('./data/sampledataEEG.mat');
EEG  = ALLEEG(1); % loading continuous data

%% Time-domain and frequency-domain convolution
% Select an electrode for example:
elec_name = 'Oz';
elec = find(strcmp({EEG.chanlocs.labels}, elec_name));
Xchan = squeeze(EEG.data(elec,:)); 

% Define Molret wavelet:
time = -1:1/EEG.srate:1;
freq = 10;
s = 5/(2*pi*freq);
molret = exp((-time.^2)/(2*s^2))/30;
 
n_conv = length(Xchan) + length(molret) - 1;
half = ceil(length(molret)/2);
 
% Time-domain convolution:
conv_time = conv(Xchan, molret, 'same');

% frequency-based convolution (convolution theorem):
fft_Xchan = fft(Xchan, n_conv);
fft_molret = fft(molret, n_conv);
ift   = ifft(fft_Xchan.*fft_molret, n_conv);
conv_FFT = real(ift(half:end-half+1)); % cut the initial and last segments.

%%
figure(1); clf; hold on;
plot(EEG.times/1000, Xchan, 'k');
plot(EEG.times/1000, conv_time, 'r', 'linewidth', 2);
xlim([10 20]);
xlabel('Time (s)');
ylabel('EEG (\mu V)');
title(['Electrode ' elec_name ': time-based convolution']);
legend('original EEG', 'wavelet-convolved EEG');
set(gca, 'plotboxaspectratio', [3 1 1]);
set(gcf, 'color', 'w');

figure(2); clf; hold on;
plot(EEG.times/1000, Xchan, 'k');
plot(EEG.times/1000, conv_FFT, 'r', 'linewidth', 2);
xlim([10 20]);
xlabel('Time (s)');
ylabel('EEG (\mu V)');
title(['Electrode ' elec_name ': FFT-based convolution']);
legend('original EEG', 'wavelet-convolved EEG');
set(gca, 'plotboxaspectratio', [3 1 1]);
set(gcf, 'color', 'w');
