%% EEG Alpha Waves dataset
% https://zenodo.org/record/2348892
% FP1, FP2, FC5, FC6, FZ, T7, CZ, T8, P7, P3, PZ, P4, P8, O1, Oz, and O2.

%% Loading data and defining some variables:
clear all; close all; clc;
subnum = 0;
% load subject_00.mat;
load(['./data/subject_' num2str(subnum, '%02d') '.mat']);
srate = 512;
Fnyquist = srate/2;
Ntime = size(SIGNAL,1);
time = [0:(Ntime-1)]/srate;
eeg_oz = SIGNAL(:,16); 
eeg_oz = detrend(eeg_oz);
x = eeg_oz';

% create wavelet
wtime = -1:1/srate:1;
freq = 10; % frequency of sine wave in Hz
s = 4.5/(2*pi*freq); 
wavelet = exp(-wtime.^2./(2*s^2)).*exp(1i*2*pi*freq.*wtime);
% half of the wavelet size, useful for chopping off edges after convolution.
half = ceil(length(wavelet)/2);
 
% convolve with data
% compute Gaussian
n_conv = length(wavelet) + length(x) - 1;
 
fft_wavelet = fft(wavelet,n_conv);
fft_e = fft(x,n_conv);
ift   = ifft(fft_e.*fft_wavelet,n_conv)*sqrt(s)/10; % sqrt... is an empirical scaling factor that works here
wavelet_conv_data = real(ift(half:end-half+1));

figure(1); clf; 
subplot(311); hold on;
plot(time, wavelet_conv_data, 'k');
onsets_eyeclose = find(SIGNAL(:,18)==1);
onsets_eyeopen = find(SIGNAL(:,19)==1);
for onset = onsets_eyeclose'
     plot([time(onset) time(onset)], [-50 50], 'r--', 'linewidth', 2);
end
for onset = onsets_eyeopen'
     plot([time(onset) time(onset)], [-50 50], 'b--', 'linewidth', 2);
end
rectangle('position', [10 -50 3 100], 'edgecolor', 'r');
rectangle('position', [20 -50 3 100], 'edgecolor', 'b');
set(gca, 'plotboxaspectratio', [3 1 1], 'xlim', [time(1) time(end)], 'ylim', [-50 50]);
xlabel('Time (s)');
ylabel('EEG (normalized)')
title('Oz (wavelet convolved), Eyes open');
box on;
set(gcf, 'color', 'w');

subplot(312); hold on;
plot(time, wavelet_conv_data, 'k');
set(gca, 'plotboxaspectratio', [3 1 1], 'xlim', [10 13], 'ylim', [-50 50]);
xlabel('Time (s)');
ylabel('EEG (normalized)')
title('Oz (wavelet convolved)');
box on;
set(gcf, 'color', 'w');

subplot(313); hold on;
plot(time, wavelet_conv_data, 'k');
set(gca, 'plotboxaspectratio', [3 1 1], 'xlim', [20 23], 'ylim', [-50 50]);
xlabel('Time (s)');
ylabel('EEG (normalized)')
title('Oz (wavelet convolved), Eyes closed');
box on;
set(gcf, 'color', 'w');