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
x = eeg_oz;

%%
win = srate;

f = 10.^(linspace(log10(2), log10(60), 30));

[xspec, fspec, tspec] = spectrogram(x, hann(win), win/2, f, srate);
figure(1); clf; hold on;
imagesc(tspec, f/2, 10*log10(abs(xspec).^2))
% set(gca, 'clim', [30 70], 'ydir', 'normal', 'plotboxaspectratio', [3 1 1], ...
%     'xtick', [0:20:120], 'ytick', f, ...
%     'fontsize', 12);
set(gca, 'clim', [30 70], 'ydir', 'normal', 'plotboxaspectratio', [3 1 1], 'fontsize', 14)
axis tight
xlabel('Time (s)', 'fontsize', 16);
ylabel('Frequency (Hz)', 'fontsize', 16);
colormap(gray)
set(gcf, 'color', 'w');