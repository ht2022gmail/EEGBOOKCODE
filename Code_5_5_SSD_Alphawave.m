%% EEG Alpha Waves dataset
% https://zenodo.org/record/2348892
% FP1, FP2, FC5, FC6, FZ, T7, CZ, T8, P7, P3, PZ, P4, P8, O1, Oz, and O2.

%% Loading data and defining some variables:
clear all; close all; clc;
subnum = 0;
% load subject_00.mat;
load(['./data/subject_' num2str(subnum, '%02d') '.mat']);
load('./data/channels_alphawavedataset.mat');
addpath('./functions/');
srate = 512;
Fnyquist = srate/2;
X = SIGNAL(:,2:17)'; 
[Nelectrodes, Ntime] = size(X);
Ntime = size(SIGNAL,1);
time = [0:(Ntime-1)]/srate;
% Bandpass filtering:
filterLength = 200;
cutoffFreq_HP = 2; % in Hz
filterFrequencySpread_HP  = 0.3; % Hz +/- the center frequency
transitionWidth_HP= 0.;
ffrequencies_HP   = ...
    [ 1 (1-transitionWidth_HP)*(cutoffFreq_HP-filterFrequencySpread_HP) (cutoffFreq_HP-filterFrequencySpread_HP)  srate/2 ]/(srate/2);
idealResponse_HP  = [ 0 0 1 1 ];
highpassfilter  = firls(filterLength,ffrequencies_HP,idealResponse_HP);
for nn=1:Nelectrodes
    X(nn,:) = filtfilt(highpassfilter, 1, X(nn,:));
end

onsets_eyeclose = find(SIGNAL(:,18)==1);
onsets_eyeopen = find(SIGNAL(:,19)==1);

X = double(X)'; % data must be in TIME x CHANNEL format.

%% SSD
freq_alpha = [8 12; 2 20; 7 13];    % Entire alpha band
freq_low= [8 10; 2 20; 7 11];    % Lower alpha
freq_high = [10 12; 2 20; 9 13];     % Higher alpha
sampling_freq = 512;
filter_order = 2;
[W_alpha, A_alpha, lambda_alpha, C_s_alpha, X_alpha_ssd] = ssd(X, freq_alpha, sampling_freq, filter_order, []);
[W_low, A_low, lambda_low, C_s_low, X_low_ssd] = ssd(X, freq_low, sampling_freq, filter_order, []);
[W_high, A_high, lambda_high, C_s_high, X_high_ssd] = ssd(X, freq_high, sampling_freq, filter_order, []);

% 
x_alpha = X_alpha_ssd(:,1)./std(X_alpha_ssd(:,1));
x_low = X_low_ssd(:,1)./std(X_low_ssd(:,1));
x_high = X_high_ssd(:,1)./std(X_high_ssd(:,1));

figure(1); clf; hold on;
for kk=1:length(onsets_eyeopen)
    t1 = time(onsets_eyeclose(kk));
    t2 = time(onsets_eyeopen(kk));
    ff = fill([t1 t2 t2 t1 t1], 4.9*[-1 -1 1 1 -1], [0.9 0.9 0.9]);
    set(ff, 'edgecolor', [1 1 1])
end
plot(time, x_alpha, 'k');
plot(time, abs(hilbert(x_alpha)), 'r--');
plot(time, -abs(hilbert(x_alpha)), 'r--');
set(gca, 'PlotBoxAspectRatio', [4 1 1], ...
    'XTick', [0:20:120], 'YTick', [-5:2.5:5]);
axis([time(1) time(end) -5 5])
xlabel('Time (s)'); ylabel('SSD component (a.u.)')
set(gcf, 'Color', 'w')

figure(2); clf; hold on;
for kk=1:length(onsets_eyeopen)
    t1 = time(onsets_eyeclose(kk));
    t2 = time(onsets_eyeopen(kk));
    ff = fill([t1 t2 t2 t1 t1], 4.9*[-1 -1 1 1 -1], [0.9 0.9 0.9]);
    set(ff, 'edgecolor', [1 1 1])
end
plot(time, x_low, 'k');
plot(time, abs(hilbert(x_low)), 'r--');
plot(time, -abs(hilbert(x_low)), 'r--');
set(gca, 'PlotBoxAspectRatio', [4 1 1], ...
    'XTick', [0:20:120], 'YTick', [-5:2.5:5]);
axis([time(1) time(end) -5 5])
xlabel('Time (s)');  ylabel('SSD component (a.u.)')
set(gcf, 'Color', 'w')

figure(3); clf; hold on;
for kk=1:length(onsets_eyeopen)
    t1 = time(onsets_eyeclose(kk));
    t2 = time(onsets_eyeopen(kk));
    ff = fill([t1 t2 t2 t1 t1], 4.9*[-1 -1 1 1 -1], [0.9 0.9 0.9]);
    set(ff, 'edgecolor', [1 1 1])
end
plot(time, x_high, 'k');
plot(time, abs(hilbert(x_high)), 'r--');
plot(time, -abs(hilbert(x_high)), 'r--');
set(gca, 'PlotBoxAspectRatio', [4 1 1], ...
    'XTick', [0:20:120], 'YTick', [-5:2.5:5]);
axis([time(1) time(end) -5 5])
xlabel('Time (s)');  ylabel('SSD component (a.u.)')
set(gcf, 'Color', 'w')

figure(4); clf;
topoplot(A_alpha(:,1), chanlocs_alphawavedataset, 'style', 'map', 'electrodes', 'ptslabels');
set(gcf, 'Color', 'w')

figure(5); clf;
topoplot(A_low(:,1), chanlocs_alphawavedataset, 'style', 'map', 'electrodes', 'ptslabels');
set(gcf, 'Color', 'w')

figure(6); clf;
topoplot(A_high(:,1), chanlocs_alphawavedataset, 'style', 'map', 'electrodes', 'ptslabels');
set(gcf, 'Color', 'w')

%%
figure(7); clf;
topoplot(W_low(:,1), chanlocs_alphawavedataset, 'style', 'map', 'electrodes', 'ptslabels');
set(gcf, 'Color', 'w')

figure(8); clf;
topoplot(W_high(:,1), chanlocs_alphawavedataset, 'style', 'map', 'electrodes', 'ptslabels');
set(gcf, 'Color', 'w')



%%
% figure(7); clf; hold on;
% % plot(time, abs(hilbert(x_alpha)), 'k');
% plot(time, abs(hilbert(x_low)), 'r');
% plot(time, -abs(hilbert(x_low)), 'r');
% plot(time, abs(hilbert(x_high)), 'b');
% plot(time, -abs(hilbert(x_high)), 'b');
% set(gca, 'PlotBoxAspectRatio', [4 1 1], ...
%     'XTick', [0:20:120], 'YTick', [-5:2.5:5]);
% axis([time(1) time(end) -5 5])

%% Welch method:
win = sampling_freq * 3;    % 3-second window
overlap = 0.5*win; % 50% overlap 
nfft = win;         % # samples for FFT
[px_alpha_welch, f_welch, CI_alpha_welch] = ...
    pwelch(X_alpha_ssd(:,1), win, overlap, nfft, sampling_freq, 'ConfidenceLevel', 0.95);
[px_low_welch, f_welch, CI_low_welch] = ...
    pwelch(X_low_ssd(:,1), win, overlap, nfft, sampling_freq, 'ConfidenceLevel', 0.95);
[px_high_welch, f_welch, CI_high_welch] = ...
    pwelch(X_high_ssd(:,1), win, overlap, nfft, sampling_freq, 'ConfidenceLevel', 0.95);
