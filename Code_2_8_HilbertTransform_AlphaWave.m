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

%% Bandpass filtering:
transition_width = 0.2; % percent
filter_low    = 8; % Hz
filter_high   = 12; % Hz
ffrequencies  = [ 0 filter_low*(1-transition_width) filter_low filter_high filter_high*(1+transition_width) Fnyquist ]/Fnyquist;
idealresponse = [ 0 0 1 1 0 0 ];
filterweights = firls(round(2*(srate/filter_low)),ffrequencies,idealresponse);
x      = filtfilt(filterweights, 1, x);
x = x/std(x);

figure(1); clf; hold on;
plot(time, x, 'k');
onsets_eyeclose = find(SIGNAL(:,18)==1);
onsets_eyeopen = find(SIGNAL(:,19)==1);
for onset = onsets_eyeclose'
     plot([time(onset) time(onset)], [-5 5], 'r--', 'linewidth', 2);
end
for onset = onsets_eyeopen'
     plot([time(onset) time(onset)], [-5 5], 'b--', 'linewidth', 2);
end
set(gca, 'plotboxaspectratio', [3 1 1], 'xlim', [time(1) time(end)]);
xlabel('Time (s)');
title('Oz (bandpass filtered [3 Hz, 30 Hz])')
box on;
set(gcf, 'color', 'w');

%% 
z = hilbert(x);
plot(time, abs(z), 'c')
plot(time, -abs(z), 'c')

figure(2); clf; hold on;
plot3(time, real(z), imag(z), 'k');
plot3(time, real(z), -5*ones(1,Ntime), 'color', [0.8 0.8 0.8]);
plot3(time, 5*ones(1,Ntime), imag(z), 'color', [0.8 0.8 0.8]);
set(gca, 'plotboxaspectratio', [3 1 1]);
xlabel('Time (s)');
ylabel('real(EEG) (\mu V)')
zlabel('imag(EEG) (\mu V)')
view([-20 35])
xlim([time(1) time(end)])
set(gcf, 'color', 'w');

