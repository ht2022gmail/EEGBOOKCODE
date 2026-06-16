%% EEG Alpha Waves dataset
% https://zenodo.org/record/2348892
% FP1, FP2, FC5, FC6, FZ, T7, CZ, T8, P7, P3, PZ, P4, P8, O1, Oz, and O2.

clear all; close all; clc;

%% Loading data and defining some variables:
subnum = 0;
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
filter_low    = 3; % Hz
filter_high   = 30; % Hz
ffrequencies  = [ 0 filter_low*(1-transition_width) filter_low filter_high filter_high*(1+transition_width) Fnyquist ]/Fnyquist;
idealresponse = [ 0 0 1 1 0 0 ];
filterweights = firls(round(2*(srate/filter_low)),ffrequencies,idealresponse);
x      = filtfilt(filterweights, 1, x);
x = x/std(x);

figure(1); clf; 
subplot(211); hold on;
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
ylabel('EEG (normalized)')
title('Oz (bandpass filtered [3 Hz, 30 Hz])')
box on;
set(gcf, 'color', 'w');

%% Periodgram
win = srate * 1;    % 2-second window
overlap = 0.5*win;  % 50% overlap 
nfft = win;         % # samples for FFT

time_win = [0:1:time(end)];
time_win = time_win(1:end-1);

Px = [];

for nn=1:length(time_win)
    idx_start = dsearchn(time', time_win(nn));
    idx_end = dsearchn(time', time_win(nn)+3);

    xseg = x(idx_start:idx_end);
    [px_welch,f_welch, pxc_welch] = ...
        pwelch(xseg, win, overlap, nfft, srate, 'ConfidenceLevel',0.95);
    Px = [Px px_welch];
end
subplot(212); hold on;
imagesc(time_win+1.5,f_welch(1:31), 10*log10(Px(1:31,:)))
set(gca, 'clim', [-5 max(max(10*log10(Px(1:31,:))))]);

for onset = onsets_eyeclose'
     plot([time(onset) time(onset)], [0 30], 'r--', 'linewidth', 2);
end
for onset = onsets_eyeopen'
     plot([time(onset) time(onset)], [0 30], 'y--', 'linewidth', 2);
end
axis([time_win(1)+1.5 time_win(end)+1.5 3 20])
xlabel('Time (s)');
ylabel('Frequency (Hz)')
title(['Periodgram for subject #' num2str(subnum,'%02d')]);
set(gca, 'plotboxaspectratio', [3 1 1]);
colormap(jet);
% colorbar; 
set(gcf, 'color', 'w');