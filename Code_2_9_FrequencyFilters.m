clear all; close all; clc;
load('./data/sampledataEEG.mat');
EEG  = ALLEEG(1); % loading continuous data

%% Frequency filters:
X = double(EEG.data);
[Nelectrodes, Ntime, Ntrials] = size(X);

filterLength = 200;

% Select an electrode for example:
elec_name = 'Oz';
elec = find(strcmp({EEG.chanlocs.labels}, elec_name));
Xchan = squeeze(X(elec,:,:)); % time * trials

% Bandpass filter:
centerFreq_BP = 10; % in Hz
filterFrequencySpread_BP  = 2; % Hz +/- the center frequency
transitionWidth_BP= 0.2;
ffrequencies_BP   = ...
    [ 0 (1-transitionWidth_BP)*(centerFreq_BP-filterFrequencySpread_BP) (centerFreq_BP-filterFrequencySpread_BP) (centerFreq_BP+filterFrequencySpread_BP) (1+transitionWidth_BP)*(centerFreq_BP+filterFrequencySpread_BP) EEG.srate/2 ]/(EEG.srate/2);
idealResponse_BP  = [ 0 0 1 1 0 0 ];
bandpassfilter  = firls(filterLength,ffrequencies_BP,idealResponse_BP);

% Highpass filter
cutoffFreq_HP = 20; % in Hz
filterFrequencySpread_HP  = 2; % Hz +/- the center frequency
transitionWidth_HP= 0.2;
% cutoffFreq_HP = 2; % in Hz
% filterFrequencySpread_HP  = 0.3; % Hz +/- the center frequency
% transitionWidth_HP= 0.;
ffrequencies_HP   = ...
    [ 1 (1-transitionWidth_HP)*(cutoffFreq_HP-filterFrequencySpread_HP) (cutoffFreq_HP-filterFrequencySpread_HP)  EEG.srate/2 ]/(EEG.srate/2);
idealResponse_HP  = [ 0 0 1 1 ];
highpassfilter  = firls(filterLength,ffrequencies_HP,idealResponse_HP);

% Lowpass filter
cutoffFreq_LP = 7; % in Hz
filterFrequencySpread_LP  = 1; % Hz +/- the center frequency
transitionWidth_LP= 0.2;
ffrequencies_LP   = ...
    [ 1 (1-transitionWidth_LP)*(cutoffFreq_LP-filterFrequencySpread_LP) (cutoffFreq_LP-filterFrequencySpread_LP)  EEG.srate/2 ]/(EEG.srate/2);
idealResponse_LP  = [ 1 1 0 0 ];
lowpassfilter  = firls(filterLength,ffrequencies_LP,idealResponse_LP);

% Applying the frequency filters to data:
Xchan_BP = filtfilt(bandpassfilter, 1, Xchan);
Xchan_HP = filtfilt(highpassfilter, 1, Xchan);
Xchan_LP = filtfilt(lowpassfilter, 1, Xchan);

% Checking the 
ff = [0:filterLength]/(filterLength+1)*EEG.srate;
fft_bandpassfilter  = abs(fft(bandpassfilter));
fft_bandpassfilter  = fft_bandpassfilter./max(fft_bandpassfilter); 
fft_highpassfilter  = abs(fft(highpassfilter));
fft_highpassfilter  = fft_highpassfilter./max(fft_highpassfilter); 
fft_lowpassfilter  = abs(fft(lowpassfilter));
fft_lowpassfilter  = fft_lowpassfilter./max(fft_lowpassfilter); 

%% Plot the results:
figure(1); clf;
subplot(3,2,1);
plot(EEG.times/1000, Xchan_BP, 'k')
xlabel('Time (s)'); ylabel('EEG (\mu V)');
title('Bandpass filtered')
xlim([100 105])
set(gca, 'plotboxaspectratio', [4 1 1]);

subplot(3,2,3)
plot(EEG.times/1000, Xchan_HP, 'k')
xlabel('Time (s)'); ylabel('EEG (\mu V)');
title('Highpass filtered')
xlim([100 105])
set(gca, 'plotboxaspectratio', [4 1 1]);

subplot(3,2,5)
plot(EEG.times/1000, Xchan_LP, 'k')
xlabel('Time (s)'); ylabel('EEG (\mu V)');
title('Lowpass filtered')
xlim([100 105])
set(gca, 'plotboxaspectratio', [4 1 1]);

subplot(3,2,2); hold on;
plot(ffrequencies_BP*EEG.srate/2, idealResponse_BP, 'k')
plot(ff, fft_bandpassfilter, 'r--')
xlabel('Frequency (Hz)'); % ylabel('EEG (\mu V)');
title('Bandpass filter response')
xlim([0 EEG.srate/2/2]); box on;
set(gca, 'plotboxaspectratio', [4 1 1]);

subplot(3,2,4); hold on;
plot(ffrequencies_HP*EEG.srate/2, idealResponse_HP, 'k')
plot(ff, fft_highpassfilter, 'r--')
xlabel('Frequency (Hz)'); % ylabel('EEG (\mu V)');
title('Highpass filter response')
xlim([0 EEG.srate/2/2]); box on;
set(gca, 'plotboxaspectratio', [4 1 1]);

subplot(3,2,6); hold on;
plot(ffrequencies_LP*EEG.srate/2, idealResponse_LP, 'k')
plot(ff, fft_lowpassfilter, 'r--')
xlabel('Frequency (Hz)'); % ylabel('EEG (\mu V)');
title('Lowpass filter response')
xlim([0 EEG.srate/2/2]); box on;
set(gca, 'plotboxaspectratio', [4 1 1]);

set(gcf, 'color', 'w');