clear all; close all; clc;
load('./data/sampledataEEG.mat');
EEG  = ALLEEG(1); % loading continuous data
addpath('./functions/');

%% Frequency filters:
X = double(EEG.data);
[Nelectrodes, Ntime, Ntrials] = size(X);

filterLength = 200;

% Select an electrode for example:
elec_name = 'Oz';
elec = find(strcmp({EEG.chanlocs.labels}, elec_name));
x = squeeze(X(elec,:,:)); % time * trials

% Lowpass filter
cutoffFreq_LP = 50; % in Hz
filterFrequencySpread_LP  = 2; % Hz +/- the center frequency
transitionWidth_LP= 0.1;
ffrequencies_LP   = ...
    [ 1 (1-transitionWidth_LP)*(cutoffFreq_LP-filterFrequencySpread_LP) (cutoffFreq_LP-filterFrequencySpread_LP)  EEG.srate/2 ]/(EEG.srate/2);
idealResponse_LP  = [ 1 1 0 0 ];
lowpassfilter  = firls(filterLength,ffrequencies_LP,idealResponse_LP);

% Applying the frequency filters to data:
x_lowpass = filtfilt(lowpassfilter, 1, x);
% x_lowpass = x;
%% Fourier transform:
% naive FFT:
% f = [0:Ntime-1]/Ntime*EEG.srate;
% fx_lowpass = fft(x_lowpass);
% pow_lowpass = 1/(EEG.srate*Ntime)*abs(fx_lowpass).^2;
% figure(1); clf; 
% plot(f, 10*log10(pow_lowpass), 'k');
% xlim([0 35])

win = EEG.srate * 1;    % 3-second window
overlap = 0.5*win;  % 50% overlap 
nfft = win;         % # samples for FFT
[pow_lowpass, f, cl_lowpass] = pwelch(x_lowpass, win, overlap, nfft, EEG.srate, 'ConfidenceLevel', 0.95);
figure(1); clf; hold on
plot(f, 10*log10(pow_lowpass), 'k', 'LineWidth', 2);
plot(f, 10*log10(pow_lowpass+cl_lowpass(:,2)), 'k--');
plot(f, 10*log10(pow_lowpass-cl_lowpass(:,2)), 'k--');
xlim([0 40]); ylim([-15 30]);
xlabel('Frequency (Hz)');
ylabel('log_{10} (EEG power)');
set(gca, 'PlotBoxAspectRatio', [1.5 1 1]);
set(gcf, 'color', 'w');

%% Downsampling 
factor_downsample = 8;
srate_downsample = EEG.srate/factor_downsample;
x_lowpass_downsample = downsample(x, factor_downsample);
times_downsample = downsample(EEG.times, factor_downsample);

% Interpolation:
dt = 1000/srate_downsample;
x_interp = Whittaker_Shannon_interp(EEG.times, dt, ...
    times_downsample, x_lowpass_downsample);
% fx_interp = fft(x_interp);
% pow_interp = 1/(EEG.srate*Ntime)*abs(fx_interp).^2;
[pow_interp, f, cl_interp] = pwelch(x_interp, win, overlap, nfft, EEG.srate, 'ConfidenceLevel', 0.95);
figure(2); clf; hold on;
plot(f, 10*log10(pow_interp), 'k', 'LineWidth', 2);
plot(f, 10*log10(pow_interp+cl_interp(:,2)), 'k--');
plot(f, 10*log10(pow_interp-cl_interp(:,2)), 'k--');
xlim([0 40]); ylim([-15 30]);
xlabel('Frequency (Hz)');
ylabel('log_{10} (EEG power)');
set(gca, 'PlotBoxAspectRatio', [1.5 1 1]);
set(gcf, 'color', 'w');

%% Visualizing time series:

figure(3); clf; hold on;
plot(EEG.times/1000, x_lowpass, 'k-');
plot(EEG.times/1000, x_lowpass, 'k.', 'Markersize', 15);
xlim([10 11]);
xlabel('Time (ms)');
ylabel('EEG (\mu V)');
set(gca, 'PlotBoxAspectRatio', [1.5 1 1], 'xtick', 10:0.2:11);
set(gcf, 'color', 'w');

figure(4); clf; hold on;
plot(EEG.times/1000, x_interp, 'r-');
plot(EEG.times/1000, x_lowpass, '--', 'color', [0.3 0.3 0.3]);
plot(times_downsample/1000, x_lowpass_downsample, 'r.', 'MarkerSize', 15);
xlim([10 11]);
xlabel('Time (ms)');
ylabel('EEG (\mu V)');
set(gca, 'PlotBoxAspectRatio', [1.5 1 1], 'xtick', 10:0.2:11);
set(gcf, 'color', 'w');

