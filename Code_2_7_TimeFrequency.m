%% EEG Alpha Waves dataset
% https://zenodo.org/record/2348892
% FP1, FP2, FC5, FC6, FZ, T7, CZ, T8, P7, P3, PZ, P4, P8, O1, Oz, and O2.

%% Loading data (alpha-wave) and defining some variables:
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

%% Bandpass filtering:
transition_width = 0.2; % percent
filter_low    = 2; % Hz
filter_high   = 40; % Hz
ffrequencies  = [ 0 filter_low*(1-transition_width) filter_low filter_high filter_high*(1+transition_width) Fnyquist ]/Fnyquist;
idealresponse = [ 0 0 1 1 0 0 ];
filterweights = firls(round(2*(srate/filter_low)),ffrequencies,idealresponse);
x      = filtfilt(filterweights, 1, x);
x = x/std(x);

%%
min_freq =  2;
max_freq = 30;
num_frex = 30;

% define wavelet parameters
wave_time = -1:1/srate:1;
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
s = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);
% s = 3./(2*pi*frex); % this line is for figure 13.14
% s = 10./(2*pi*frex); % this line is for figure 13.14

% definte convolution parameters
n_wavelet = length(wave_time);
n_data = length(x);
n_convolution = n_wavelet+n_data-1;
n_conv_pow2 = pow2(nextpow2(n_convolution));
half = (n_wavelet-1)/2;

% get FFT of data
eegfft = fft(x, n_conv_pow2);

% initialize
eegpower = zeros(num_frex, n_data); % frequencies X time X trials


% loop through frequencies and compute synchronization
for ii=1:num_frex
    
    wavelet = fft( sqrt(1/(s(ii)*sqrt(pi))) * exp(2*1i*pi*frex(ii).*wave_time) .* exp(-wave_time.^2./(2*(s(ii)^2))) , n_conv_pow2 );
    
    % convolution
    eegconv = ifft(wavelet.*eegfft);
    eegconv = eegconv(1:n_convolution);
    eegconv = eegconv(half+1:end-half);
    
    % Average power over trials (this code performs baseline transform,
    % which you will learn about in chapter 18)
    eegpower(ii,:) = 10*log10(abs(eegconv).^2);
end


%%
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
% leg = legend({'Eye close', 'Eye open'} , 'fontsize', 16)
set(gca, 'plotboxaspectratio', [3 1 1], 'xlim', [time(1) time(end)], 'fontsize', 14);
xlabel('Time (s)', 'fontsize', 16);
ylabel('EEG (normalized)', 'fontsize', 16)
% title('Oz (bandpass filtered [3 Hz, 30 Hz])')
box on;
set(gcf, 'color', 'w');

figure(2); 
contourf(time,frex,eegpower,40,'linecolor','none')
set(gca, 'clim', [-20 10], 'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),...
    'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10,...
    'ydir', 'normal')
% colormap(gray)
ylim([min_freq max_freq])
xlabel('Time (s)', 'fontsize', 16);
ylabel('Frequency (Hz)', 'fontsize', 16);
set(gca, 'plotboxaspectratio', [3 1 1], 'fontsize', 14);
set(gcf, 'color', 'w');
