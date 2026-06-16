clear all; close all; clc;
load('./data/sampledataEEG.mat');
EEG  = ALLEEG(2); % loading epoched data

%%
X = double(EEG.data);
[Nelectrodes, Ntime, Ntrials] = size(X);

% Select Fz channnel for example:
elec_name = 'Pz';
elec = find(strcmp({EEG.chanlocs.labels}, elec_name));
x = squeeze(X(elec,:,:));
X = reshape(x, 1, Ntime*Ntrials); 

% Define Molret wavelet:
time = -1:1/EEG.srate:1;
freq = 4;
s = 5/(2*pi*freq);
molret = exp((-time.^2)/(2*s^2)).*exp(2*pi*1i*freq*time);
 
n_conv = length(X) + length(molret)-1;
half = ceil(length(molret)/2);
 
% frequency-based convolution (convolution theorem):
fft_X = fft(X, n_conv);
fft_molret = fft(molret, n_conv);
ift   = ifft(fft_X.*fft_molret, n_conv);
conv_x = reshape(ift(half:(end-half+1)), [], Ntrials);
% conv_x = real(ift(half:end-half+1)); % cut the initial and last segments.

ITPC = abs(mean(exp(1i*angle(conv_x)),2));

figure
plot(EEG.times, ITPC, 'k');
xlim([-200 1500])
