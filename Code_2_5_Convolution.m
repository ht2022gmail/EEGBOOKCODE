clear all; close all; clc;
load('./data/sampledataEEG.mat');
EEG  = ALLEEG(1); % loading continuous data

%% Time-domain and frequency-domain convolution
% Select an electrode for example:
elec_name = 'Oz';
elec = find(strcmp({EEG.chanlocs.labels}, elec_name));
x = squeeze(EEG.data(elec,:)); 

freqs = logspace(log10(3), log10(30), 5);
xwaves = zeros(length(freqs), length(x));
for nn=1:length(freqs)
    xwaves(nn,:) = wavelet_transform(x, EEG.srate, freqs(nn));
end

%%
figure(1); clf; hold on;
plot(EEG.times/1000, x, 'k');
plot(EEG.times/1000, xwaves(5,:), 'r', 'linewidth', 2);
xlim([10 20]);
xlabel('Time (s)');
ylabel('EEG (\mu V)');
title(['Electrode ' elec_name ': time-based convolution']);
legend('original EEG', 'wavelet-convolved EEG');
set(gca, 'plotboxaspectratio', [3 1 1]);
set(gcf, 'color', 'w');

%% 
function xwave = wavelet_transform(x, srate, freq)

    % Define Molret wavelet:
    time = -1:1/srate:1;
    % freq = 10;
    s = 5/(2*pi*freq);
    morlret = exp(2*1i*pi*frequency.*time) .* exp(-time.^2./(2*(4/(2*pi*frequency))^2));

    n_conv = length(x) + length(molret) - 1;
    half = ceil(length(molret)/2);
    
    % frequency-based convolution (convolution theorem):
    fft_Xchan = fft(x, n_conv);
    fft_molret = fft(molret, n_conv);
    ift   = ifft(fft_Xchan.*fft_molret, n_conv);
    xwave = real(ift(half:end-half+1)); % cut the initial and last segments.

end
