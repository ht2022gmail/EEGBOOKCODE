clear all; close all; clc;
load('./data/sampledataEEG.mat');
EEG  = ALLEEG(2); % loading epoched data

%% Loading MNN data:
X = EEG.data;
[Nelectrodes, Ntime, Ntrials] = size(X);

% Select Fz channnel for example:
elec_name = 'Fz';
elec = find(strcmp({EEG.chanlocs.labels}, elec_name));

standardTrials = unique([EEG.event(find(strcmp({EEG.event.type}, '1'))).epoch]);
oddballTrials = unique([EEG.event(find(strcmp({EEG.event.type}, '3'))).epoch]);
x = squeeze(reshape(X(elec, :, :), 1, 1, Ntime*Ntrials))';

%%
min_freq =  2;
max_freq = 80;
num_frex = 30;

% define wavelet parameters
time = -1:1/EEG.srate:1;
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
s = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);
% s = 3./(2*pi*frex); % this line is for figure 13.14
% s = 10./(2*pi*frex); % this line is for figure 13.14

% definte convolution parameters
n_wavelet            = length(time);
n_data               = EEG.pnts*EEG.trials;
n_convolution        = n_wavelet+n_data-1;
n_conv_pow2          = pow2(nextpow2(n_convolution));
half_of_wavelet_size = (n_wavelet-1)/2;

% get FFT of data
eegfft = fft(x, n_conv_pow2);

% initialize
eegpower = zeros(num_frex, EEG.pnts); % frequencies X time X trials
eegpower_trials = zeros(num_frex, EEG.pnts, EEG.trials);
baseidx = dsearchn(EEG.times', [-100 -50]');

% loop through frequencies and compute synchronization
for fi=1:num_frex
    
    wavelet = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2 );
    
    % convolution
    eegconv = ifft(wavelet.*eegfft);
    eegconv = eegconv(1:n_convolution);
    eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
    
    % Average power over trials (this code performs baseline transform,
    % which you will learn about in chapter 18)
    eegpower_trials(fi,:,:) = abs(reshape(eegconv,EEG.pnts,EEG.trials)).^2;
    baselinepower = squeeze(mean(eegpower_trials(fi, baseidx(1):baseidx(2), :),2));
    eegpower_trials(fi,:,:) = 10*log10(squeeze(eegpower_trials(fi,:,:))./repmat(baselinepower', EEG.pnts, 1));

    temppower = mean(abs(reshape(eegconv,EEG.pnts,EEG.trials)).^2,2);
    eegpower(fi,:) = 10*log10(temppower./mean(temppower(baseidx(1):baseidx(2))));
end

figure(1); clf;
contourf(EEG.times,frex,eegpower,40,'linecolor','none')
set(gca,'clim', [-3 3], 'xlim',[-200 1000],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
xlabel('Time (ms)');
ylabel('Frequency *Hz');
% title('Logarithmic frequency scaling')
set(gca, 'plotboxaspectratio', [2 1 1]);
set(gcf, 'color', 'w');
colormap(gray); colorbar;
% 
figure(2); clf;
for kk=1:6
    subplot(3,2,kk);
    contourf(EEG.times,frex,eegpower_trials(:,:,kk),40,'linecolor','none')
    set(gca, 'xlim',[-200 1000],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
    xlabel('Time (ms)');
    ylabel('Frequency *Hz');
    title(['Trial ' num2str(kk)])
    set(gca, 'plotboxaspectratio', [2 1 1]);
    colormap(gray);
end
set(gcf, 'color', 'w');

%% CANDECOMP
addpath('./toolboxes/nway320');

model = parafac(eegpower_trials, 2, [0 0 0 10 1000], [2 0 2]);

[frequencies, times, trials] = fac2let(model);

figure(3); clf;
contourf(EEG.times,frex,frequencies(:,1)*times(:,1)',40,'linecolor','none');
set(gca, 'xlim',[-200 1000],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
xlabel('Time (ms)');
ylabel('Frequency *Hz');
title('Component 1');
set(gca, 'plotboxaspectratio', [2 1 1]);
colormap(gray); colorbar;

figure(4); clf;
contourf(EEG.times,frex,frequencies(:,2)*times(:,2)',40,'linecolor','none');
set(gca, 'xlim',[-200 1000],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
xlabel('Time (ms)');
ylabel('Frequency *Hz');
title('Component 2');
set(gca, 'plotboxaspectratio', [2 1 1]);
colormap(gray); colorbar;
