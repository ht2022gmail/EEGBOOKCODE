clear all; close all; clc;
load('./data/sampledataEEG.mat');
EEG  = ALLEEG(2); % loading epoched data

%% Woody's method
X = EEG.data;
[Nelectrodes, Ntime, Ntrials] = size(X);

% Select Pz channnel for example:
elec_name = 'Fz';
elec = find(strcmp({EEG.chanlocs.labels}, elec_name));
Xchan = squeeze(X(elec,:,:)); % time * trials
Xchan0 = Xchan; % Original data 

% Woody's algorithm
lags = zeros(Ntrials,1); % total lags
lags0 = zeros(Ntrials,1);

err = 1;
iter = 1;
while err>1e-3 && iter<1000
    lags_old = lags;
    Xav = mean(Xchan,2);
    for kk=1:Ntrials
        [c, lag] = xcorr(Xav, Xchan(:,kk));
        [~, lagidx] = max(c); % find the lag by maximizing cross correlation
        lags0(kk) = lag(lagidx);
        Xchan(:,kk) = circshift(Xchan(:,kk), lags0(kk)); % circshift by the lag
    end
    lags = lags + lags0;
    err = norm(lags-lags_old)/norm(lags);
    iter = iter + 1;
end

%% Visualization 
figure(1); clf;
subplot(2,2,1);
imagesc(EEG.times, [1:Ntrials], Xchan0');
xlabel('Time (ms)')
ylabel('Trials')
title(['Electrode ' elec_name ', Original data'])
set(gca, 'plotboxaspectratio', [2 1 1])
subplot(2,2,3);
plot(EEG.times, mean(Xchan0,2), 'k')
axis([EEG.times(1) EEG.times(end) -10 50])
xlabel('Time (ms)')
ylabel('EEG (\mu V)')
title(['Electrode ' elec_name ', Original data'])
set(gca, 'plotboxaspectratio', [2 1 1])

subplot(2,2,2);
imagesc(EEG.times, [1:Ntrials], Xchan');
xlabel('Time (ms)')
ylabel('Trials')
title(['Electrode ' elec_name ', Lag-adjusted data'])
set(gca, 'plotboxaspectratio', [2 1 1])
subplot(2,2,4);
plot(EEG.times, mean(Xchan,2), 'k')
axis([EEG.times(1) EEG.times(end) -10 50])
xlabel('Time (ms)')
ylabel('EEG (\mu V)')
title(['Electrode ' elec_name ', Lag-adjusted data'])
set(gca, 'plotboxaspectratio', [2 1 1])

set(gcf, 'color', 'w');