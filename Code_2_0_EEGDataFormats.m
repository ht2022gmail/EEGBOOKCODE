clear all; close all; clc;
eeglab; close; % Adding EEGLAB path:
load('./data/sampledataEEG.mat');
EEG  = ALLEEG(2); % loading epoched data
addpath('functions\');

%%
X = EEG.data;
[Nelectrodes, Ntime, Ntrials] = size(X);

% Select Pz channnel for example:
elec_name = 'Pz';
trial = 1;
elec = find(strcmp({EEG.chanlocs.labels}, elec_name));
xelec = squeeze(X(elec, :, :)); % time x trials
xtrial = squeeze(X(:, :, trial)); % electrodes x time

%% Vector
% time series:
figure(1); clf;
plot(EEG.times, xelec(:,1), 'k');
xlim([0 1000]);
xlabel('Time (ms)');
ylabel('EEG (\mu V)');
set(gca, 'PlotBoxAspectRatio', [2 1 1]);
set(gcf, 'color', 'w');

% topomap
t0idx = dsearchn(EEG.times', 300);
figure(2); clf; hold on; 
currentx = mean(xtrial(:,t0idx-5:t0idx+5), 2);
topoplot(currentx, EEG.chanlocs);
title([num2str(EEG.times(t0idx)) 'ms'], 'fontsize', 16)
set(gcf, 'color', 'w');

% power spectrum:
x0 = xelec(:);
win = EEG.srate * 3;    % 3-second window
overlap = 0.5*win;  % 50% overlap 
nfft = win;         % # samples for FFT
[px_welch, f_welch, cl_welch] = ...
    pwelch(x0, win, overlap, nfft, EEG.srate, 'ConfidenceLevel', 0.95);
figure(3); 
plot(f_welch, 10*log10(px_welch), 'k');
xlim([0 40])
xlabel('Frequency (Hz)');
ylabel('10 $\log_{10}$ (power)', 'Interpreter','latex');
set(gca, 'PlotBoxAspectRatio', [2 1 1]);
set(gcf, 'color', 'w');

%% Matrix:
% time x trial
figure(4); clf; hold on;
plot_time_trial(EEG.times, xelec);

% time x elecrode
figure(5); clf; hold on;
plot_time_electrode(EEG.times, xtrial);

% time x electrode (contourf)
% figure(6); clf; hold on;
% contourf(EEG.times, 1:Ntrials, xelec', 'linecolor', 'none');
% axis tight;
% set(gca, 'PlotBoxAspectRatio', [2 1 1], ...
%     'YTick', []);
% xlabel('Time (ms)');
% ylabel('Electrode');
% set(gcf, 'color', 'w');
% colorbar

% time x frequency:
figure(6); clf;
[tfmap, freqs] = timefrequencymap(xelec, EEG.times/1000, EEG.pnts, EEG.trials, 2, 60, 30);
plot_time_frequency(EEG.times, tfmap, freqs);

%% Tensor:
% elecrode x time x trial
elec_FPz = find(strcmp({EEG.chanlocs.labels}, 'FPz'));
elec_Fz = find(strcmp({EEG.chanlocs.labels}, 'Fz'));
elec_Pz = find(strcmp({EEG.chanlocs.labels}, 'Pz'));
elec_Cz = find(strcmp({EEG.chanlocs.labels}, 'Cz'));
elec_Oz = find(strcmp({EEG.chanlocs.labels}, 'Oz'));
xelec_FPz = squeeze(X(elec_FPz, :, :));
xelec_Fz = squeeze(X(elec_Fz, :, :));
xelec_Pz = squeeze(X(elec_Pz, :, :));
xelec_Cz = squeeze(X(elec_Cz, :, :));
xelec_Oz = squeeze(X(elec_Oz, :, :));
figure; hold on; plot_time_trial(EEG.times, xelec_FPz);
figure; hold on; plot_time_trial(EEG.times, xelec_Fz);
figure; hold on; plot_time_trial(EEG.times, xelec_Pz);
figure; hold on; plot_time_trial(EEG.times, xelec_Cz);
figure; hold on; plot_time_trial(EEG.times, xelec_Oz);

% time x frequency x electrodeP
[tfmap_FPz, freqs] = timefrequencymap(xelec_FPz, EEG.times/1000, EEG.pnts, EEG.trials, 2, 60, 30);
[tfmap_Fz, freqs] = timefrequencymap(xelec_Fz, EEG.times/1000, EEG.pnts, EEG.trials, 2, 60, 30);
[tfmap_Cz, freqs] = timefrequencymap(xelec_Cz, EEG.times/1000, EEG.pnts, EEG.trials, 2, 60, 30);
[tfmap_Pz, freqs] = timefrequencymap(xelec_Pz, EEG.times/1000, EEG.pnts, EEG.trials, 2, 60, 30);
[tfmap_Oz, freqs] = timefrequencymap(xelec_Oz, EEG.times/1000, EEG.pnts, EEG.trials, 2, 60, 30);
figure; plot_time_frequency(EEG.times, tfmap_FPz, freqs);
figure; plot_time_frequency(EEG.times, tfmap_Fz, freqs);
figure; plot_time_frequency(EEG.times, tfmap_Cz, freqs);
figure; plot_time_frequency(EEG.times, tfmap_Pz, freqs);
figure; plot_time_frequency(EEG.times, tfmap_Oz, freqs);

%% plot_time_trial function:
function plot_time_trial(times, x)
    for nn=1:size(x,2)
        offset = (nn-1)*100;
        plot(times, x(:,nn)+offset, 'k');
    end
    xlim([0 1000]);
    ylim([-100 100*size(x,2)]);
    set(gca, 'PlotBoxAspectRatio', [2 1 1], ...
        'YTick', []);
    xlabel('Time (ms)');
    ylabel('Trial');
    box on;
    set(gcf, 'color', 'w');
end

%% plot_time_electrode function:
function plot_time_electrode(times, x)
    for nn=1:size(x,1)
        offset = (nn-1)*100;
        plot(times, x(nn,:)+offset, 'k');
    end
    xlim([0 1000]);
    ylim([-100 100*size(x,1)]);
    set(gca, 'PlotBoxAspectRatio', [2 1 1], ...
        'YTick', []);
    xlabel('Time (ms)');
    ylabel('Electrode');
    box on;
    set(gcf, 'color', 'w');
end

%% plot_time_frequency function:
function plot_time_frequency(times, tfmap, freqs)
    contourf(times, freqs, tfmap, 40, 'linecolor', 'none')
    set(gca,'clim',[-3 3],'xlim',[-00 1000],'yscale','log','ytick',logspace(log10(2),log10(60),6), ...
        'yticklabel',round(logspace(log10(2),log10(60),6)*10)/10)
    set(gca, 'PlotBoxAspectRatio', [2 1 1]);
    xlabel('Time (ms)');
    ylabel('Frequency (Hz)');
    set(gcf, 'color', 'w');
end