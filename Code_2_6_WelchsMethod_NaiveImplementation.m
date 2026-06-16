clear all; close all; clc;
load('./data/sampledataEEG.mat');
EEG  = ALLEEG(1); % loading continuous data

%% Fourier transform
X = EEG.data;
[Nelectrodes, Ntime, Ntrials] = size(X);

f = [0:(Ntime-1)]/Ntime*EEG.srate;

% Select an electrode for example:
elec_name = 'Oz';
elec = find(strcmp({EEG.chanlocs.labels}, elec_name));
Xchan = squeeze(X(elec,:,:)); % time * trials
Xchan = detrend(Xchan);

fXchan = fft(Xchan);
Xpower = 10*log10(1/(EEG.srate*Ntime)*abs((fXchan).^2)*2);

%% Welch method (naive implementation):
win = EEG.srate * 3;    % 3-second window
overlap = 0.5*win;  % 50% overlap 
nfft = win;         % # samples for FFT
windows = 1:overlap:(length(fXchan)-win);
Nwindows = length(windows);
hann_window = hann(win);

f_welch = [0:(win-1)]/win*EEG.srate;
Xchan_win = zeros(Nwindows, win);
for kk=1:Nwindows
%     Xchan_win(kk,:) = Xchan(windows(kk):windows(kk)+win-1);
    Xchan_win(kk,:) = Xchan(windows(kk):windows(kk)+win-1).*hann_window';
end
fXchan_win = fft(Xchan_win')';
px_welch_win = abs(fXchan_win).^2/(win*EEG.srate); 
px_welch = mean(abs(fXchan_win).^2)/(win*EEG.srate);

%% Visualizing the results:
figure(1); clf; hold on;
% plot(f, Xpower, 'k');
plot(f_welch, 10*log10(px_welch_win), 'color', [0.8 0.8 0.8]);
plot(f_welch, 10*log10(px_welch), 'k', 'linewidth', 2);
xlabel('Frequency (Hz)');
ylabel('log_{10} (EEG power)');
set(gca, 'PlotBoxAspectRatio', [1.5 1 1]);
xlim([0 40]); ylim([-30 30])
box on
set(gcf, 'color', 'w');

%%
figure(2); clf;
plot(EEG.times(1:win)/1000, Xchan_win(5,:), 'k-');
xlabel('Time (s)');
ylabel('EEG (\mu V)');
xlim([0 3]); ylim([-60 60])
set(gca, 'PlotBoxAspectRatio', [1.5 1 1], 'xtick', 0:3, 'ytick', -50:25:50);
set(gcf, 'color', 'w');

figure(3); clf; 
plot(f_welch, 10*log10(px_welch_win(5,:)), 'k-');
xlabel('Frequency (Hz)');
ylabel('log_{10} (EEG power)');
set(gca, 'PlotBoxAspectRatio', [1.5 1 1]);
xlim([0 40]); ylim([-30 30])
set(gcf, 'color', 'w');

figure(4); clf; 
plot(EEG.times(1:win)/1000, Xchan_win(15,:), 'k-');
xlabel('Time (s)');
ylabel('EEG (\mu V)');
xlim([0 3]); ylim([-60 60])
set(gca, 'PlotBoxAspectRatio', [1.5 1 1], 'xtick', 0:3, 'ytick', -50:25:50);
set(gcf, 'color', 'w');

figure(5); clf; 
plot(f_welch, 10*log10(px_welch_win(15,:)), 'k-');
xlabel('Frequency (Hz)');
ylabel('log_{10} (EEG power)');
set(gca, 'PlotBoxAspectRatio', [1.5 1 1]);
xlim([0 40]); ylim([-30 30])
set(gcf, 'color', 'w');

figure(7); clf; 
plot(EEG.times(1:win)/1000, Xchan_win(25,:), 'k-');
xlabel('Time (s)');
ylabel('EEG (\mu V)');
xlim([0 3]); ylim([-60 60])
set(gca, 'PlotBoxAspectRatio', [1.5 1 1], 'xtick', 0:3, 'ytick', -50:25:50);
set(gcf, 'color', 'w');

figure(8); clf; 
plot(f_welch, 10*log10(px_welch_win(25,:)), 'k-');
xlabel('Frequency (Hz)');
ylabel('log_{10} (EEG power)');
set(gca, 'PlotBoxAspectRatio', [1.5 1 1]);
xlim([0 40]); ylim([-30 30])
set(gcf, 'color', 'w');
