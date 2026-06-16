clear all; close all; clc;
load('./data/sampledataEEG.mat');
EEG  = ALLEEG(1); % loading continuous data

%% Fourier transform:
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

%% Define EEG segments:
win = EEG.srate * 3;    % 3-second window
overlap = 0.5*win;  % 50% overlap 
nfft = win;         % # samples for FFT
windows = 1:overlap:(length(fXchan)-win);
Nwindows = length(windows);
f_welch = [0:(win-1)]/win*EEG.srate;
Xchan_win = zeros(Nwindows, win);
for kk=1:Nwindows
    Xchan_win(kk,:) = Xchan(windows(kk):windows(kk)+win-1);
end

%% Define window functions:
rect_window = ones(win,1);
hann_window = hann(win);
hamming_window = hamming(win);
gauss_window = gauss(win, 2)';

% Naive Welch's method:
fXchan_win = fft(Xchan_win')';
px_welch_win = abs(fXchan_win).^2/(win*EEG.srate); 
px_welch = mean(abs(fXchan_win).^2)/(win*EEG.srate);

% Hann window:
Xchan_win_hann = Xchan_win.*repmat(hann_window', Nwindows, 1);
fXchan_win_hann = fft(Xchan_win_hann')';
px_welch_win_hann = abs(fXchan_win_hann).^2/(win*EEG.srate); 
px_welch_hann = mean(abs(fXchan_win_hann).^2)/(win*EEG.srate);

% Hamming window:
Xchan_win_hamming = Xchan_win.*repmat(hamming_window', Nwindows, 1);
fXchan_win_hamming = fft(Xchan_win_hamming')';
px_welch_win_hamming = abs(fXchan_win_hamming).^2/(win*EEG.srate); 
px_welch_hamming = mean(abs(fXchan_win_hamming).^2)/(win*EEG.srate);

% Gauss window:
Xchan_win_gauss = Xchan_win.*repmat(gauss_window', Nwindows, 1);
fXchan_win_gauss = fft(Xchan_win_gauss')';
px_welch_win_gauss = abs(fXchan_win_gauss).^2/(win*EEG.srate); 
px_welch_gauss = mean(abs(fXchan_win_gauss).^2)/(win*EEG.srate);

%% Visualizing a single window:

figure(1); clf; hold on;
plot(EEG.times(1:win)/1000, Xchan_win(55,:), 'k-');
xlabel('Time (s)');
ylabel('EEG (\mu V)');
xlim([0 3]); ylim([-50 50])
set(gca, 'PlotBoxAspectRatio', [1.5 1 1], 'xtick', 0:3, 'ytick', -40:20:40);
set(gcf, 'color', 'w');
box on;

% Hann window
figure(2); clf; hold on;
ff_x = [EEG.times(1:win) EEG.times(win:(-1):1) EEG.times(1)]/1000;
ff_y = 50*[-hann_window; hann_window(end:(-1):1); -hann_window(1)]';
ff = fill(ff_x, ff_y, [0.9 0.9 0.9]); 
plot(EEG.times(1:win)/1000, Xchan_win_hann(55,:), 'k-');
plot(EEG.times(1:win)/1000, 50*hann_window, 'k--');
plot(EEG.times(1:win)/1000, -50*hann_window, 'k--');
set(ff, 'edgecolor', 'w');
xlabel('Time (s)');
ylabel('EEG (\mu V)');
xlim([0 3]); ylim([-50 50])
set(gca, 'PlotBoxAspectRatio', [1.5 1 1], 'xtick', 0:3, 'ytick', -40:20:40);
set(gcf, 'color', 'w');
box on

figure(3); clf; 
plot(f_welch, 10*log10(px_welch_win_hann(55,:)), 'k-');
xlabel('Frequency (Hz)');
ylabel('log_{10} (EEG power)');
set(gca, 'PlotBoxAspectRatio', [1.5 1 1]);
xlim([0 40]); ylim([-40 20])
set(gcf, 'color', 'w');

% % Hamming window
% figure(4); clf;
% plot(EEG.times(1:win)/1000, Xchan_win_hamming(55,:), 'k-');
% xlabel('Time (s)');
% ylabel('EEG (\mu V)');
% xlim([0 3]); ylim([-40 40])
% set(gca, 'PlotBoxAspectRatio', [1.5 1 1], 'xtick', 0:3, 'ytick', -50:25:50);
% set(gcf, 'color', 'w');
% 
% figure(5); clf; 
% plot(f_welch, 10*log10(px_welch_win_hamming(55,:)), 'k-');
% xlabel('Frequency (Hz)');
% ylabel('log_{10} (EEG power)');
% set(gca, 'PlotBoxAspectRatio', [1.5 1 1]);
% xlim([0 40]); ylim([-40 20])
% set(gcf, 'color', 'w');
% 
% % Gauss window
% figure(6); clf;
% plot(EEG.times(1:win)/1000, Xchan_win_gauss(55,:), 'k-');
% xlabel('Time (s)');
% ylabel('EEG (\mu V)');
% xlim([0 3]); ylim([-40 40])
% set(gca, 'PlotBoxAspectRatio', [1.5 1 1], 'xtick', 0:3, 'ytick', -50:25:50);
% set(gcf, 'color', 'w');
% 
% figure(7); clf; 
% plot(f_welch, 10*log10(px_welch_win_gauss(55,:)), 'k-');
% xlabel('Frequency (Hz)');
% ylabel('log_{10} (EEG power)');
% set(gca, 'PlotBoxAspectRatio', [1.5 1 1]);
% xlim([0 40]); ylim([-40 20])
% set(gcf, 'color', 'w');

%% Visualizing spectra:
% figure(8); clf; hold on;
% % Hann window
% % plot(f_welch, 10*log10(px_welch), 'k', 'linewidth', 2);
% plot(f_welch, 10*log10(px_welch_hann), 'k', 'linewidth', 2);
% plot(f_welch, 10*log10(px_welch_hamming), 'k', 'linewidth', 2);
% plot(f_welch, 10*log10(px_welch_gauss), 'k', 'linewidth', 2);
% xlabel('Frequency (Hz)');
% ylabel('log_{10} (EEG power)');
% set(gca, 'PlotBoxAspectRatio', [1.5 1 1]);
% xlim([0 40]); ylim([-30 30])
% box on
% set(gcf, 'color', 'w');
