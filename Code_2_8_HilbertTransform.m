clear all; close all; clc;
load('./data/sampledataEEG.mat');
EEG  = ALLEEG(1); % loading continuous data

%% Fourier transform
X = double(EEG.data);
[Nelectrodes, Ntime, Ntrials] = size(X);

f = [0:(Ntime-1)]/Ntime*EEG.srate;

% Select an electrode for example:
elec_name = 'Oz';
elec = find(strcmp({EEG.chanlocs.labels}, elec_name));
Xchan = squeeze(X(elec,:,:)); % time * trials

% Perform bandpass filter aroung 10 Hz:
center_freq = 10; % in Hz
filter_frequency_spread  = 2; % Hz +/- the center frequency
transition_width = 0.2;
ffrequencies   = ...
    [ 0 (1-transition_width)*(center_freq-filter_frequency_spread) (center_freq-filter_frequency_spread) (center_freq+filter_frequency_spread) (1+transition_width)*(center_freq+filter_frequency_spread) EEG.srate/2 ]/(EEG.srate/2);
idealresponse  = [ 0 0 1 1 0 0 ];
filterweights  = firls(200,ffrequencies,idealresponse);
Xchan = filtfilt(filterweights, 1, Xchan);

% Perform Hilbert transform:
fXchan = fft(Xchan);
hilbert_filter = [1 2*ones(1,floor(Ntime/2)-1) zeros(1,floor(Ntime/2))];
hilbert_fXchan = hilbert_filter.*fXchan;
hilbert_Xchan = ifft(hilbert_fXchan);

% Or use MATLAB build-in function:
% z = hilbert(Xchan);

%% Plot the results:
figure(4); clf; hold on;
plot(EEG.times/1000, real(hilbert_Xchan), 'k');
plot(EEG.times/1000, imag(hilbert_Xchan), 'r');
xlim([100 110])
xlabel('Time (s)');
ylabel('EEG (\mu V)')
legend('Real', 'Imaginary')
set(gca, 'plotboxaspectratio', [2 1 1]);
set(gcf, 'color', 'w');

figure(5); clf; hold on;
plot3(EEG.times/1000, real(hilbert_Xchan), imag(hilbert_Xchan), 'b')
plot3(EEG.times/1000, real(hilbert_Xchan), -50*ones(1,Ntime), 'k')
plot3(EEG.times/1000, 50*ones(1,Ntime), imag(hilbert_Xchan), 'r')
xlim([100 110])
set(gca, 'plotboxaspectratio', [3 1 1]);
xlabel('Time (s)');
ylabel('real(EEG) (\mu V)')
zlabel('imag(EEG) (\mu V)')
view([-20 35])
set(gcf, 'color', 'w');

figure(6); clf;
subplot(1,2,1);
plot(EEG.times/1000, abs(hilbert_Xchan), 'k');
xlabel('Time (s)');
ylabel('Amplitude (\mu V)')
title('Amplitude $| z(t) |$', 'Interpreter', 'latex');
set(gca, 'plotboxaspectratio', [2 1 1]);
xlim([100 105])

subplot(1,2,2);
plot(EEG.times/1000, angle(hilbert_Xchan), 'k');
xlabel('Time (s)');
ylabel('Angle (rad)')
title('Phase $\mathrm{arg}(z(t))$', 'Interpreter', 'latex');
set(gca, 'plotboxaspectratio', [2 1 1]);
xlim([100 105])
set(gcf, 'color', 'w');


