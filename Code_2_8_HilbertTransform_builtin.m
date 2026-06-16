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
% Using the MATLAB builtin "bandpass" function:
x = bandpass(Xchan, [8 12], EEG.srate);

% Perform Hilbert transform:
% Or use MATLAB build-in function:
z = hilbert(x);

%% Plot the results:
figure(1); clf; hold on;
plot(EEG.times/1000, real(z), 'k');
plot(EEG.times/1000, imag(z), 'r');
xlim([100 110])
xlabel('Time (s)');
ylabel('EEG (\mu V)')
legend('Real', 'Imaginary')
set(gca, 'plotboxaspectratio', [2 1 1]);
set(gcf, 'color', 'w');

figure(2); clf; hold on;
plot3(EEG.times/1000, real(z), imag(z), 'b')
plot3(EEG.times/1000, real(z), -50*ones(1,Ntime), 'k')
plot3(EEG.times/1000, 50*ones(1,Ntime), imag(z), 'r')
xlim([100 110])
set(gca, 'plotboxaspectratio', [3 1 1]);
xlabel('Time (s)');
ylabel('real(EEG) (\mu V)')
zlabel('imag(EEG) (\mu V)')
view([-20 35])
set(gcf, 'color', 'w');

figure(3); clf;
subplot(1,2,1);
plot(EEG.times/1000, abs(z), 'k');
xlabel('Time (s)');
ylabel('Amplitude (\mu V)')
title('Amplitude $| z(t) |$', 'Interpreter', 'latex');
set(gca, 'plotboxaspectratio', [2 1 1]);
xlim([100 105])

subplot(1,2,2);
plot(EEG.times/1000, angle(z), 'k');
xlabel('Time (s)');
ylabel('Angle (rad)')
title('Phase $\mathrm{arg}(z(t))$', 'Interpreter', 'latex');
set(gca, 'plotboxaspectratio', [2 1 1]);
xlim([100 105])
set(gcf, 'color', 'w');


