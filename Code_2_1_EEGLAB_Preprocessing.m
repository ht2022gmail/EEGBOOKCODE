clear all; close all; clc;
% eeglab; close; % Adding EEGLAB path:

% Load raw data:
EEG = pop_loadset( 'filename', 'eeglab_data.set', 'filepath', './data/');
EEG = eeg_checkset( EEG );
rawdata = EEG.data;

% Load PREP preprocessed data:
EEG = pop_loadset( 'filename', 'eeglab_data_PREP_preprocessed.set', 'filepath', './data/');
EEG = eeg_checkset( EEG );
PREPdata = EEG.data;

%%
figure(1); clf; hold on;
plot_time_channels(EEG.times/1000, rawdata, 'k');
plot_time_channels(EEG.times/1000, PREPdata, 'r');
xlim([100 110])
xlabel('Time (ms)');
ylabel('EEG (\mu V)');
set(gca, 'PlotBoxAspectRatio', [2 1 1], 'YTick', []);
set(gcf, 'color', 'w');

[px_rawdata, f_welch] = compute_average_power_spectrum(rawdata, EEG.srate, 0.5);
[px_PREPdata, f_welch] = compute_average_power_spectrum(PREPdata, EEG.srate, 0.5);
figure(2); clf; hold on;
plot(f_welch, mean(10*log10(px_rawdata),2), 'k', 'LineWidth', 2);
plot(f_welch, mean(10*log10(px_PREPdata),2), 'r', 'LineWidth', 2);
axis([0 64 -20 30])
xlabel('Frequency (Hz)');
ylabel('10 $\log_{10}$ (power)', 'Interpreter','latex');
set(gca, 'PlotBoxAspectRatio', [1.5 1 1]);
set(gcf, 'color', 'w');

%% 
function plot_time_channels(times, x, clr)

    for nn=1:size(x,1)
        offset = (nn-1)*100;
        plot(times, x(nn,:)+offset, 'Color', clr);
    end

end

function [px_welch_all, f_welch] = compute_average_power_spectrum(x, srate, overlap)
    px_welch_all = [];
    win = srate * 3;    % 3-second window
    overlap = overlap*win;  % 50% overlap 
    nfft = win;         % # samples for FFT
    for nn=1:size(x,1)
        x0 = x(nn, :);
        [px_welch, f_welch] = ...
            pwelch(x0, win, overlap, nfft, srate, 'ConfidenceLevel', 0.95);
        px_welch_all =[px_welch_all px_welch];
    end
end
