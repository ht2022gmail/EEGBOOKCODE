clear all; close all; clc;
eeglab; close; % Adding EEGLAB path:

% Load epoched data:
EEG = pop_loadset( 'filename', 'eeglab_data_PREP_preprocessed_epoched.set', 'filepath', './data/');
EEG = eeg_checkset( EEG );
x_epoched = EEG.data;
times_epoched = EEG.times;

% Load continuous data:
EEG = pop_loadset( 'filename', 'eeglab_data_PREP_preprocessed.set', 'filepath', './data/');
EEG = eeg_checkset( EEG );
x_continous = EEG.data;
times_continous = EEG.times;


%% plot continuous data:
figure(1); clf; hold on;
plot_time_channels(times_continous/1000, x_continous, 'k');
ymin = -50; ymax = EEG.nbchan*50;
win = round(100/EEG.srate);
latency = cell2mat({EEG.event.latency});
square_idx = find(strcmp({EEG.event.type}, 'square'));
latency = latency(square_idx);
for nn=1:length(latency)
    t0 = times_continous(round(latency(nn)))/1000;
    plot([t0 t0],[-80 32*80], 'r');
    x0 = [t0-0.5 t0+0.5 t0+0.5 t0-0.5 t0-0.5];
    y0 = [ymin ymin ymax ymax ymin];
    f = fill(x0, y0, 'r');
    set(f, 'FaceAlpha', 0.2, 'EdgeColor', 'w');
end
axis([99 112 ymin ymax]);
xlabel('Time (s)');
ylabel('EEG (\mu V)');
set(gca, 'PlotBoxAspectRatio', [2 1 1], 'YTick', []);
set(gcf, 'color', 'w');

%% plot epoched data:
for nn=1:7
    figure(nn+1); clf; hold on;
    plot_time_channels(times_epoched/1000, x_epoched(:,:,nn), 'k');
    plot([0 0],[-50 32*50], 'r');
    axis([-0.5 0.5 ymin ymax]);
    xlabel('Time (s)');
    ylabel('EEG (\mu V)');
    set(gca, 'PlotBoxAspectRatio', [1 1 1], 'YTick', []);
    set(gcf, 'color', 'w');
end

%% 
function plot_time_channels(times, x, clr)
    for nn=1:size(x,1)
        offset = (nn-1)*50;
        plot(times, x(nn,:)+offset, 'Color', clr);
    end
end


