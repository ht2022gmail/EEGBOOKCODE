clear all; close all; clc;

fns = {'Kuhlman1970_Figure6', 'Trimper2013_Figure1', 'Ouedraogo2016_Figure4', 'Tort2010_Figure5'};

srate = 1000;
dt = 1/srate;

for nn=1:length(fns)
    a = readtable([fns{nn} '.csv']);
    a = table2array(a);
    [~, idx] = sort(a(:,1), 'ascend');

    % sort time in ascending order:
    t0 = a(idx, 1);
    x0 = a(idx, 2);

    % remove duplicate values:
    [t0, idx] = unique(t0);
    x0 = x0(idx);

    % interpolate to 1000 Hz:
    t = t0(1):dt:t0(end); % sampling in 1000 Hz
    x = interp1(t0, x0, t, 'linear'); % interpolate x0 values at the new time vector t
    
    t = t-t(1);
    x = x - mean(x);

    
    figure(nn); 
    plot(t,x,'k', 'LineWidth',2);
    xlim([t(1) t(end)]);
    title(fns{nn})
    set(gca, 'PlotBoxAspectRatio', [4 1 1])
    set(gcf, 'color', 'w');
    xlabel('Time (s)');
    ylabel('EEG (\mu V)')
    box off

endCode_5_8_NonSinusoidalWaveForms.m