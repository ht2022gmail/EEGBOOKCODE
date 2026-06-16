clear all; close all; clc
rng();

%% Three generative models of ERPs:
t = [-500:1000]/1000;
T = length(t);
mask_task = (t>=0 & t<=0.2);
mask_nontask = (t<0 | t>0.2);
mask_postrask = (t>0);
s = sin(2*pi*5*t)'.*mask_task';
Ntrials = 100;

% Additive model;
x_additive = repmat(s, 1, Ntrials) + 20*pinknoise(T, Ntrials); 

% Phase-reset model:
phi = 2*pi*rand(Ntrials, 1) - pi;
background = sin(2*pi*5*t+phi)'.*mask_nontask';
x_phasereset = s + background + 10*pinknoise(T, Ntrials);

% baseline-shift model:
phi = 2*pi*rand(Ntrials, 1) - pi;
background = sin(2*pi*5*t+phi)'+1;
x_baselineshift = (2+s).*background + (-2) + 10*pinknoise(T, Ntrials);


%% Visualizing ERPs:
figure(1); clf; hold on;
plot(t, x_additive, 'color', [0.8 0.8 0.8]);
plot(t, mean(x_additive,2), 'k', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('EEG (\mu V)');
xlim([-0.2 0.4]); ylim([-3 5])
set(gca, 'PlotBoxAspectRatio', [1.5 1 1]);
set(gcf, 'color', 'w');

figure(2); clf; hold on;
plot(t, x_phasereset, 'color', [0.8 0.8 0.8]);
plot(t, mean(x_phasereset,2), 'k', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('EEG (\mu V)');
xlim([-0.2 0.4]); ylim([-3 5])
set(gca, 'PlotBoxAspectRatio', [1.5 1 1]);
set(gcf, 'color', 'w');

figure(3); clf; hold on;
plot(t, x_baselineshift, 'color', [0.8 0.8 0.8]);
plot(t, mean(x_baselineshift,2), 'k', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('EEG (\mu V)');
xlim([-0.2 0.4]); ylim([-3 5])
set(gca, 'PlotBoxAspectRatio', [1.5 1 1]);
set(gcf, 'color', 'w');

% ERP images:
figure(4); clf;
imagesc(t, 1:Ntrials, x_additive');
set(gca, 'ydir', 'reverse', 'clim', [-3 5], 'PlotBoxAspectRatio', [1.5 1 1]);
xlabel('Time (s)'); ylabel('Trial');
xlim([-0.2 0.4]);
colormap(hot)
set(gcf, 'color', 'w');

figure(5); clf;
imagesc(t, 1:Ntrials, x_phasereset');
set(gca, 'ydir', 'reverse', 'clim', [-3 5], 'PlotBoxAspectRatio', [1.5 1 1]);
xlabel('Time (s)'); ylabel('Trial');
xlim([-0.2 0.4]);
colormap(hot)
set(gcf, 'color', 'w');

figure(6); clf;
imagesc(t, 1:Ntrials, x_baselineshift');
set(gca, 'ydir', 'reverse', 'clim', [-3 5], 'PlotBoxAspectRatio', [1.5 1 1]);
xlabel('Time (s)'); ylabel('Trial');
xlim([-0.2 0.4]);
colormap(hot)
set(gcf, 'color', 'w');