clear all; close all; clc
load('./data/G00001JiS.mat'); % MMN data

%% Loading MNN data:
X = EEG.data;
[Nelectrodes, Ntime, Ntrials] = size(X);

% Select Fz channnel for example:
elec_name = 'Fz';
elec = find(strcmp({EEG.chanlocs.labels}, elec_name));

standardTrials = unique([EEG.event(find(strcmp({EEG.event.type}, '1'))).epoch]);
oddballTrials = unique([EEG.event(find(strcmp({EEG.event.type}, '3'))).epoch]);

%% Trial avaraging (oddball trials):
Xstandard = squeeze(X(elec, :, standardTrials));
Ntrials = [10 30 50 100];
Niters = 100;
ERPs = zeros(Niters, Ntime, length(Ntrials)); 
for kk=1:length(Ntrials)
    for iter=1:Niters
        idx_rand = randperm(length(standardTrials));
        idx_rand = idx_rand(1:Ntrials(kk));
        ERPs(iter, :, kk) = mean(Xstandard(:, idx_rand),2)';
    end
end

mERPs = squeeze(mean(ERPs,1));
sdERPs = squeeze(std(ERPs,[],1));

%% Visualizing the results:
xx = [EEG.times/1000 EEG.times(end:(-1):1)/1000 EEG.times(1)/1000];

figure(1); clf; hold on;
% plot(EEG.times/1000, squeeze(ERPs(:, :, 1)), 'color', [0.8 0.8 0.8]);
area = [mERPs(:,1)-2*sdERPs(:,1); mERPs(end:(-1):1,1)+2*sdERPs(end:(-1):1,1); mERPs(1,1)-2*sdERPs(1,1)];
ff = fill(xx, area, [0.85 0.85 0.85]); 
set(ff, 'edgecolor', 'w');
plot(EEG.times/1000, mERPs(:, 1)', 'k', 'linewidth', 2);
axis([-0.1 0.4 -5 5]);
xlabel('Time (s)', 'fontsize', 16); ylabel('EEG (\mu V)', 'fontsize', 16);
title(['Averaging ' num2str(Ntrials(1)) ' trials'], 'fontsize', 16);
set(gca, 'plotboxaspectratio', [2 1 1], 'xtick', [-0.1:0.1:0.5], 'ytick', [-5:5]);
set(gcf, 'color', 'w');

figure(2); clf; hold on;
% plot(EEG.times/1000, squeeze(ERPs(:, :, 2)), 'color', [0.8 0.8 0.8]);
area = [mERPs(:,2)-2*sdERPs(:,2); mERPs(end:(-1):1,2)+2*sdERPs(end:(-1):1,2); mERPs(1,2)-2*sdERPs(1,2)];
ff = fill(xx, area, [0.85 0.85 0.85]); 
set(ff, 'edgecolor', 'w');
plot(EEG.times/1000, mERPs(:, 2)', 'k', 'linewidth', 2);
axis([-0.1 0.4 -5 5]);
xlabel('Time (s)', 'fontsize', 16); ylabel('EEG (\mu V)', 'fontsize', 16);
title(['Averaging ' num2str(Ntrials(2)) ' trials'], 'fontsize', 16);
set(gca, 'plotboxaspectratio', [2 1 1], 'xtick', [-0.1:0.1:0.5], 'ytick', [-5:5]);
set(gcf, 'color', 'w');

figure(3); clf; hold on;
% plot(EEG.times/1000, squeeze(ERPs(:, :, 3)), 'color', [0.8 0.8 0.8]);
area = [mERPs(:,3)-2*sdERPs(:,3); mERPs(end:(-1):1,3)+2*sdERPs(end:(-1):1,3); mERPs(1,3)-2*sdERPs(1,3)];
ff = fill(xx, area, [0.85 0.85 0.85]); 
set(ff, 'edgecolor', 'w');
plot(EEG.times/1000, mERPs(:, 3)', 'k', 'linewidth', 2);
axis([-0.1 0.4 -5 5]);
xlabel('Time (s)', 'fontsize', 16); ylabel('EEG (\mu V)', 'fontsize', 16);
title(['Averaging ' num2str(Ntrials(3)) ' trials'], 'fontsize', 16);
set(gca, 'plotboxaspectratio', [2 1 1], 'xtick', [-0.1:0.1:0.5], 'ytick', [-5:5]);
set(gcf, 'color', 'w');

figure(4); clf; hold on;
% plot(EEG.times/1000, squeeze(ERPs(:, :, 4)), 'color', [0.8 0.8 0.8]);
area = [mERPs(:,4)-2*sdERPs(:,4); mERPs(end:(-1):1,4)+2*sdERPs(end:(-1):1,4); mERPs(1,4)-2*sdERPs(1,4)];
ff = fill(xx, area, [0.85 0.85 0.85]); 
set(ff, 'edgecolor', 'w');
plot(EEG.times/1000, mERPs(:, 4)', 'k', 'linewidth', 2);
axis([-0.1 0.4 -5 5]);
xlabel('Time (s)', 'fontsize', 16); ylabel('EEG (\mu V)', 'fontsize', 16);
title(['Averaging ' num2str(Ntrials(4)) ' trials'], 'fontsize', 16);
set(gca, 'plotboxaspectratio', [2 1 1], 'xtick', [-0.1:0.1:0.5], 'ytick', [-5:5]);
set(gcf, 'color', 'w');

%% Scaling law:

Ntrials = [5:200];
Niters = 100;
ERPs = zeros(Niters, Ntime, length(Ntrials)); 
for kk=1:length(Ntrials)
    for iter=1:Niters
        idx_rand = randperm(length(standardTrials));
        idx_rand = idx_rand(1:Ntrials(kk));
        ERPs(iter, :, kk) = mean(Xstandard(:, idx_rand),2)';
    end
end

mERPs = squeeze(mean(ERPs,1));
sdERPs = squeeze(std(ERPs,[],1));

figure(5); clf; 
% plot(log10(Ntrials), log10(mean(sdERPs)), 'k');
loglog(Ntrials, mean(sdERPs), 'k', 'linewidth', 3);
hold on;
p = polyfit(log10(Ntrials), log10(mean(sdERPs)), 1);
loglog(Ntrials, Ntrials.^p(1)*10^p(2), 'r--', 'linewidth', 2);
xlabel('Trials');
ylabel('Standard deviation (\mu V)'); 
xlim([0 200]); ylim([0 3]);
set(gca, 'plotboxaspectratio', [1.5 1 1], 'xtick', [10 50:50:200], 'ytick', [0:1:3]);
leg = legend('data', ['sd = $ 5.34 \times K^{-0.52}$']);
set(leg, 'interpreter', 'latex', 'fontsize', 14);
set(gcf, 'color', 'w');