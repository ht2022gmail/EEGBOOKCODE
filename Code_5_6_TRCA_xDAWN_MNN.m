clear all; close all; clc
addpath('./functions')
load('./data/G00001JiS.mat'); % MMN data

%% 
X = double(EEG.data);
[Nelectrodes, Ntime, Ntrials] = size(X);

standardTrials = unique([EEG.event(find(strcmp({EEG.event.type}, '1'))).epoch]);
oddballTrials = unique([EEG.event(find(strcmp({EEG.event.type}, '3'))).epoch]);
Xcont = squeeze(reshape(X, Nelectrodes, Ntime*Ntrials, []));

%% TRCA
tau = Ntime;
times = EEG.times;

% standard trials
tidx_standard = (standardTrials-1)*Ntime+1;
[Wstandard, Lambdastandard, Ystandard, Xbstandard] = TRCA(Xcont, tidx_standard, tau);
[~, ybstandard] = blocky(Xcont, Wstandard(:,1), tidx_standard, tau);
ystandard = mean(ybstandard); 
ystandard = ystandard - mean(ystandard);
ystandard = ystandard/std(ystandard);
astandard = mean(Xbstandard,3)*ystandard'/Ntime;

% oddball trials
tidx_oddball = (oddballTrials-1)*Ntime+1;
[Woddball, Lambdaoddball, Yoddball, Xboddball] = TRCA(Xcont, tidx_oddball, tau);
[~, yboddball] = blocky(Xcont, Woddball(:,1), tidx_oddball, tau);
yoddball = mean(yboddball);
yoddball = yoddball - mean(yoddball);
yoddball = yoddball / std(yoddball);
aoddball = mean(Xboddball,3)*yoddball'/Ntime;

% Plotting
figure(1); clf; hold on
plot(times, ystandard,'k','LineWidth',2); 
plot(times, yoddball,'r','LineWidth',2);
xlabel('Time (ms)');
ylabel('Task-related component (a.u.)');
set(gca, 'PlotBoxAspectRatio', [3 1 1]);
set(gcf, 'color', 'w');

figure(2); clf; 
topoplot(astandard, EEG.chanlocs, 'style', 'map', 'electrodes', 'ptslabels');
set(gcf, 'Color', 'w')

figure(3); clf; 
topoplot(aoddball, EEG.chanlocs, 'style', 'map', 'electrodes', 'ptslabels');
set(gcf, 'Color', 'w')


%% xDAWN

% standard trials
[Ustandard, ~, Astandard] =  xDAWN(Xcont, tidx_standard, tau);
Zstandard = Ustandard'*mean(X(:,:,standardTrials),3);
zstandard = Zstandard(1,:) - mean(Zstandard(1,:));
zstandard = zstandard/std(zstandard);
astandard = mean(X(:,:,standardTrials),3)*zstandard'/size(Xcont,2);

% oddball trials
[Uoddball, ~, Aoddball] =  xDAWN(Xcont, tidx_oddball, tau);
Zoddball = Uoddball'*mean(X(:,:,oddballTrials),3);
zoddball = Zoddball(1,:) - mean(Zoddball(1,:));
zoddball = zoddball/std(zoddball);
aoddball = mean(X(:,:,oddballTrials),3)*zoddball'/size(Xcont,2);

% Plotting
figure(4); clf; hold on
plot(times, zoddball,'r','LineWidth',2);
plot(times, -zstandard, 'k', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('xDAWN component (a.u.)');
set(gca, 'PlotBoxAspectRatio', [3 1 1]);
set(gcf, 'color', 'w');

figure(5); clf; 
topoplot(-astandard, EEG.chanlocs, 'style', 'map', 'electrodes', 'ptslabels');
set(gcf, 'Color', 'w')

figure(6); clf; 
topoplot(aoddball, EEG.chanlocs, 'style', 'map', 'electrodes', 'ptslabels');
set(gcf, 'Color', 'w')