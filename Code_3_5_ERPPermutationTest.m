clear all; close all; clc
load('./data/G00001JiS.mat'); % MMN data

%% Loading MNN data:
X = double(EEG.data);
[Nelectrodes, Ntime, Ntrials] = size(X);
t = EEG.times/1000;

% Select Fz channnel for example:
elec_name = 'Fz';
elec = find(strcmp({EEG.chanlocs.labels}, elec_name));

% 1700 standard trials, 192 oddball trials
standardTrials = unique([EEG.event(find(strcmp({EEG.event.type}, '1'))).epoch]);
oddballTrials = unique([EEG.event(find(strcmp({EEG.event.type}, '3'))).epoch]);

xstandard = squeeze(X(elec,:,standardTrials))';
xoddball = squeeze(X(elec,:,oddballTrials))';

figure(1); clf; hold on;
plot(t, mean(xstandard), 'k');
plot(t, mean(xoddball), 'r');

%% Permutation:
tspan = dsearchn(t', [0.160;0.240]);
ve = mean(mean(xstandard(:,tspan)))-mean(mean(xoddball(:,tspan)));

Nresampling = 1000;
NstandardTrials = 1700;
NoddballTrials = 192;

vnull = perm_stats_amp(Nresampling, Ntrials, NstandardTrials, NoddballTrials, [xstandard; xoddball], tspan);

%% Computing p-values:
% count-based p-value
pval1 = sum(vnull>ve)/Nresampling

figure(2); clf; hold on;
histogram(vnull, 'Normalization','pdf', 'FaceColor', [0.4 0.4 0.4]);
plot(ve, 0, 'r*')
xlabel('stats'); ylabel('PDF');
set(gca, 'PlotBoxAspectRatio', [2 1 1], 'XTick', [-2:2]);
set(gcf, 'color', 'w');

% z-transform-based p-value
ze = (ve-mean(vnull))/std(vnull);
pval2 = 1-normcdf(ze)

figure(3); clf; hold on;
x0 = -5:0.01:5;
y0 = normpdf(x0);
plot(x0, y0, 'k');
plot(ze, 0, 'r*')
xlabel('z-value'); ylabel('PDF');
set(gca, 'PlotBoxAspectRatio', [2 1 1], 'XTick', [-5:2.5:5]);
set(gcf, 'color', 'w');

%% function for permutation test:
function vnull = perm_stats_amp(Nresampling, Ntrials, NstandardTrials, NoddballTrials, x, tspan)
    vnull = zeros(Nresampling, 1);
    for nn=1:Nresampling
        perm_idx = randperm(Ntrials);
        perm_standardTrials = perm_idx(1:NstandardTrials);
        perm_oddballTrials = perm_idx(NstandardTrials+1:Ntrials);
        perm_xstandard = x(perm_standardTrials,:);
        perm_xoddball = x(perm_oddballTrials,:);
        vnull(nn) = mean(mean(perm_xstandard(:,tspan)))-mean(mean(perm_xoddball(:,tspan)));
    end
end