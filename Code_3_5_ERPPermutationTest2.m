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

NstandardTrials = 1700;
NoddballTrials = 192;

Nresampling_all = round(logspace(1, 5, 10));
Niterations = 30;

tic;
m_vnull = zeros(length(Nresampling_all), Niterations);
std_vnull = zeros(length(Nresampling_all), Niterations); 
for mm = 1:length(Nresampling_all)
    Nresampling = Nresampling_all(mm);
    for nn=1:Niterations
        vnull = ...
            perm_stats_amp(Nresampling, Ntrials, NstandardTrials, NoddballTrials, [xstandard; xoddball], tspan);
        m_vnull(mm,nn) = mean(vnull);
        std_vnull(mm,nn) = std(vnull);
    end
    toc
    pause(0.1)
end
toc
beep
%%
figure(2); clf;
% plot(log10(Nresampling_all), m_vnull', 'k.');

boxplot(m_vnull', log10(Nresampling_all));
figure(3); clf;
plot(log10(Nresampling_all), m_vnull', 'k.');

% violinplot(m_vnull)

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