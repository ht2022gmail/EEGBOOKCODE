clear all; close all; clc
load('./data/G00001JiS.mat'); % MMN data

%% Loading MNN data:
X = double(EEG.data);
[Nelectrodes, Ntime, Ntrials] = size(X);
t = EEG.times/1000;

% 1700 standard trials, 192 oddball trials
standardTrials = unique([EEG.event(find(strcmp({EEG.event.type}, '1'))).epoch]);
oddballTrials = unique([EEG.event(find(strcmp({EEG.event.type}, '3'))).epoch]);

m = mean(X(:,:,standardTrials),3)-mean(X(:,:,oddballTrials),3);

figure(1); clf;
pcolor(t, 1:Nelectrodes, m); colorbar;
set(gca, 'yTick', 1:36, 'yTickLabel', {EEG.chanlocs.labels})
set(gcf, 'color', 'w');

%% 
Nperms = 5000;
NstandardTrials = 1700;
NoddballTrials = 192;

tic;
m_null = zeros(Nelectrodes, Ntime, Nperms);
for nn=1:Nperms
    idx = randperm(Ntrials);
    m_null(:,:,nn) = mean(X(:,:,idx(1:NstandardTrials)),3) ...
        - mean(X(:,:,idx((NstandardTrials+1):end)),3);
end
toc
%%
mean_m_null = mean(m_null,3);
std_m_null = std(m_null,[],3);

z = (m-mean_m_null)./std_m_null;
z_null = (m_null-repmat(mean_m_null,1,1,Nperms))...
    ./repmat(std_m_null,1,1,Nperms);

figure(2); clf;
pcolor(t, 1:Nelectrodes, z); colorbar;
xlabel('Time (s)'); ylabel('Electrodes');
set(gca, 'yTick', 1:36, 'yTickLabel', {EEG.chanlocs.labels})
set(gcf, 'color', 'w');


z_max = max(squeeze(reshape(z_null, Nelectrodes*Ntime,1,Nperms)));
z_min = min(squeeze(reshape(z_null, Nelectrodes*Ntime,1,Nperms)));

th_max = prctile(z_max, 97.5);
th_min = prctile(z_min, 2.5);

mask = (z>=th_max | z<=th_min);
figure(3); clf;
pcolor(t, 1:Nelectrodes, mask.*z); colorbar;
xlabel('Time (s)'); ylabel('Electrodes');
set(gca, 'yTick', 1:36, 'yTickLabel', {EEG.chanlocs.labels})
set(gcf, 'color', 'w');