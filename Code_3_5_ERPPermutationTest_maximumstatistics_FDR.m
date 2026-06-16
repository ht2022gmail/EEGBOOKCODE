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
Nperms = 10000;
NstandardTrials = 1700;
NoddballTrials = 192;

x = [xstandard; xoddball];

m = mean(xstandard) - mean(xoddball);
m_null = zeros(Nperms, Ntime);
for nn=1:Nperms
    idx = randperm(Ntrials);
    m_null(nn,:) = ...
        mean(x(idx(1:NstandardTrials),:)) - mean(x(idx((NstandardTrials+1):end),:));    
end

% z-transform
mean_m_null = mean(m_null);
std_m_null = std(m_null);

z = (m-mean_m_null)./std_m_null;
z_null = (m_null-repmat(mean_m_null,Nperms,1))./repmat(std_m_null,Nperms,1);

%% maximum statistics:
z_max = max(z_null,[],2);
z_min= min(z_null,[],2);

th_max = prctile(z_max, 97.5);
th_min = prctile(z_min, 2.5);

figure(1); clf; hold on; 
histogram(z_max, 'Normalization', 'pdf'); 
histogram(z_min, 'Normalization', 'pdf');
plot([th_min th_min], [0 1], 'k--', 'linewidth', 2);
plot([th_max th_max], [0 1], 'k--', 'linewidth', 2);

figure(2); clf; hold on;
plot(t, z, 'k', 'LineWidth', 2);
plot([min(t) max(t)], [th_min th_min], 'k--');
plot([min(t) max(t)], [th_max th_max], 'k--');

sig_idx = find(z<=th_min | z>=th_max);
plot(t(sig_idx), z(sig_idx), 'r', 'LineWidth', 2);

%% FDR (BH):

p = 1-normcdf(z);
[psorted, idxsorted] = sort(p, 'ascend');

% plot(psorted)

[pID, pN] = fdr(p, 0.05)
sig_idx_FDR = find(p<=pN);

figure(3); clf;
semilogy(t, p); hold on
semilogy([min(t) max(t)], [pN pN], 'b--')
semilogy([min(t) max(t)], [pID pID], 'k--')