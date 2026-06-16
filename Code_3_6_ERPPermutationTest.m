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
x = mean(X(elec, :, oddballTrials),3) - mean(X(elec, :, standardTrials),3);

%% Permutation testing:
Nresampling = 3000; % typically >1000

% Random permutation:
xresample = zeros(Nresampling, Ntime);
for nn=1:Nresampling
    random_idx = randperm(Ntrials);
    random_standardTrials = random_idx(1:length(standardTrials));
    random_oddballTrials = random_idx(length(standardTrials)+1:end);
    xresample(nn,:) = mean(X(elec, :, random_oddballTrials),3) - mean(X(elec, :, random_standardTrials),3);
end

% Computing the threshold:
xresample_mean = mean(xresample);
xresample_std = std(xresample);
threshold = xresample_std*2;

% Compute the null disrtibution of cluster size:
nullDist = zeros(Nresampling,1);

for nn=1:Nresampling
    Cresample = bwconncomp(abs(xresample(nn,:))>=threshold);
    numPixels_resample = cellfun(@numel,Cresample.PixelIdxList);
    if ~isempty(numPixels_resample)
        nullDist(nn) = max(numPixels_resample);
    else
         nullDist(nn) = 0;
    end
end

threshold_cluster = quantile(nullDist, 0.95);
C = bwconncomp(abs(x)>=threshold);
numPixels = cellfun(@numel,C.PixelIdxList);
idx_significant = find(numPixels>threshold_cluster);

%% Visualize the results:
% ERP (trial average):
figure(1); clf; hold on;
for kk=1:length(idx_significant)
    idx_cluster = cell2mat(C.PixelIdxList(idx_significant));
    xfill = EEG.times([idx_cluster(1) idx_cluster(end) idx_cluster(end) idx_cluster(1) idx_cluster(1)]);
    yfill = 4*[-1 -1 1 1 -1];
    ff = fill(xfill, yfill, [1 .8 .8]); 
    set(ff, 'edgecolor', 'none');
end
plot(EEG.times, mean(X(elec, :, standardTrials),3), 'k', 'linewidth', 2);
plot(EEG.times, mean(X(elec, :, oddballTrials),3), 'r', 'linewidth', 2);
axis([0 400 -4.05 2]);
xlabel('Time (ms)')
ylabel('EEG (\mu V)')
title(['Electrode ' EEG.chanlocs(elec).labels])
legend('', 'Standard', 'Oddball');
set(gca, 'plotboxaspectratio', [2 1 1]);
set(gcf, 'color', 'w')

figure(2); clf; hold on;
histogram(nullDist, 'facecolor', [0 0 0]);
plot(threshold_cluster*[1 1], [0 500], 'r--', 'linewidth', 2);
xlim([-1 20])
xlabel('Maximum cluster size')
ylabel('Frequency')
set(gca, 'plotboxaspectratio', [1 1 1]);
set(gcf, 'color', 'w')
