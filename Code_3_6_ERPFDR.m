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
ve = mean(xstandard)-mean(xoddball);

Nresampling = 5000;
NstandardTrials = 1700;
NoddballTrials = 192;

vnull = perm_stats_amp_time(Nresampling, Ntrials, NstandardTrials, NoddballTrials, [xstandard; xoddball]);

%% Computing p-values:

% z-transform-based p-value
ze = (ve-mean(vnull))./std(vnull);
pvals = 1-normcdf(ze);
[pvals_sorted, idx_sorted] = sort(pvals, 'ascend'); 

M = length(pvals); % # of hypotheses
alpha_FDR = 0.05;

% Benjamini-Hochberg:
qBH = [1:M]/M*alpha_FDR;
kBH = find(pvals_sorted<=qBH, 1, 'last' );
idx_BH = sort(idx_sorted(1:kBH),'ascend'); % indices for statistical significant hypotheses.

% Benjamini-Yekutieli:
qBY = [1:M]/M/sum(1./[1:M])*alpha_FDR;
kBY = find(pvals_sorted<=qBY, 1, 'last' );
idx_BY = sort(idx_sorted(1:kBY),'ascend');

% qq plot:
figure(2); clf; hold on;
plot(-log10([1:M]/M), -log10(qBH), 'k--');
plot(-log10([kBH+1:M]/M), -log10(pvals_sorted(kBH+1:M)), 'ko', 'MarkerFaceColor', 'k');
plot(-log10([1:kBH]/M), -log10(pvals_sorted(1:kBH)), 'ko', 'MarkerSize', 5);
xlim([0 log10(M)]);
xlabel(' $ -\log_{10} \left( q_i \right) $', 'Interpreter', 'latex');
ylabel(' $-\log_{10} \left(  p_{(i)} \right) $', 'Interpreter','latex');
set(gca, 'PlotBoxAspectRatio', [1.5 1 1]);
set(gcf, 'Color', 'w');

figure(3); clf; hold on;
plot(-log10([1:M]/M), -log10(qBY), 'k--');
plot(-log10([kBY+1:M]/M), -log10(pvals_sorted(kBY+1:M)), 'ko', 'MarkerFaceColor', 'k');
% plot(-log10([1:kBH]/M), -log10(pvals_sorted(1:kBH)), 'ko', 'MarkerSize', 5);
plot(-log10([1:kBY]/M), -log10(pvals_sorted(1:kBY)), 'ko', 'MarkerSize', 5);
xlim([0 log10(M)]);
xlabel(' $ -\log_{10} \left( q_i \right) $', 'Interpreter', 'latex');
ylabel(' $-\log_{10} \left(  p_{(i)} \right) $', 'Interpreter','latex');
set(gca, 'PlotBoxAspectRatio', [1.5 1 1]);
set(gcf, 'Color', 'w');

%
figure(4); clf; hold on;
plot(t, -log10(pvals), 'k');
plot(t(idx_BH), -log10(pvals(idx_BH)), 'r', 'LineWidth', 2);
plot([min(t) max(t)], -log10([qBH(kBH) qBH(kBH)]), 'r--')
axis([min(t) max(t) 0 10]);
xlabel('Time (s)'); ylabel('$ -\log_{10} \left( p \right) $', 'Interpreter', 'latex');
set(gca, 'PlotBoxAspectRatio', [2 1 1]);
set(gcf, 'Color', 'w');

figure(5); clf; hold on;
plot(t, -log10(pvals), 'k');
plot(t(idx_BY), -log10(pvals(idx_BY)), 'r', 'LineWidth', 2);
plot([min(t) max(t)], -log10([qBY(kBY) qBY(kBY)]), 'r--')
axis([min(t) max(t) 0 10]);
xlabel('Time (s)'); ylabel('$ -\log_{10} \left( p \right) $', 'Interpreter', 'latex');
set(gca, 'PlotBoxAspectRatio', [2 1 1]);
set(gcf, 'Color', 'w');

%% function for permutation test:
function vnull = perm_stats_amp_time(Nresampling, Ntrials, NstandardTrials, NoddballTrials, x)
    vnull = zeros(Nresampling, size(x,2));
    for nn=1:Nresampling
        perm_idx = randperm(Ntrials);
        perm_standardTrials = perm_idx(1:NstandardTrials);
        perm_oddballTrials = perm_idx(end-NoddballTrials+1:end);
        perm_xstandard = x(perm_standardTrials,:);
        perm_xoddball = x(perm_oddballTrials,:);
        vnull(nn,:) = mean(perm_xstandard)-mean(perm_xoddball);
    end
end