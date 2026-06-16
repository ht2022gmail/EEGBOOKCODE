clear all; close all; clc
addpath('./functions')
load('./data/G00001JiS.mat'); % MMN data

%% 
X = double(EEG.data);
[Nelectrodes, Ntime, Ntrials] = size(X);

% for nn=1:Ntrials
%     X(:,:,nn) = bsxfun(@minus, X(:,:,nn), mean(X(:,:,nn),2));
% end

standardTrials = unique([EEG.event(find(strcmp({EEG.event.type}, '1'))).epoch]);
oddballTrials = unique([EEG.event(find(strcmp({EEG.event.type}, '3'))).epoch]);
Xcont = cat(3, X(:,:,standardTrials), X(:,:,oddballTrials));
Xcont = reshape(Xcont, Nelectrodes, Ntime*Ntrials, []);

%% Time-domain and frequency-domain convolution

D = [repmat(zeros(Ntime), 1, length(standardTrials)) repmat(eye(Ntime), 1, length(oddballTrials))];

A = Xcont*pinv(D);
% figure(1); clf;
% plot(EEG.times, A);

CS = (A*D)*(A*D)'/Ntime/Ntrials;
CR = Xcont*Xcont'/Ntime/Ntrials;

[U, Lambda] = eig(CS, CR);
Lambda = diag(Lambda);
[Lambda, idx] = sort(Lambda, 'descend');
U = U(:,idx);

Ystandard = U'*mean(X(:,:,standardTrials),3);
Yoddball = U'*mean(X(:,:,oddballTrials),3);

figure(1); clf; hold on;
plot(EEG.times, Ystandard(1,:), 'k', 'LineWidth', 2);
plot(EEG.times, Yoddball(1,:), 'r', 'LineWidth', 2);
xlim([-100 500]);
set(gca, 'plotboxaspectratio', [2 1 1]);
xlabel('Time (ms)');
ylabel('EEG (normalized)');
set(gcf, 'color', 'w')

figure(2); clf; hold on;
plot(EEG.times, Ystandard(2,:), 'k', 'LineWidth', 2);
plot(EEG.times, Yoddball(2,:), 'r', 'LineWidth', 2);
xlim([-100 500]);
set(gca, 'plotboxaspectratio', [2 1 1]);
xlabel('Time (ms)');
ylabel('EEG (normalized)');
set(gcf, 'color', 'w')

figure(3); clf; hold on;
plot(EEG.times, Ystandard(3,:), 'k', 'LineWidth', 2);
plot(EEG.times, Yoddball(3,:), 'r', 'LineWidth', 2);
xlim([-100 500]);
set(gca, 'plotboxaspectratio', [2 1 1]);
xlabel('Time (ms)');
ylabel('EEG (normalized)');
set(gcf, 'color', 'w')

%%
Astandard = Ystandard(2,:)*mean(X(:,:,standardTrials),3)';
Aoddball = Yoddball(2,:)*mean(X(:,:,oddballTrials),3)';
figure
topoplot(Astandard, EEG.chanlocs)
figure
topoplot(Aoddball, EEG.chanlocs)