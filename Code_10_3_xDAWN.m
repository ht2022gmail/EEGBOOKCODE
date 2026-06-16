%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                          %%%%%%
%%%%%         sample classification ofthe P300 test data       %%%%%%
%%%%%                                                          %%%%%%
%%%%%                BCI Competition III Challenge             %%%%%%
%%%%%                                                          %%%%%%
%%%%%        (C) Dean Krusienski and Gerwin Schalk 2004        %%%%%%
%%%%%                 Wadsworth Center/NYSDOH                  %%%%%%
%%%%%                                                          %%%%%%
%%%%%               function calls: topoplotEEG.m              %%%%%%
%%%%%                                                          %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all; clc
addpath('./functions/');
load('./data/Subject_A_Train.mat') % load data file
% load 'Subject_B_Train.mat' % load data file

% convert to double precision
Signal=double(Signal);
Flashing=double(Flashing);
StimulusCode=double(StimulusCode);
StimulusType=double(StimulusType);

%%
fs = 240; % Sampling frequency
tau = 100;

trial = 1;
t = [0:(tau-1)]/fs;
X = squeeze(Signal(trial, :, :))';
X = bsxfun(@minus, X, mean(X,2));
X = bsxfun(@rdivide, X, std(X,[],2));
[Nchannels, Ntime] = size(X);

% xdawn filter for target:
onsets_target = onset_detection(StimulusType(trial, :));
Ntargets = sum(onsets_target>0);
% conditions = StimulusCode(trial, onsets_target);
% target_character = TargetChar(trial);
% column_row = sort(unique(conditions), 'ascend');
[U_target, A_target, S_target, Lambda_target] = xdawn_filter(X, onsets_target, tau);

% xdawn filter for distractors:
onsets_distractors = onset_detection(StimulusCode(trial, :));
onsets_distractors(onsets_target==1) = 0;
Ndistractors = sum(onsets_distractors>0);
[U_distractors, A_distractors, S_distractors, Lambda_distractors] = xdawn_filter(X, onsets_distractors, tau);

%%
Xepoch_target = epoch(X, find(onsets_target), tau);
Xepoch_distractors = epoch(X, find(onsets_distractors), tau);

f_target = zeros(Nchannels, Ntargets);
for kk=1:Ntargets
    f_target(:,kk) = diag((U_target'*Xepoch_target(:,:, kk))*(S_target)')/tau;
end
f_distractors = zeros(Nchannels, Ndistractors);
for kk=1:Ndistractors
    f_distractors(:,kk) = diag((U_distractors'*Xepoch_distractors(:,:, kk))*(S_distractors)')/tau;
end


%% Visualizing the results:

figure(1); clf; hold on;
plot(t, S_distractors(1, :)', 'k', 'LineWidth', 2);
plot(t, S_target(1, :)', 'r', 'LineWidth', 2);
xlim([t(1) t(end)]);
xlabel('Time (s)', 'FontSize', 14); ylabel('xDAWN Component', 'FontSize', 14);
set(gca, 'plotboxaspectratio', [1.5 1 1]);
box on;
set(gcf, 'color', 'w');

figure(2); clf; hold on;
plot(t, S_distractors(2, :)', 'k', 'LineWidth', 2);
plot(t, S_target(2, :)', 'r', 'LineWidth', 2);
xlim([t(1) t(end)]);
xlabel('Time (s)', 'FontSize', 14); ylabel('xDAWN Component', 'FontSize', 14);
set(gca, 'plotboxaspectratio', [1.5 1 1],'xtick', 0:0.1:1, 'ytick', -1:0.2:1);
box on;set(gcf, 'color', 'w');

figure(3); clf; hold on;
plot(t, S_distractors(3, :)', 'k', 'LineWidth', 2);
plot(t, S_target(3, :)', 'r', 'LineWidth', 2);
xlim([t(1) t(end)]);
xlabel('Time (s)', 'FontSize', 14); ylabel('xDAWN Component', 'FontSize', 14);
set(gca, 'plotboxaspectratio', [1.5 1 1]);
box on;set(gcf, 'color', 'w');

figure(4); clf; hold on
plot(f_target(1,:), f_target(2,:), 'rx', 'MarkerSize', 10, 'linewidth' ,1.5);
plot(f_distractors(1,:), f_distractors(2,:), 'k.', 'MarkerSize', 15);
xlabel('xDAWN component 1', 'FontSize', 14); 
ylabel('xDAWN component 2', 'FontSize', 14); 
set(gca, 'plotboxaspectratio', [1 1 1], 'xtick', -1:0.5:2, 'YTick', -1:0.5:1);
box on;
set(gcf, 'color', 'w');

% figure(5); clf; hold on
% plot3(f_target(1,:), f_target(2,:), f_target(3,:), 'rx', 'MarkerSize', 10, 'linewidth' ,1.5);
% plot3(f_distractors(1,:), f_distractors(2,:), f_distractors(3,:), 'k.', 'MarkerSize', 15);
% xlabel('xDAWN component 1', 'FontSize', 14); 
% ylabel('xDAWN component 2', 'FontSize', 14); 
% zlabel('xDAWN component 3', 'FontSize', 14); 
% set(gca, 'plotboxaspectratio', [1 1 1], 'xtick', -1:0.5:2, 'YTick', -1:0.5:1);
% box on;
% set(gcf, 'color', 'w');

%% LDA
mu_target = mean(f_target(1:2,:), 2);
mu_distractors = mean(f_distractors(1:2,:), 2);

C_target = cov(f_target(1:2,:)');
C_distractors = cov(f_distractors(1:2,:)');

w = inv(C_target+C_distractors)*(mu_target+mu_distractors);
b = w'*(mu_target+mu_distractors)/2;

figure(4);
x0 = linspace(-1, 1.5, 100);
plot(x0, -(w(1)*x0-b)/w(2), 'k--')
% plot(mu_target(1), mu_target(2), 'rx', 'markersize', 20);
% plot(mu_distractors(1), mu_distractors(2), 'kx', 'markersize', 20);