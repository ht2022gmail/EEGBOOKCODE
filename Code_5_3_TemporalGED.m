%% EEG Alpha Waves dataset
% https://zenodo.org/record/2348892
% FP1, FP2, FC5, FC6, FZ, T7, CZ, T8, P7, P3, PZ, P4, P8, O1, Oz, and O2.

%% Loading data and defining some variables:
clear all; close all; clc;
subnum = 0;
% load subject_00.mat;
load(['./data/subject_' num2str(subnum, '%02d') '.mat']);
srate = 512;
Fnyquist = srate/2;
Ntime = size(SIGNAL,1);
time = [0:(Ntime-1)]/srate;
% eeg_oz = SIGNAL(:,16); 
eeg_oz = SIGNAL(:,17); 
eeg_oz = detrend(eeg_oz);
x = eeg_oz;

onsets_eyeclose = find(SIGNAL(:,18)==1);
onsets_eyeopen = find(SIGNAL(:,19)==1);

%% Temporal GED
p = 512/4+1; % Corresponding to 500 ms window.
T = length(x);
timewin = [-p:p]/srate;

% Computing reference-data matrix:
X = zeros(2*p+1, T-2*p);
for tt=1:(T-2*p)
    X(:,tt) = [x(tt:tt+(2*p))'];
end

% Computing signal-data matrix
idx = [onsets_eyeclose(1):onsets_eyeopen(1) onsets_eyeclose(2):onsets_eyeopen(2) ...
    onsets_eyeclose(3):onsets_eyeopen(3) onsets_eyeclose(4):onsets_eyeopen(4) ...
    onsets_eyeclose(5):onsets_eyeopen(5)];
idx = idx - p;
XS = X(:,idx);

% Computing the reference and signal covariance matrices:
CR = X*X'/size(X,2);
CS = XS*XS'/size(XS,2);

% Solving temporal GED:
[W, Lambda] = eig(CS, CR);
Lambda = diag(Lambda);
[Lambda, idx] = sort(Lambda, 'descend');
W = W(:,idx);
Y = W'*X;

%% Visualizing the results:
figure(1); clf; hold on;
for kk=1:length(onsets_eyeopen)
    t1 = time(onsets_eyeclose(kk));
    t2 = time(onsets_eyeopen(kk));
    ff = fill([t1 t2 t2 t1 t1], 108*[-1 -1 1 1 -1], [0.9 0.9 0.9]);
    set(ff, 'edgecolor', [1 1 1])
end
plot(time, x, 'k');
xlim([time(p+1) time(T-p)]); ylim([-110 110])
xlabel('Time (s)');
ylabel('EEG (\mu V)');
set(gca, 'plotboxaspectratio', [3 1 1]);
box on;
set(gcf, 'color', 'w');

% PC1
figure(2); clf; hold on;
for kk=1:length(onsets_eyeopen)
    t1 = time(onsets_eyeclose(kk));
    t2 = time(onsets_eyeopen(kk));
    ff = fill([t1 t2 t2 t1 t1], 4.9*[-1 -1 1 1 -1], [0.9 0.9 0.9]);
    set(ff, 'edgecolor', [1 1 1])
end
plot(time(p+1:T-p), Y(1,:), 'k');
xlim([time(p+1) time(T-p)]); ylim([-5 5]);
rectangle('Position', [20 -5 5 10], 'EdgeColor', 'r')
xlabel('Time (s)'); 
ylabel('EEG (normalized)');
set(gca, 'plotboxaspectratio', [3 1 1]);
box on;
set(gcf, 'color', 'w');

% PC2
figure(3); clf; hold on;
for kk=1:length(onsets_eyeopen)
    t1 = time(onsets_eyeclose(kk));
    t2 = time(onsets_eyeopen(kk));
    ff = fill([t1 t2 t2 t1 t1], 4.9*[-1 -1 1 1 -1], [0.9 0.9 0.9]);
    set(ff, 'edgecolor', [1 1 1])
end
plot(time(p+1:T-p), Y(2,:), 'k');
rectangle('Position', [20 -5 5 10], 'EdgeColor', 'r')
xlim([time(p+1) time(T-p)]); ylim([-5 5]);
xlabel('Time (s)'); 
ylabel('EEG (normalized)');
set(gca, 'plotboxaspectratio', [3 1 1]);
box on;
set(gcf, 'color', 'w');

% PC1 & PC2;
figure(4); clf; hold on;
plot(time(p+1:T-p), Y(1,:), 'k');
plot(time(p+1:T-p), Y(2,:), 'r');
xlim([20 25]); ylim([-5 5]);
xlabel('Time (s)'); 
ylabel('EEG (normalized)');
legend('PC 1', 'PC 2', 'Fontsize', 20);
set(gca, 'plotboxaspectratio', [3 1 1]);
box on;
set(gcf, 'color', 'w');

figure(5); clf; hold on;
plot(timewin, W(:,1), 'k');
plot(timewin, W(:,2), 'r');
xlim([timewin(1) timewin(end)])
ylim([-2E-3 2E-3])
xlabel('Time (s)');
ylabel('Filter (A.U.)');
legend('PC 1', 'PC 2', 'Fontsize', 20);
set(gca, 'plotboxaspectratio', [3 1 1]);
box on;
set(gcf, 'color', 'w');