%% EEG Alpha Waves dataset
% https://zenodo.org/record/2348892
% FP1, FP2, FC5, FC6, FZ, T7, CZ, T8, P7, P3, PZ, P4, P8, O1, Oz, and O2.

%% Loading data and defining some variables:
clear all; close all; clc;
subnum = 0;
% load subject_00.mat;
load(['./data/subject_' num2str(subnum, '%02d') '.mat']);
load("./data/channels_alphawavedataset.mat");
srate = 512;
Fnyquist = srate/2;
X = SIGNAL(:,2:17)'; 
[Nelectrodes, Ntime] = size(X);
Ntime = size(SIGNAL,1);
time = [0:(Ntime-1)]/srate;
% Bandpass filtering:
filterLength = 200;
cutoffFreq_HP = 2; % in Hz
filterFrequencySpread_HP  = 0.3; % Hz +/- the center frequency
transitionWidth_HP= 0.;
ffrequencies_HP   = ...
    [ 1 (1-transitionWidth_HP)*(cutoffFreq_HP-filterFrequencySpread_HP) (cutoffFreq_HP-filterFrequencySpread_HP)  srate/2 ]/(srate/2);
idealResponse_HP  = [ 0 0 1 1 ];
highpassfilter  = firls(filterLength,ffrequencies_HP,idealResponse_HP);
for nn=1:Nelectrodes
    X(nn,:) = filtfilt(highpassfilter, 1, X(nn,:));
end

onsets_eyeclose = find(SIGNAL(:,18)==1);
onsets_eyeopen = find(SIGNAL(:,19)==1);

%% Spatial GED:

% Computing signal-data matrix
idx = [onsets_eyeclose(1):onsets_eyeopen(1) onsets_eyeclose(2):onsets_eyeopen(2) ...
    onsets_eyeclose(3):onsets_eyeopen(3) onsets_eyeclose(4):onsets_eyeopen(4) ...
    onsets_eyeclose(5):onsets_eyeopen(5)];
XS = X(:, idx);

% Computing the reference and signal covariance matrices:
CR = X*X'/size(X,2);
CS = XS*XS'/size(XS,2);

% Solving temporal GED:
[W, Lambda] = eig(CS, CR);
Lambda = diag(Lambda);
[Lambda, idx_W] = sort(Lambda, 'descend');
W = W(:,idx_W);  % weight matrix
Y = W'*X;   % component matrix
A = inv(W)';  % scalp map matrix

%% 
figure(1); clf; hold on;
for kk=1:length(onsets_eyeopen)
    t1 = time(onsets_eyeclose(kk));
    t2 = time(onsets_eyeopen(kk));
    ff = fill([t1 t2 t2 t1 t1], 3.9*[-1 -1 1 1 -1], [0.9 0.9 0.9]);
    set(ff, 'edgecolor', [1 1 1])
end
plot(time, Y(1,:), 'k');
xlim([time(1) time(end)]); ylim([-4 4])
xlabel('Time (s)'); 
ylabel('EEG (normalized)');
set(gca, 'plotboxaspectratio', [3 1 1]);
box on;
set(gcf, 'color', 'w');


figure(2); clf;
topoplot(A(:,1), chanlocs_alphawavedataset, 'style', 'map', 'electrodes', 'ptslabels');
set(gcf, 'Color', 'w');

figure(3); clf;
topoplot(W(:,1), chanlocs_alphawavedataset, 'style', 'map', 'electrodes', 'ptslabels');
set(gcf, 'Color', 'w');