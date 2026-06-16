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
eeg_oz = SIGNAL(:,16); 
eeg_oz = detrend(eeg_oz);
x = eeg_oz;
x = bandpass(x, [1 40], srate);

%% Cut eye-open and eye-close segments:
onsets_eyeclose = find(SIGNAL(:,18)==1);
onsets_eyeopen = find(SIGNAL(:,19)==1);

idx_eyeopen = [1:onsets_eyeclose(1)-1 ...
    onsets_eyeopen(1):onsets_eyeclose(2)-1 ...
    onsets_eyeopen(2):onsets_eyeclose(3)-1 ...
    onsets_eyeopen(3):onsets_eyeclose(4)-1 ...
    onsets_eyeopen(4):onsets_eyeclose(5)-1 ...
    onsets_eyeopen(5):length(x)];

idx_eyeclose = [ ...
    onsets_eyeclose(1):onsets_eyeopen(1)-1 ...
    onsets_eyeclose(2):onsets_eyeopen(2)-1 ...
    onsets_eyeclose(3):onsets_eyeopen(3)-1 ...
    onsets_eyeclose(4):onsets_eyeopen(4)-1 ...
    onsets_eyeclose(5):onsets_eyeopen(5)-1];

x_eyeopen = x(idx_eyeopen);
x_eyeclose = x(idx_eyeclose);

x_eyeopen = x_eyeopen(1:length(x_eyeopen)-mod(length(x_eyeopen), srate));
x_eyeclose = x_eyeclose(1:length(x_eyeclose)-mod(length(x_eyeclose), srate));

% Reshape to time*blocks:
x_eyeopen = reshape(x_eyeopen, srate, length(x_eyeopen)/srate);
x_eyeclose= reshape(x_eyeclose, srate, length(x_eyeclose)/srate);

%% Computing spectra:
f = (0:srate-1); % frequency for one-second period

% computer spectra:
fx_eyeopen = fft(x_eyeopen);
fx_eyeclose = fft(x_eyeclose);
px_eyeopen = 20*log10(abs(fx_eyeopen).^2/srate^2);
px_eyeclose = 20*log10(abs(fx_eyeclose).^2/srate^2);

m_px_eyeopen = mean(px_eyeopen,2);
m_px_eyeclose = mean(px_eyeclose,2);
std_px_eyeopen = std(px_eyeopen,[],2);
std_px_eyeclose = std(px_eyeclose,[],2);

figure(1); clf; hold on;
fill_alpha = fill([8 12 12 8 8], [-25 -25 30 30 -25], [0.9 0.9 0.9], 'EdgeColor','none');
fill_beta = fill([20 30 30 20 20], [-25 -25 30 30 -25], [0.9 0.9 0.9], 'EdgeColor','none');
plot(f, m_px_eyeclose, 'r', 'LineWidth', 2);
plot(f, m_px_eyeopen, 'k', 'LineWidth', 2);
axis([3 32 -25 30])
set(gca, 'PlotBoxAspectRatio', [2 1 1]);
set(gcf, 'Color', 'w');
xlabel('Frequency (Hz)'); ylabel('10 log_{10} (PSD) (dB)');


%% Violin plots:
idx_alpha = dsearchn(f', [8; 12]);
px_eyeopen_alpha = mean(px_eyeopen(idx_alpha,:)); 
px_eyeclose_alpha = mean(px_eyeclose(idx_alpha,:)); 

idx_beta = dsearchn(f', [20; 30]);
px_eyeopen_beta = mean(px_eyeopen(idx_beta,:)); 
px_eyeclose_beta = mean(px_eyeclose(idx_beta,:)); 

figure(2); clf; 
v_alpha = violinplot(categorical({'Open', 'Close'}), [px_eyeopen_alpha; [px_eyeclose_alpha NaN(1,16)]]');
v_alpha(1).FaceColor = [0.3 0.3 0.3];
v_alpha(2).FaceColor = [1 0.3 0.3];
ylabel('10 log_{10} (PSD) (dB)');
set(gca, 'PlotBoxAspectRatio', [2 1 1]);
set(gcf, 'Color', 'w');

figure(3); clf; 
v_beta = violinplot(categorical({'Open', 'Close'}), [px_eyeopen_beta; [px_eyeclose_beta NaN(1,16)]]');
v_beta(1).FaceColor = [0.3 0.3 0.3];
v_beta(2).FaceColor = [1 0.3 0.3];
ylabel('10 log_{10} (PSD) (dB)');
set(gca, 'PlotBoxAspectRatio', [2 1 1]);
set(gcf, 'Color', 'w');

%% Permutation test:
px_alpha = [px_eyeopen_alpha px_eyeclose_alpha];
px_beta = [px_eyeopen_beta px_eyeclose_beta];
Nperms = 10000;
N_eyeopen = length(px_eyeopen_alpha);
N_eyeclose = length(px_eyeclose_alpha);
N = N_eyeclose+N_eyeopen;

v_alpha = mean(px_eyeclose_alpha) - mean(px_eyeopen_alpha);
v_beta = mean(px_eyeclose_beta) - mean(px_eyeopen_beta);
v_alpha_null = zeros(Nperms,1);
v_beta_null = zeros(Nperms,1);
for nn=1:Nperms
    idx = randperm(N);
    v_alpha_null(nn) = mean(px_alpha(idx(1:N_eyeclose))) - mean(px_alpha(idx(N_eyeclose+1:N))) ;    
    v_beta_null(nn) = mean(px_beta(idx(1:N_eyeclose))) - mean(px_beta(idx(N_eyeclose+1:N))) ;    
end

% Calculate the p-value based on the null distribution
p_value_alpha = sum(v_alpha_null >= v_alpha) / Nperms;
p_value_beta = sum(v_beta_null >= v_beta) / Nperms;

figure(4); clf; hold on;
histogram(v_alpha_null, 'Normalization', 'pdf', 'FaceColor', [0.6 0.6 0.6]);
plot(v_alpha, 0, 'rx', 'MarkerSize', 15, 'LineWidth', 2);
xlabel('Power difference (dB)'); ylabel('PDF');
set(gca, 'PlotBoxAspectRatio', [2 1 1]);
set(gcf, 'Color', 'w');

figure(5); clf; hold on;
histogram(v_beta_null, 'Normalization', 'pdf', 'FaceColor', [0.6 0.6 0.6]);
plot(v_beta, 0, 'rx', 'MarkerSize', 15, 'LineWidth', 2);
xlabel('Power difference (dB)'); ylabel('PDF');
set(gca, 'PlotBoxAspectRatio', [2 1 1]);
set(gcf, 'Color', 'w');
