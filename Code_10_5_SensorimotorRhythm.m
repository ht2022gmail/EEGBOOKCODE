% This code adopted from the following article:
% Delorme, A., Kothe, C., Vankov, A., Bigdely-Shamlo, N., Oostenveld, R., Zander, T. O., & Makeig, S. (2010). MATLAB-based tools for BCI research. In Brain-computer interfaces (pp. 241-259). Springer, London.
% https://sccn.ucsd.edu/wiki/Minimalist_BCI

% Data downloaded from:
% BCI Competition III data set IVb
% http://bbci.de/competition/iii/
% http://bbci.de/competition/iii/
% data_set_ivb 
% data_set_IVb_al_train.mat , data_set_IVb_al_test.mat, true_labels.mat

clear all; close all; clc;
load('./data/data_set_IVb_al_train');


% training: 
flt = @(f)(f>7&f<30).*(1-cos((f-(7+30)/2)/(7-30)*pi*4));

EEG = single(cnt);
Fs = nfo.fs;
mrk = sparse(1,mrk.pos,(mrk.y+3)/2);
wnd0 = [0.5 3.5];
nof = 3;
n = 200;

% frequency filtering and temporal filter estimation
[t,c] = size(EEG); idx = reshape(1:t*c-mod(t*c,n),n,[]);
FLT = real(ifft(fft(EEG).*repmat(flt(Fs*(0:t-1)/t)',1,c)));
T = FLT(idx)/EEG(idx);

% data epoching, class-grouping and CSP
wnd = round(Fs*wnd0(1)):round(Fs*wnd0(2));
for k = 1:2
    EPO{k} = FLT(repmat(find(mrk==k),length(wnd),1) + repmat(wnd',1,nnz(mrk==k)),:);
end
[V,D] = eig(cov(EPO{2}),cov(EPO{1})+cov(EPO{2}));
W = V(:,[1:nof end-nof+1:end]);
 
% log-variance feature extraction
for k = 1:2
    X{k} = squeeze(log(var(reshape(EPO{k}*W, length(wnd),[],2*nof))));
end
% LDA wight and bias:
w = ((mean(X{2})-mean(X{1}))/(cov(X{1})+cov(X{2})))';
b = (mean(X{1})+mean(X{2}))*w/2;

tt = 0:1/Fs:(t-1)/Fs;

%% Visualizing EEG data:
eC3 = 52; eC4 = 56;
C3 = EEG(:,eC3);
C3fft = fft(C3);

figure(1); clf; 
subplot(1,3,1);
plot(Fs*(0:t-1)/t,abs(C3fft), 'k');
xlim([0 40]); ylim([0 10^6])
set(gca, 'plotboxaspectratio', [2 1 1]);
xlabel('Frequency (Hz)'); ylabel('Fourier coefficient (\mu V)')

subplot(1,3,2);
plot(Fs*(0:t-1)/t,abs(C3fft).*flt(Fs*(0:t-1)/t)', 'k');
xlim([0 40]); ylim([0 10^6])
set(gca, 'plotboxaspectratio', [2 1 1]);
xlabel('Frequency (Hz)'); ylabel('Fourier coefficient (\mu V)')

subplot(1,3,3);
plot(Fs*(0:t-1)/t,flt(Fs*(0:t-1)/t), 'k');
xlim([0 40]); ylim([0 3])
set(gca, 'plotboxaspectratio', [2 1 1]);
xlabel('Frequency (Hz)'); ylabel('Filter')

set(gcf, 'color', 'w');

figure(2); clf; hold on;

subplot(2,1,1);
plot(tt, EEG(:,eC3), 'k');
xlim([150 160]);
title('C3 (raw)')
xlabel('Time (s)'); ylabel('EEG (\mu V)');
set(gca, 'plotboxaspectratio', [3 1 1]);

subplot(2,1,2);
plot(tt, FLT(:,eC3), 'k');
xlim([150 160]);
title('C3 (bandpass-filtered)')
xlabel('Time (s)'); ylabel('EEG (\mu V)');
set(gca, 'plotboxaspectratio', [3 1 1]);

set(gcf, 'color', 'w');

%% Visualizing CSPs:
class1 = find(mrk==1); class2 = find(mrk==2);
tt = 0:1/Fs:(t-1)/Fs;

CSP = FLT*W;
figure(3); clf; hold on;
std6 = std(CSP(:,1));
XX = [0 300 300 0 0]; YY = [1 1 -1 -1 1]*std6*4; 
for k=1:length(class1)
    ff1 = fill((XX+class1(k))/Fs, YY, [1.0 0.4 0.4]); set(ff1, 'edgecolor', [1 1 1]);
    ff2 = fill((XX+class2(k))/Fs, YY, [0.0 0.4 1.0]); set(ff1, 'edgecolor', [1 1 1]);
end
plot(tt, CSP(:,1), 'k');
xlim([620 670]); ylim([-3 3]);
xlabel('Time (s)'); ylabel('CSP1');
set(gca, 'plotboxaspectratio', [3 1 1]);
set(gcf, 'color', 'w');

figure(4); clf; hold on;
std6 = std(CSP(:,6));
XX = [0 300 300 0 0]; YY = [1 1 -1 -1 1]*std6*5; 
for k=1:length(class1)
    ff1 = fill((XX+class1(k))/Fs, YY, [1.0 0.4 0.4]); set(ff1, 'edgecolor', [1 1 1]);
    ff2 = fill((XX+class2(k))/Fs, YY, [0.0 0.4 1.0]); set(ff1, 'edgecolor', [1 1 1]);
end
plot(tt, CSP(:,6), 'k');
xlim([620 670]); ylim([-3 3]);
xlabel('Time (s)'); ylabel('CSP6');
set(gca, 'plotboxaspectratio', [3 1 1]);
set(gcf, 'color', 'w');

%% Comparing electrode and CSPs:
% class1 = find(mrk==1); class2 = find(mrk==2);
% eC3 = 52;
% eC4 = 56;
% XX = [0 300 300 0 0]; YY0 = [1 1 -1 -1 1]; 
% 
% figure(5); clf; hold on;
% subplot(3,1,1); hold on;
% for k=1:length(class1)
%     ff1 = fill((XX+class1(k))/Fs, YY0*max(abs(EEG(:,eC3))), [1.0 0.4 0.4]);
%     set(ff1, 'edgecolor', [1 1 1]);
%     ff2 = fill((XX+class2(k))/Fs, YY0*max(abs(EEG(:,eC3))), [0.0 0.4 1.0]);
%     set(ff1, 'edgecolor', [1 1 1]);
% end
% plot(tt, EEG(:,eC3), 'k-'); xlim([250 350]); ylim([-1000 1000])
% set(gca, 'plotboxaspectratio', [5 1 1]);
% xlabel('time (s)'); ylabel('EEG (\mu V)');
% title('C3 (raw)');
% 
% subplot(3,1,2); hold on;
% for k=1:length(class1)
%     ff1 = fill((XX+class1(k))/Fs, YY0*max(abs(FLT(:,eC3))), [1.0 0.4 0.4]); 
%     set(ff1, 'edgecolor', [1 1 1]);
%     ff2 = fill((XX+class2(k))/Fs, YY0*max(abs(FLT(:,eC3))), [0.0 0.4 1.0]); 
%     set(ff1, 'edgecolor', [1 1 1]);
% end
% plot(tt, FLT(:,eC3), 'k-'); xlim([250 350]); ylim([-300 300])
% set(gca, 'plotboxaspectratio', [5 1 1]);
% xlabel('time (s)'); ylabel('EEG (\mu V)');
% title('C3 (bandpass filtered)');
% 
% FLT1 = FLT*W;
% subplot(3,1,3); hold on;
% for k=1:length(class1)
%     ff1 = fill((XX+class1(k))/Fs, YY0*max(abs(FLT1(:,6))), [1.0 0.4 0.4]); 
%     set(ff1, 'edgecolor', [1 1 1]);
%     ff2 = fill((XX+class2(k))/Fs, YY0*max(abs(FLT1(:,6))), [0.0 0.4 1.0]); 
%     set(ff1, 'edgecolor', [1 1 1]);
% end
% plot(tt, FLT1(:,6), 'k-'); xlim([250 350]); ylim([-3 3])
% set(gca, 'plotboxaspectratio', [5 1 1]);
% xlabel('time (s)'); ylabel('EEG (normalized)');
% title('CSP 6 (raw)');
% set(gcf, 'color', 'w');


%% Visualizing the feature vectors:
EPO1 = EPO{1}*W;
EPO2 = EPO{2}*W;
xx1 = zeros(6,105);
xx2 = zeros(6,105);

for n=1:105
    trialIndex = (length(wnd))*(n-1)+1:(length(wnd))*n;
    xx1(:,n) = var(EPO1(trialIndex,:),[],1)';
    xx2(:,n) = var(EPO2(trialIndex,:),[],1)';
end

figure(6); clf; hold on;
plot(xx1(1,:), xx1(6,:), 'kx');
plot(xx2(1,:), xx2(6,:), 'r+');
xlabel('Var(CSP1)');
ylabel('Var(CSP6)');
set(gca, 'plotboxaspectratio', [1 1 1]);
box on
set(gcf, 'color', 'w');

figure(7); clf; hold on;
plot(log(xx1(1,:)), log(xx1(6,:)), 'kx', 'markersize', 10);
plot(log(xx2(1,:)), log(xx2(6,:)), 'r+', 'markersize', 10);
xlabel('log(Var(CSP1))', 'fontsize', 16);
ylabel('log(Var(CSP6))', 'fontsize', 16);
% leg
set(gca, 'plotboxaspectratio', [1 1 1]);
box on
set(gcf, 'color', 'w');

w16 = -w(1)/w(6);
b16 = b/w(6);
x16 = linspace(-2, 0.5, 100);
plot(x16, w16*x16+b16, 'b--');
axis([-2 0.5 -4 1]);
set(gca, 'xtick', -2:1, 'ytick', -4:1);
