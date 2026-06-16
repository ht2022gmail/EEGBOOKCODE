clear all; close all; clc;
load('./data/sampledataEEG.mat');
EEG  = ALLEEG(1); % loading continuous data

%% Principal component analysis (PCA)
X = double(EEG.data);
[Nelectrodes, Ntime, Ntrials] = size(X);

% Three ways to subtract mean values:
% onevec = ones(Ntime,1);
% X = X*(eye(Ntime)-onevec*onevec'/Ntime); % (this is memory consuming)
% X = X - repmat(mean(X,2), 1, Ntime);
X = bsxfun(@minus, X, mean(X,2));

% PCA by SVD:
[U,S,V] = svd(X, 'econ');
S = diag(S);

% Fourier transform:
f = [0:(Ntime-1)]/Ntime*EEG.srate;
fV = zeros(size(V));
for nn=1:Nelectrodes
    fV(:,nn) = fft(V(:,nn));
end
% fV = fft(V);


%%
figure(1); clf;
for nn=1:3
    subplot(3, 5, 5*(nn-1)+1);
    topoplot(U(:,nn), EEG.chanlocs);
    subplot(3, 5, 5*(nn-1)+[2 3]);
    plot(EEG.times/1000, V(:,nn), 'k');
    xlim([10 20]);
    xlabel('Time (s)'); ylabel('EEG (\mu V)');
    subplot(3, 5, 5*(nn-1)+[4 5]);
    plot(f, 20*log10(abs(fV(:,nn))), 'k');
    xlabel('Frequency (Hz)'); ylabel('10 log (Power)');
    xlim([0 30]); ylim([-30 30])
end
set(gcf, 'color', 'w');

figure(2); clf;
subplot(1, 2, 1); 
semilogy(1:Nelectrodes, S, 'k.-', 'linewidth' ,2, 'markersize', 15);
set(gca, 'plotboxaspectratio', [1.5 1 1]); xlim([0 33]);
xlabel('Order'); ylabel('log(eigenvalue)');

subplot(1, 2, 2); 
plot(cumsum(S)/sum(S), 'k.-', 'linewidth' ,2, 'markersize', 15);
set(gca, 'plotboxaspectratio', [1.5 1 1]); xlim([0 33]);
xlabel('Order'); ylabel('Normalized cumulative sum');

r = find(cumsum(S)/sum(S)>0.9, 1, 'first')