clear all; close all; clc;
load('./data/sampledataEEG.mat');
EEG  = ALLEEG(1); % loading continuous data

%% Independent component analysis (ICA)
X = double(EEG.data);
[Nelectrodes, Ntime, Ntrials] = size(X);

% Three ways to subtract mean values:
% onevec = ones(Ntime,1);
% X = X*(eye(Ntime)-onevec*onevec'/Ntime); % (this is memory consuming)
% X = X - repmat(mean(X,2), 1, Ntime);
X = bsxfun(@minus, X, mean(X,2));
 
eta = 0.001;
Wold = eye(Nelectrodes);
Niter = 2000;
err = 1;
 
iter = 1;
while iter<Niter && err>10e-6
    Y = Wold*X;
    phi = tanh(Y);
    Wnew = Wold + eta*(eye(Nelectrodes)-phi*Y'/Ntime)*Wold;
    err = norm(Wnew-Wold,'fro')/norm(Wold,'fro');
    Wold = Wnew;
    iter = iter+1;
end
 
Y = Wnew*X;
W = Wnew;
A = inv(W);

%%
% Fourier transform:
f = [0:(Ntime-1)]/Ntime*EEG.srate;
fY = zeros(size(Y));
for nn=1:Nelectrodes
    fY(nn,:) = fft(Y(nn,:));
end

figure(1); clf;
for nn=1:5
    subplot(5, 5, 5*(nn-1)+1);
    topoplot(A(:,nn), EEG.chanlocs);
    subplot(5, 5, 5*(nn-1)+[2 3]);
    plot(EEG.times/1000, Y(nn,:), 'k');
    xlim([0 100]);
    xlabel('Time (s)'); ylabel('EEG (\mu V)');
    subplot(5, 5, 5*(nn-1)+[4 5]);
    plot(f, 20*log10(abs(fY(nn,:))), 'k');
    xlabel('Frequency (Hz)'); ylabel('10 log (Power)');
    xlim([0 30]); ylim([-0 80])
end
set(gcf, 'color', 'w');