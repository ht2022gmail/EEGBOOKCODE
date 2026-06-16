clear all; close all; clc;
load('./data/sampledataEEG.mat');
EEG  = ALLEEG(1); % loading continuous data

%% Time-domain and frequency-domain convolution
X = double(EEG.data);
[Nelectrodes, Ntime, Ntrials] = size(X);
% X = X';

Nblock = EEG.srate*1.5;
time = linspace(-0.5, 1, Nblock);

rt_idx = round([ EEG.event(find(strcmp({EEG.event.type}, 'square'))).latency ]);
D = zeros(Nblock, Ntime);
for kk=1:length(rt_idx)
    idx = rt_idx(kk)-EEG.srate/2;
    D(:,idx:idx+Nblock-1) = eye(Nblock);
end

A = X*pinv(D);
figure(1); 
plot(time, A);
xlim([-0.25 1])

%
CS = (A*D)*(A*D)'/Ntime;
CR = X*X'/Ntime;

[U, Lambda] = eig(CS, CR);
Lambda = diag(Lambda);
[Lambda, idx] = sort(Lambda, 'descend');
U = U(:,idx);

Y = U'*A;
figure(2); 
plot(time, Y(1:3,:))