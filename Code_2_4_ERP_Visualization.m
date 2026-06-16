% eeglab; close; % Adding EEGLAB path:
clear all; close all; clc;
load('./data/sampledataEEG.mat');
EEG  = ALLEEG(2); % loading epoched data

%% 2.4 ERP visualization:
X = EEG.data;
[Nelectrodes, Ntime, Ntrials] = size(X);
Xav = mean(X,3);
Xsd = std(X,[], 3);

% Select Fz channnel for example:
elec_name = 'Fz';
elec = find(strcmp({EEG.chanlocs.labels}, elec_name));

% ERP (trial average):
figure(1); clf; hold on;
erp_sd_x = [EEG.times EEG.times(end:(-1):1) EEG.times(1)];
erp_sd_y = [Xav(elec,:)-1*Xsd(elec,:) Xav(elec,end:(-1):1)+1*Xsd(elec,end:(-1):1) Xav(elec,1)-1*Xsd(elec,1)];
ff = fill(erp_sd_x, erp_sd_y, [0.85 0.85 0.85]); 
set(ff, 'edgecolor', 'w');
plot(EEG.times, Xav(elec,:), 'k', 'linewidth', 2);
xlabel('Time (ms)', 'fontsize', 14);
ylabel('EEG (\mu V)', 'fontsize', 14)
% title(['Electrode ' EEG.chanlocs(elec).labels ': ERP'])
xlim([-500 1000]); ylim([-30 50]);
set(gca, 'plotboxaspectratio', [2 1 1], 'xtick', [-500 0:150:750 1000]);
set(gcf, 'color', 'w');

% topomap:
t0 = [0 150 300 450 600 750];
t0idx = dsearchn(EEG.times', t0')';
figure(2); clf; hold on; 
for nn=1:length(t0)
    subplot(1, 6, nn);
    currentX = mean(X(:,t0idx(nn)-5:t0idx(nn)+5), 2);
    topoplot(currentX, EEG.chanlocs);
    title([num2str(t0(nn)) 'ms'], 'fontsize', 16)
end
set(gcf, 'color', 'w');

% ERP image:
Xelec = squeeze(X(elec, :, :));
figure(3); clf; hold on;
imagesc(EEG.times, 1:Ntrials, squeeze(X(elec, :, :))')
plot([0 0], [0 80], 'k-', 'linewidth', 2);
xlabel('Time (ms)', 'fontsize', 14)
ylabel('Trials', 'fontsize', 14)
axis([-500 1000 1 80]);
set(gca, 'plotboxaspectratio', [2 1 1], 'xtick', [-500 0:150:750 1000]);
% colormap('gray')

set(gcf, 'color', 'w');


% ERP image using pop_erpimage:
figure(4); clf; 
pop_erpimage(EEG, 1, [4],[[]],'Fz',10, 1,{},[],'' ,...
    'yerplabel','\muV','erp','off','cbar','off','topo', { [4] EEG.chanlocs EEG.chaninfo } );
axis([-500 1000 1 80]);
set(gca, 'plotboxaspectratio', [2 1 1], 'xtick', -500:500:1000);
% colormap(gray)
set(gcf, 'color', 'w');
