clear all; close all; clc;
load('./data/EmptyEEG.mat');
lf.GainN = bsxfun(@times,squeeze(lf.Gain(:,1,:)),lf.GridOrient(:,1)') ...
    + bsxfun(@times,squeeze(lf.Gain(:,2,:)),lf.GridOrient(:,2)') + ...
    bsxfun(@times,squeeze(lf.Gain(:,3,:)),lf.GridOrient(:,3)');

[Nelectrodes, Nsources] = size(lf.GainN);

%%
dip = 109;
% lf.GridLoc: 3d positions of 2004 cortical sources:
figure(1); clf;hold on;
plot3(lf.GridLoc(:,1), lf.GridLoc(:,2), lf.GridLoc(:,3),...
    'k.', 'Markersize', 10);
plot3(lf.GridLoc(dip,1), lf.GridLoc(dip,2), lf.GridLoc(dip,3),...
    'r.', 'Markersize', 30);
axis equal off;
view([-90 90]);
set(gcf, 'color', 'w')
% lf.GainN: Corresponding scalp maps:
figure(2); clf;
topoplot(lf.GainN(:,dip), EEG.chanlocs,'numcontour',0, ...
    'electrodes','off','shading','interp');
set(gcf, 'color', 'w')

%%
electrodeNames = {EEG.chanlocs.labels};
X = [EEG.chanlocs.X]/1000;
Y = [EEG.chanlocs.Y]/1000;
Z = [EEG.chanlocs.Z]/1000;
elec = 31; % Pz
gain = lf.GainN(elec, :);
gain = (gain-mean(gain))/std(gain);
gain = -3*(gain<=-3) + 3*(gain>=3) + gain.*(abs(gain)<3);
Nbins = 100;
bins = linspace(-3, 3, Nbins);
clrs = jet(Nbins);
figure(3); clf; hold on;
for nn=1:Nsources
    idx = dsearchn(bins', gain(nn));
    plot3(lf.GridLoc(nn,1), lf.GridLoc(nn,2), lf.GridLoc(nn,3), ...
        '.', 'MarkerSize', 15, ...
        'color', clrs(idx,:));
end
view([-90 90]);
axis equal off;
set(gcf, 'color', 'w');
% plot3(X(31), Y(31), Z(31), 'k.', 'MarkerSize', 30);

figure(4); clf;
% topoplot(zeros(Nelectrodes,1), EEG.chanlocs,'numcontour',0, ...
%     'electrodes','on','shading','interp');
topoplot(zeros(Nelectrodes,1), EEG.chanlocs,'numcontour',0, ...
    'electrodes','labels','shading','interp');
set(gcf, 'color', 'w')