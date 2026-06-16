clear all; close all; clc;

BCIsummary = readtable('Summary_BCI_Performance_Maslova2023_mod.csv', 'Delimiter', 'comma');

BCIparadigm = table2cell(BCIsummary(:,1));
authors = table2cell(BCIsummary(:,2));
bitrate = table2array(BCIsummary(:,3));
accuracy = table2array(BCIsummary(:,4));

figure(1); clf;
b_bitrate = bar(authors, bitrate);
b_bitrate.FaceColor = 'flat';
for nn=1:4
    b_bitrate.CData(nn,:) = [1.0 0.3 0.3];
end
for nn=5:6
    b_bitrate.CData(nn,:) = [0.3 1.0 0.3];
end
for nn=7:12
    b_bitrate.CData(nn,:) = [0.3 0.3 1.0];
end
for nn=13:15
    b_bitrate.CData(nn,:) = [1.0 1.0 0.3];
end
ylabel('bit rate (bit/min)');
set(gca, 'PlotBoxAspectRatio', [2 1 1]);
set(gcf, 'color', 'w');

figure(2); clf;
b_accuracy = bar(authors, accuracy);
b_accuracy.FaceColor = 'flat';
for nn=1:4
    b_accuracy.CData(nn,:) = [1.0 0.3 0.3];
end
for nn=5:6
    b_accuracy.CData(nn,:) = [0.3 1.0 0.3];
end
for nn=7:12
    b_accuracy.CData(nn,:) = [0.3 0.3 1.0];
end
for nn=13:15
    b_accuracy.CData(nn,:) = [1.0 1.0 0.3];
end
% set(b_accuracy, 'FaceColor', [0.6 0.6 0.6])
ylabel('bit rate (bit/min)');
set(gca, 'PlotBoxAspectRatio', [2 1 1]);
set(gcf, 'color', 'w');

%%
figure(1); clf; hold on;
b_bitrate = bar(authors(1:12), bitrate(1:12));
b_bitrate.FaceColor = 'flat';
for nn=1:4
    b_bitrate.CData(nn,:) = [1.0 0.3 0.3];
end
for nn=5:6
    b_bitrate.CData(nn,:) = [0.3 1.0 0.3];
end
for nn=7:12
    b_bitrate.CData(nn,:) = [0.3 0.3 1.0];
end
plot([1.5 12.5], [1000 1000], 'r--', 'linewidth', 2)
ylabel('Bit rate (bit/min)');
set(gca, 'PlotBoxAspectRatio', [2 1 1]);
set(gcf, 'color', 'w');
% xlim([0.5 12.5]); 
ylim([0 1200]);
grid on;

figure(2); clf; hold on;
b_accuracy = bar(authors(1:12), accuracy(1:12));
b_accuracy.FaceColor = 'flat';
for nn=1:4
    b_accuracy.CData(nn,:) = [1.0 0.3 0.3];
end
for nn=5:6
    b_accuracy.CData(nn,:) = [0.3 1.0 0.3];
end
for nn=7:12
    b_accuracy.CData(nn,:) = [0.3 0.3 1.0];
end
plot([0.5 12.5], [100 100], 'r--', 'linewidth', 2)
ylabel('Accuracy (%)');
set(gca, 'PlotBoxAspectRatio', [2 1 1]);
set(gcf, 'color', 'w');
grid on;
ylim([80 110])