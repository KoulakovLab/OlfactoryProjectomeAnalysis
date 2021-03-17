% Auxillary stuff that I want to look at;
%   * OB injection general strength for pattern alignment
%   * OB injection IPR histogram
% ! * Conprob for OB vs regions
figaux = gobjects(3, 1);

% AUX.A) Projection strength OB
% Check general level of projections accross regions
figaux(1) = figure;
set(figaux(1), 'name', 'AUX.A: Projection patterns');
sto = aliOB.getDataset(1:(aliOB.nData-1));
hold('on');
for i = 1:length(sto)
    semilogy(imgaussfilt(mean(sto(i).prjImg, 1), 1), 'LineWidth', 3);
    hold('on');
end
hold('off');
legend({sto.name});
xlabel('Source slice');
ylabel('Mean projection');
title('Brain based slice means');

% AUX.B) Projection strength OB
% Check general level of projections accross regions
figaux(2) = figure;
set(figaux(2), 'name', 'AUX.B: IPR histogram');
histogram(obdata.brcPrjIpr);
xlabel('Inverse Participation Ratio (IPR)');
ylabel('Frequency (barcodes)');
title('IPR Histogram');

% AUX.C) Histograms of different slices in different datasets of OB data