% Clustering analysis
load('data/classifier.mat');
load('data/processed.mat');

ELIM_AMNT = .85;
COLS = lines(3);
% Get the corner points
[~, ~, CORNERS] = aux.genSimplex(2);
% Calculate the edge points opposite this corner
EDGES = zeros(3, 2);
for c = 1:3
  c1 = mod(c + 1 - 1, 3) + 1;
  c2 = mod(c + 2 - 1, 3) + 1;
  EDGES(c, :) = aux.getPerpToPts(...
    obdata.data.classifier.pts([c1, c2], :), ...
    CORNERS([c1, c2], :));
end
figure6 = gobjects(25, 1);

%% Overview
figure6(1) = figure;
set(figure6(1), 'name', 'Overview');

sto_ax(2) = subplot(1, 2, 1);
hold('on');
for c = 1:obdata.data.temp.nIds
  sto_ids = obdata.data.temp.realId == c;
  scatter( ...
    obdata.brcCor(sto_ids, 1), ...
    obdata.brcCor(sto_ids, 2), ...
    100, COLS(c, :), 'filled', 'MarkerEdgeColor', 'k');
end
hold('off');
xticks([]);
yticks([]);
legend(obdata.data.temp.names);
title('Template barcodes');
axis('image');

sto_ax(1) = subplot(1, 2, 2);
colormap(sto_ax(1), aux.viridis);
scatter(obdata.brcCor(:, 1), obdata.brcCor(:, 2), ...
  35, obdata.brcPrjIpr, 'filled', 'MarkerEdgeColor', 'k');
sto_bar = colorbar;
title(sto_bar, 'IPR');
xticks([]);
yticks([]);
title('All barcodes');
axis('image');
linkaxes(sto_ax);

% Different brains image
figure6(2) = figure;
set(figure6(2), 'name', 'Brains');
COL_BR = aux.pmkmp(length(obdata.brcName{1}), 'IsoL');
hold('on');
for i = 1:length(obdata.brcName{1})
  sto_ind = obdata.brcId(:, 1) == i;
  scatter(obdata.brcCor(sto_ind, 1), obdata.brcCor(sto_ind, 2), ...
  20, COL_BR(i, :), 'filled', 'MarkerEdgeColor', 'k');
end
hold('off');
legend(obdata.brcName{1});
xticks([]);
yticks([]);
title('All barcodes; different brains.');
axis('image');
sgtitle('Barcodes in tSNE space');

savefig(figure6(1:2), 'data/figures/classifier_overview.fig');

%% Neural network

% Neural Network; plotting
figure6(3) = figure;
set(figure6(3), 'name', 'Classifier on Templates');

% TSNE space plots for the template data
subplot(2, 3, 1);
hold('on')
sto_leg = cell(obdata.data.temp.nIds);
sto_plt = gobjects(obdata.data.temp.nIds);
% Draw circles with the classifier results
for cR = 1:obdata.data.temp.nIds
  for cC = 1:obdata.data.temp.nIds
    sto_ids = ...
      (obdata.data.temp.classId == cC) & (obdata.data.temp.realId == cR);
    if sum(sto_ids) == 0
      data = nan(1, 2);
    else
      data = obdata.data.temp.tsneCor(sto_ids, :);
    end
    sto_plt(cR, cC) = scatter(data(:, 1), data(:, 2), ...
      30, 'filled', 'LineWidth', 1, ...
      'MarkerFaceColor', COLS(cR, :), 'MarkerEdgeColor', COLS(cC, :));
    sto_leg{cR, cC} = [ ...
      'True: ',  obdata.data.temp.names{cR}, ' ', ...
      'Class: ', obdata.data.temp.names{cC}];
  end
end
hold('off');
legend(sto_plt(:), sto_leg(:), ...
  'Location', 'southoutside', 'NumColumns', obdata.data.temp.nIds);
xticks([]);
yticks([]);
title('tSNE Space');

% Template data in probability space
subplot(2, 3, 4);
hold('on');
sto_leg = cell(3, 2);
% Draw decision boundaries corresponding to each class
for c = 1:obdata.data.temp.nIds
  % Draw this diamong
  sto_cor = [ ...
    CORNERS(c, :); ...
    EDGES(mod(c + 1 - 1, 3) + 1, :); ...
    obdata.data.classifier.center; ...
    EDGES(mod(c + 2 - 1, 3) + 1, :)];
  fill(sto_cor(:, 1), sto_cor(:, 2), COLS(c, :), 'FaceAlpha', .25);
  sto_leg{c, 1} = ['Decision region: ', obdata.data.temp.names{c}];
end
% Put points in this probability space
for c = 1:obdata.data.temp.nIds
  sto_ids = obdata.data.temp.realId == c;
  sto_cor = obdata.data.temp.conv_fromBaryo( ...
    obdata.data.temp.prob(sto_ids, :));
  scatter(sto_cor(:, 1), sto_cor(:, 2), ...
    15, COLS(c, :), 'filled');
  sto_leg{c, 2} = ['Template, ', obdata.data.temp.names{c}];
end
% Show the control points
% for c = 1:obdata.data.temp.nIds
%   scatter( ...
%     obdata.data.classifier.pts(c, 1), obdata.data.classifier.pts(c, 2), ...
%     150, COLS(c, :), 'diamond', 'filled', ...
%     'MarkerEdgeColor', .5 * COLS(c, :), 'MarkerFaceAlpha', .5);
%   sto_leg{c, 3} = [obdata.data.temp.names{c}, ' control pt.'];
% end
hold('off');
legend(sto_leg(:), 'Location', 'southoutside', 'NumColumns', ...
  size(sto_leg, 2));
xticks([]);
yticks([]);
title('Probability Simplex');
axis('image');

% ROC curve plot for individual classes
subplot(2, 3, 2);
hold('on');
sto_leg = cell(3, 1);
for c = 1:obdata.data.temp.nIds
  plot(obdata.data.temp.fpr(:, c), obdata.data.temp.tpr(:, c), ...
    'Color', COLS(c, :), 'LineWidth', 3);
  sto_leg{c} = [obdata.data.temp.names{c}, ' (AUC: ', ...
    num2str(obdata.data.temp.auc(c), 4), ')'];
end
legend(sto_leg, 'Location', 'southeast', 'AutoUpdate', 'off');
plot([0, 1], [0, 1], ':', 'LineWidth', 1, 'Color', 'k');
hold('off');
xlabel('False positive rate');
ylabel('True positive rate');
title('ROC Curve for bin. class.');
axis('image');

% Confusion matrix
subplot(2, 3, 5);
confusionchart(obdata.data.temp.confusionMatrix, obdata.data.temp.names);
title('Confusion matrix');

% Template loadings
sto_spl = subplot(2, 3, 3);
% Create a loading image ranking;
% - First rank ones that are mitral in mitral score
% - Then rank ones that are tufted in tufted score
sto_chosenScore = obdata.data.temp.prob(sub2ind( ...
  [obdata.data.temp.nBrc, obdata.data.temp.nIds], ...
  (1:obdata.data.temp.nBrc)', obdata.data.temp.classId));
[~, sto_ord] = sortrows([obdata.data.temp.classId, (1 - sto_chosenScore)]);
% Find seperation indices
nClass = sum(obdata.data.temp.classId == 1:3, 1);
% Get image, normalize and sort
% sto_img = obdata.data.temp.prjReg ./ sum(obdata.data.temp.prjReg, 2);
sto_img = obdata.data.temp.prjReg ./ max(obdata.data.temp.prjReg, [], 2);
sto_img = sto_img(sto_ord, :);
sto_ids = obdata.data.temp.realId(sto_ord);
sto_cim = zeros([size(sto_img), 3]);
% Create image with white background, and real ID foreground
for c = 1:obdata.data.temp.nIds
  sto_this = (sto_ids == c);
  % Convert the matrix into a value index matrix; 1-256
  sto_vin = floor(sto_img(sto_this, :) * 255) + 1;
  % Create a colormap from white to COLMAP color
  for r = 1:3
    sto_cmap = linspace(1, COLS(c, r), 256);
    sto_cim(sto_this, :, r) = sto_cmap(sto_vin);
  end
end
image(sto_cim);
% Draw seperation regions to split different classes
hold('on');
for c = 2:obdata.data.temp.nIds
  plot( ...
    .5 + [0, obdata.nPrjReg], ...
    ones(1, 2) .* (.5 + sum(nClass(1:(c - 1)))), ...
    'k', 'LineWidth', 3.5);
end
hold('off');
xlabel('Regions');
xticks(1:obdata.nPrjReg);
xticklabels(obdata.prjRegName);
ylabel('Templates');
yticks(cumsum([0, nClass(1:(end - 1))]) + .5 * nClass);
yticklabels(cellfun(@(x) ['Class. ', x], obdata.data.temp.names, ...
  'UniformOutput', false));
ytickangle(90);
title('Template Loadings per Class');

% x-y location of templates
subplot(2, 3, 6);
hold('on');
for c = 1:obdata.data.temp.nIds
  sto_this = (sto_ids == c);
  scatter( ...
    obdata.data.temp.location(sto_this, 1), ...
    obdata.data.temp.location(sto_this, 2), ...
    10, COLS(c, :), 'filled');
end
hold('off');
legend(obdata.data.temp.names);
ylabel('Mitral layer (- to inside)');
xlabel('Glomeruli layer (+ to inside)');
title('Distance of cells to layers (um)');

sgtitle('Template Classification');

savefig(figure6(3), 'data/figures/classifier_template.fig');

%% Results of real classification
figure6(4) = figure;
set(figure6(4), 'name', 'Classifier on Barcodes');

% tSNE space plot for barcodes
subplot(2, 2, 1);
hold('on');
for c = 1:obdata.data.temp.nIds
  sto_ids = obdata.brcId(:, 2) == c;
  scatter( ...
    obdata.brcCor(sto_ids, 1), ...
    obdata.brcCor(sto_ids, 2), ...
    15, COLS(c, :), 'filled');
end
hold('off');
legend(cellfun(@(x) [x , ' Barcodes'], aliOB.data.temp.names, ...
  'UniformOutput', false), 'Location', 'southoutside');
xticks([]);
yticks([]);
title('tSNE: Barcodes');

% tSNE space for templates
subplot(2, 2, 2);
hold('on');
% Templates
for c = 1:obdata.data.temp.nIds
  sto_ids = obdata.data.temp.realId == c;
  scatter( ...
    obdata.data.temp.tsneCor(sto_ids, 1), ...
    obdata.data.temp.tsneCor(sto_ids, 2), ...
    15, COLS(c, :), 'filled');
end
hold('off');
legend(cellfun(@(x) [x , ' Templates'], aliOB.data.temp.names, ...
  'UniformOutput', false), 'Location', 'southoutside');
xticks([]);
yticks([]);
title('tSNE: Barcodes');

% Bar plot of frequencies
subplot(2, 2, 3);
% Get number of barcodes; per brain; in each class
sto_barplot = zeros(obdata.data.temp.nIds, length(obdata.brcName{1}));
for b = 1:length(obdata.brcName{1})
  % Get barcodes within this brain
  sto_ind = obdata.brcId(:, 1) == b;
  % Get the counts of types ID'd in this, normalized by total count
  sto_barplot(:, b) = sum( ...
    obdata.brcId(sto_ind, 2) == (1:obdata.data.temp.nIds), 1)' ./ ...
    sum(sto_ind);
end
bar(sto_barplot);
xticklabels(obdata.data.temp.names);
legend(obdata.brcName{1});
title('Barcode count freq.');
ylabel('Perc. of barcodes');

% Barcodes In probability simplex
subplot(2, 2, 4);
hold('on');
sto_leg = cell(3, 2);
sto_prb = aliOB.data.brc_nnetProb{end};
% Draw decision boundaries corresponding to each class
for c = 1:obdata.data.temp.nIds
  % Draw this diamond
  sto_cor = [ ...
    CORNERS(c, :); ...
    EDGES(mod(c + 1 - 1, 3) + 1, :); ...
    obdata.data.classifier.center; ...
    EDGES(mod(c + 2 - 1, 3) + 1, :)];
  fill(sto_cor(:, 1), sto_cor(:, 2), COLS(c, :), 'FaceAlpha', .25);
  sto_leg{c, 1} = ['Decision region: ', obdata.data.temp.names{c}];
end
% Put points in this probability space
for c = 1:obdata.data.temp.nIds
  sto_ids = aliOB.brcId{end}(:, 2) == c;
  sto_cor = aliOB.data.temp.conv_fromBaryo(sto_prb(sto_ids, :));
  scatter(sto_cor(:, 1), sto_cor(:, 2), ...
    15, COLS(c, :), 'filled');
  sto_leg{c, 2} = ['Barcodes, ', obdata.data.temp.names{c}];
end
hold('off');
legend(sto_leg(:), 'Location', 'southoutside', 'NumColumns', 2);
xticks([]);
yticks([]);
title('Probability Simplex');
axis('image');

sgtitle('Classified Barcodes');

savefig(figure6(4), 'data/figures/classifier_performance.fig');


%% Result distribution of networks
% FIGURE 6D - Loadings
figure6(5) = figure;
set(figure6(5), 'name', 'Loadings');
% Get barcode probabilities
sto_prb = obdata.data.classifier.net(obdata.prjRegSum')';
sto_sco = sto_prb(sub2ind([obdata.nBrc, obdata.data.temp.nIds], ...
  (1:obdata.nBrc)', obdata.brcId(:, 2)));
[~, sto_ord] = sortrows([obdata.brcId(:, 2), (1 - sto_sco)]);
% Find seperation indices
nClass = sum(obdata.brcId(:, 2) == 1:3, 1);
% Get image, normalize and sort
sto_img = obdata.prjRegSum ./ max(obdata.prjRegSum, [], 2);
sto_img = sto_img(sto_ord, :);
sto_ids = obdata.brcId(sto_ord, 2);
% Create image with white background, and real ID foreground
for c = 1:obdata.data.temp.nIds
  % Get the image for this class
  sto_ind = obdata.brcId(:, 2) == c;
  sto_img = obdata.prjRegSum(sto_ind, :);
  sto_prb = obdata.data.classifier.net(sto_img')';
  % Sort according to classifier scores
  [sto_sco, sto_ord] = sort(sto_prb(:, c), 'descend');
  sto_img = sto_img(sto_ord, :);
  sto_img = sto_img ./ max(sto_img, [], 2);
  % Convert the matrix into a value index matrix; 1-256
  sto_vin = floor(sto_img	* 255) + 1;
  % Create a colormap from white to COLMAP color
  sto_cim = zeros([size(sto_vin), 3]);
  for r = 1:3
    sto_cmap = linspace(1, COLS(c, r), 256);
    sto_cim(:, :, r) = sto_cmap(sto_vin);
  end
  subplot(1, obdata.data.temp.nIds, c);
  image(sto_cim);
  hold('on');
  sto_x = xlim;
  % Draw cutoff line
  sto_cut = find(sto_sco < ELIM_AMNT, 1);
  plot(sto_x, ones(1, 2) .* (.5 + sto_cut), ':k', 'LineWidth', 3);
  hold('off');
  xlabel('Regions');
  xticks(1:obdata.nPrjReg);
  xticklabels(obdata.prjRegName);
  ylabel('Barcodes');
  title(['Putative ', obdata.data.temp.names{c}, ' Cells']);
  set(gca, 'YDir','reverse');
  set(gca, 'FontSize', 20);
end
savefig(figure6(5), 'data/figures/classifier_loadings.fig');

%% Network results
figure6(6) = figure;
set(figure6(6), 'name', 'Network Results');
s_fpr = zeros(size(nn_sets_test));
s_tpr = zeros(size(nn_sets_test));
for t = 1:3
  subplot(2, 6, (1:2) + 2*(t-1));
  aux.nnmatPlotNEW(obdata_mitral, t);
  s_fpr(:, t) = 1 - ([nn_sets_test(:, t).mitral_test_tnr]');
  s_tpr(:, t) = 0 + ([nn_sets_test(:, t).mitral_test_tpr]');
end
subplot(2, 6, 7:9);
s_plot = @(x) aux.distributionPlot(x, 'histOpt', 1, 'divFactor', 1, 'showMM', false);
s_plot(s_fpr);
title('False Positive Rate');
xticklabels({obdata_mitral.data.typenet.name});
xlabel('Different Networks');
ylim([0, 1]);
subplot(2, 6, 10:12);
s_plot(s_tpr);
title('True Positive Rate');
xticklabels({obdata_mitral.data.typenet.name});
xlabel('Different Networks');
ylim([0, 1]);

% Network test results
ARCNAME = {obdata_mitral.data.typenet(:).name};
figure6(7) = figure;
set(figure6(7), 'name', 'Confusion (all arch.)');
subplot(2, 5, 1);
sto_conf = round(mean(cat(3, nn_sets(:).mitral_train_conf), 3));
confusionchart(sto_conf, {'Negative', 'Positive'});
title('100%');
subplot(2, 5, 6);
sto_conf = round(mean(cat(3, nn_sets(:).train_conf), 3));
confusionchart(sto_conf, obdata.data.temp.names);
title('100%');
subplot(2, 5, 2);
sto_conf = round(mean(cat(3, nn_sets_test(:).mitral_train_conf), 3));
confusionchart(sto_conf, {'Negative', 'Positive'});
title('90%');
subplot(2, 5, 7);
sto_conf = round(mean(cat(3, nn_sets_test(:).train_conf), 3));
confusionchart(sto_conf, obdata.data.temp.names);
title('90%');
subplot(2, 5, 3);
sto_conf = round(mean(cat(3, nn_sets_test(:).mitral_test_conf), 3));
confusionchart(sto_conf, {'Negative', 'Positive'});
title('10%');
subplot(2, 5, 8);
sto_conf = round(mean(cat(3, nn_sets_test(:).test_conf), 3));
confusionchart(sto_conf, obdata.data.temp.names);
title('10%');
subplot(2, 5, 4);
sto_conf = round(mean(cat(3, resper2(:).mitral_train_conf), 3));
confusionchart(sto_conf, {'Negative', 'Positive'});
title('80%');
subplot(2, 5, 9);
sto_conf = round(mean(cat(3, resper2(:).train_conf), 3));
confusionchart(sto_conf, obdata.data.temp.names);
title('80%');
subplot(2, 5, 5);
sto_conf = round(mean(cat(3, resper2(:).mitral_test_conf), 3));
confusionchart(sto_conf, {'Negative', 'Positive'});
title('20%');
subplot(2, 5, 10);
sto_conf = round(mean(cat(3, resper2(:).test_conf), 3));
confusionchart(sto_conf, obdata.data.temp.names);
title('20%');
sgtitle('Comparison of network performances; all architectures');
for c = 1:3
  figure6(7 + c) = figure;
  set(figure6(7 + c), 'name', ...
    ['Confusion (', ARCNAME{c}, ')']);
  subplot(2, 5, 1);
  sto_conf = round(mean(cat(3, nn_sets(:, c).mitral_train_conf), 3));
  confusionchart(sto_conf, {'Negative', 'Positive'});
  title('100%');
  subplot(2, 5, 6);
  sto_conf = round(mean(cat(3, nn_sets(:, c).train_conf), 3));
  confusionchart(sto_conf, obdata.data.temp.names);
  title('100%');
  subplot(2, 5, 2);
  sto_conf = round(mean(cat(3, nn_sets_test(:, c).mitral_train_conf), 3));
  confusionchart(sto_conf, {'Negative', 'Positive'});
  title('90%');
  subplot(2, 5, 7);
  sto_conf = round(mean(cat(3, nn_sets_test(:, c).train_conf), 3));
  confusionchart(sto_conf, obdata.data.temp.names);
  title('90%');
  subplot(2, 5, 3);
  sto_conf = round(mean(cat(3, nn_sets_test(:, c).mitral_test_conf), 3));
  confusionchart(sto_conf, {'Negative', 'Positive'});
  title('10%');
  subplot(2, 5, 8);
  sto_conf = round(mean(cat(3, nn_sets_test(:, c).test_conf), 3));
  confusionchart(sto_conf, obdata.data.temp.names);
  title('10%');
  subplot(2, 5, 4);
  sto_conf = round(mean(cat(3, resper2(:, c).mitral_train_conf), 3));
  confusionchart(sto_conf, {'Negative', 'Positive'});
  title('80%');
  subplot(2, 5, 9);
  sto_conf = round(mean(cat(3, resper2(:, c).train_conf), 3));
  confusionchart(sto_conf, obdata.data.temp.names);
  title('80%');
  subplot(2, 5, 5);
  sto_conf = round(mean(cat(3, resper2(:, c).mitral_test_conf), 3));
  confusionchart(sto_conf, {'Negative', 'Positive'});
  title('20%');
  subplot(2, 5, 10);
  sto_conf = round(mean(cat(3, resper2(:, c).test_conf), 3));
  confusionchart(sto_conf, obdata.data.temp.names);
  title('20%');
  sgtitle(['Comparison of network performances; ', ARCNAME{c}]);
end

% The comparison with the chosen network
CC = {1:3, 1, 2, 3};
CCT = {'All arch', obdata_mitral.data.typenet(:).name};
for c = 1:length(CC)
  C = CC{c};
  T = CCT{c};
  figure6(11 + 3 * (c - 1) + 0) = figure;
  set(figure6(11 + 3 * (c - 1) + 0), 'name', ['Comparison ', T, ' (thr.)']);
  subplot(2, 3, 1);
  histogram([nn_sets(:, C).mapseq_mitral_tpr]);
  hold('on');
  histogram([nn_sets(:, C).mapseq_mitral_fpr]);
  hold('off');
  xlim([0, 1]);
  xlabel('Rate');
  ylabel('Count');
  avgRho = mean([nn_sets(:, C).mapseq_mitral_corr]);
  title(['Training 100%, Mapseq; \rho=', num2str(avgRho)]);
  legend({'TPR', 'FPR'}, 'Location', 'north');
  subplot(2, 3, 2);
  histogram([nn_sets_test(:, C).mapseq_mitral_tpr]);
  hold('on');
  histogram([nn_sets_test(:, C).mapseq_mitral_fpr]);
  hold('off');
  xlim([0, 1]);
  xlabel('Rate');
  ylabel('Count');
  avgRho = mean([nn_sets_test(:, C).mapseq_mitral_corr]);
  title(['Training 90%, Mapseq; \rho=', num2str(avgRho)]);
  legend({'TPR', 'FPR'}, 'Location', 'north');
  subplot(2, 3, 3);
  histogram([resper2(:, C).mapseq_mitral_tpr]);
  hold('on');
  histogram([resper2(:, C).mapseq_mitral_fpr]);
  hold('off');
  xlim([0, 1]);
  xlabel('Rate');
  ylabel('Count');
  avgRho = mean([resper2(:, C).mapseq_mitral_corr]);
  title(['Training 80%, Mapseq; \rho=', num2str(avgRho)]);
  legend({'TPR', 'FPR'}, 'Location', 'north');
  subplot(2, 3, 4);
  histogram([nn_sets(:, C).barseq_mitral_tpr]);
  hold('on');
  histogram([nn_sets(:, C).barseq_mitral_fpr]);
  hold('off');
  xlim([0, 1]);
  xlabel('Rate');
  ylabel('Count');
  avgRho = mean([nn_sets(:, C).barseq_mitral_corr]);
  title(['Training 100%, Barseq; \rho=', num2str(avgRho)]);
  legend({'TPR', 'FPR'}, 'Location', 'north');
  subplot(2, 3, 5);
  histogram([nn_sets_test(:, C).barseq_mitral_tpr]);
  hold('on');
  histogram([nn_sets_test(:, C).barseq_mitral_fpr]);
  hold('off');
  xlim([0, 1]);
  xlabel('Rate');
  ylabel('Count');
  avgRho = mean([nn_sets_test(:, C).barseq_mitral_corr]);
  title(['Training 90%, Barseq; \rho=', num2str(avgRho)]);
  legend({'TPR', 'FPR'}, 'Location', 'north');
  subplot(2, 3, 6);
  histogram([resper2(:, C).barseq_mitral_tpr]);
  hold('on');
  histogram([resper2(:, C).barseq_mitral_fpr]);
  hold('off');
  xlim([0, 1]);
  xlabel('Rate');
  ylabel('Count');
  avgRho = mean([resper2(:, C).barseq_mitral_corr]);
  title(['Training 80%, Barseq; \rho=', num2str(avgRho)]);
  legend({'TPR', 'FPR'}, 'Location', 'north');
  sgtitle('Comparison wrt the chosen network; thresholding');
  
  figure6(11 + 3 * (c - 1) + 1) = figure;
  set(figure6(11 + 3 * (c - 1) + 1), 'name', ['Comparison, ', T, ' (max)']);
  subplot(2, 3, 1);
  histogram(arrayfun(@(x) x.mapseq_tpr(1), nn_sets(:, C)), linspace(0, 1, 20));
  hold('on');
  histogram(arrayfun(@(x) x.mapseq_tpr(2), nn_sets(:, C)), linspace(0, 1, 20));
  histogram(arrayfun(@(x) x.mapseq_tpr(3), nn_sets(:, C)), linspace(0, 1, 20));
  histogram(arrayfun(@(x) x.mapseq_fpr(1), nn_sets(:, C)), linspace(0, 1, 20));
  histogram(arrayfun(@(x) x.mapseq_fpr(2), nn_sets(:, C)), linspace(0, 1, 20));
  histogram(arrayfun(@(x) x.mapseq_fpr(3), nn_sets(:, C)), linspace(0, 1, 20));
  hold('off');
  xlim([0, 1]);
  xlabel('Rate');
  ylabel('Count');
  title('Training 100%, Mapseq');
  legend({'TPR pMC', 'TPR pTC', 'TPR pDC', 'FPR pMC', 'FPR pTC', 'FPR pDC'}, 'Location', 'north');
  subplot(2, 3, 2);
  histogram(arrayfun(@(x) x.mapseq_tpr(1), nn_sets_test(:, C)), linspace(0, 1, 20));
  hold('on');
  histogram(arrayfun(@(x) x.mapseq_tpr(2), nn_sets_test(:, C)), linspace(0, 1, 20));
  histogram(arrayfun(@(x) x.mapseq_tpr(3), nn_sets_test(:, C)), linspace(0, 1, 20));
  histogram(arrayfun(@(x) x.mapseq_fpr(1), nn_sets_test(:, C)), linspace(0, 1, 20));
  histogram(arrayfun(@(x) x.mapseq_fpr(2), nn_sets_test(:, C)), linspace(0, 1, 20));
  histogram(arrayfun(@(x) x.mapseq_fpr(3), nn_sets_test(:, C)), linspace(0, 1, 20));
  hold('off');
  xlim([0, 1]);
  xlabel('Rate');
  ylabel('Count');
  title('Training 90%, Mapseq');
  legend({'TPR pMC', 'TPR pTC', 'TPR pDC', 'FPR pMC', 'FPR pTC', 'FPR pDC'}, 'Location', 'north');
  subplot(2, 3, 3);
  histogram(arrayfun(@(x) x.mapseq_tpr(1), resper2(:, C)), linspace(0, 1, 20));
  hold('on');
  histogram(arrayfun(@(x) x.mapseq_tpr(2), resper2(:, C)), linspace(0, 1, 20));
  histogram(arrayfun(@(x) x.mapseq_tpr(3), resper2(:, C)), linspace(0, 1, 20));
  histogram(arrayfun(@(x) x.mapseq_fpr(1), resper2(:, C)), linspace(0, 1, 20));
  histogram(arrayfun(@(x) x.mapseq_fpr(2), resper2(:, C)), linspace(0, 1, 20));
  histogram(arrayfun(@(x) x.mapseq_fpr(3), resper2(:, C)), linspace(0, 1, 20));
  hold('off');
  xlim([0, 1]);
  xlabel('Rate');
  ylabel('Count');
  title('Training 80%, Mapseq');
  legend({'TPR pMC', 'TPR pTC', 'TPR pDC', 'FPR pMC', 'FPR pTC', 'FPR pDC'}, 'Location', 'north');
  subplot(2, 3, 4);
  histogram(arrayfun(@(x) x.barseq_tpr(1), nn_sets(:, C)), linspace(0, 1, 20));
  hold('on');
  histogram(arrayfun(@(x) x.barseq_tpr(2), nn_sets(:, C)), linspace(0, 1, 20));
  histogram(arrayfun(@(x) x.barseq_tpr(3), nn_sets(:, C)), linspace(0, 1, 20));
  histogram(arrayfun(@(x) x.barseq_fpr(1), nn_sets(:, C)), linspace(0, 1, 20));
  histogram(arrayfun(@(x) x.barseq_fpr(2), nn_sets(:, C)), linspace(0, 1, 20));
  histogram(arrayfun(@(x) x.barseq_fpr(3), nn_sets(:, C)), linspace(0, 1, 20));
  hold('off');
  xlim([0, 1]);
  xlabel('Rate');
  ylabel('Count');
  title('Training 100%, Barseq');
  legend({'TPR pMC', 'TPR pTC', 'TPR pDC', 'FPR pMC', 'FPR pTC', 'FPR pDC'}, 'Location', 'north');
  subplot(2, 3, 5);
  histogram(arrayfun(@(x) x.barseq_tpr(1), nn_sets_test(:, C)), linspace(0, 1, 20));
  hold('on');
  histogram(arrayfun(@(x) x.barseq_tpr(2), nn_sets_test(:, C)), linspace(0, 1, 20));
  histogram(arrayfun(@(x) x.barseq_tpr(3), nn_sets_test(:, C)), linspace(0, 1, 20));
  histogram(arrayfun(@(x) x.barseq_fpr(1), nn_sets_test(:, C)), linspace(0, 1, 20));
  histogram(arrayfun(@(x) x.barseq_fpr(2), nn_sets_test(:, C)), linspace(0, 1, 20));
  histogram(arrayfun(@(x) x.barseq_fpr(3), nn_sets_test(:, C)), linspace(0, 1, 20));
  hold('off');
  xlim([0, 1]);
  xlabel('Rate');
  ylabel('Count');
  title('Training 90%, Barseq');
  legend({'TPR pMC', 'TPR pTC', 'TPR pDC', 'FPR pMC', 'FPR pTC', 'FPR pDC'}, 'Location', 'north');
  subplot(2, 3, 6);
  histogram(arrayfun(@(x) x.barseq_tpr(1), resper2(:, C)), linspace(0, 1, 20));
  hold('on');
  histogram(arrayfun(@(x) x.barseq_tpr(2), resper2(:, C)), linspace(0, 1, 20));
  histogram(arrayfun(@(x) x.barseq_tpr(3), resper2(:, C)), linspace(0, 1, 20));
  histogram(arrayfun(@(x) x.barseq_fpr(1), resper2(:, C)), linspace(0, 1, 20));
  histogram(arrayfun(@(x) x.barseq_fpr(2), resper2(:, C)), linspace(0, 1, 20));
  histogram(arrayfun(@(x) x.barseq_fpr(3), resper2(:, C)), linspace(0, 1, 20));
  hold('off');
  xlim([0, 1]);
  xlabel('Rate');
  ylabel('Count');
  title('Training 80%, Barseq');
  legend({'TPR pMC', 'TPR pTC', 'TPR pDC', 'FPR pMC', 'FPR pTC', 'FPR pDC'}, 'Location', 'north');
  sgtitle('Comparison wrt the chosen network; max probability');
  
  figure6(11 + 3 * (c - 1) + 2) = figure;
  set(figure6(11 + 3 * (c - 1) + 2), 'name', ['Comparison, ', T, ' (max, corr)']);
  subplot(2, 3, 1);
  imagesc(mean(cat(3, nn_sets(:, C).mapseq_corr), 3, 'omitnan'));
  colorbar;
  axis('image');
  xlabel('Chosen network label');
  ylabel('Batch network label');
  xticks(1:3);
  yticks(1:3);
  xticklabels(obdata.data.temp.names);
  yticklabels(obdata.data.temp.names);
  title('Training 100%, Mapseq');
  subplot(2, 3, 4);
  imagesc(mean(cat(3, nn_sets(:, C).barseq_corr), 3, 'omitnan'));
  colorbar;
  axis('image');
  xlabel('Chosen network label');
  ylabel('Batch network label');
  xticks(1:3);
  yticks(1:3);
  xticklabels(obdata.data.temp.names);
  yticklabels(obdata.data.temp.names);
  title('Training 100%, Barseq');
  subplot(2, 3, 2);
  imagesc(mean(cat(3, nn_sets_test(:, C).mapseq_corr), 3, 'omitnan'));
  colorbar;
  axis('image');
  xlabel('Chosen network label');
  ylabel('Batch network label');
  xticks(1:3);
  yticks(1:3);
  xticklabels(obdata.data.temp.names);
  yticklabels(obdata.data.temp.names);
  title('Training 90%, Mapseq');
  subplot(2, 3, 5);
  imagesc(mean(cat(3, nn_sets_test(:, C).barseq_corr), 3, 'omitnan'));
  colorbar;
  axis('image');
  xlabel('Chosen network label');
  ylabel('Batch network label');
  xticks(1:3);
  yticks(1:3);
  xticklabels(obdata.data.temp.names);
  yticklabels(obdata.data.temp.names);
  title('Training 90%, Barseq');
  subplot(2, 3, 3);
  imagesc(mean(cat(3, resper2(:, C).mapseq_corr), 3, 'omitnan'));
  colorbar;
  axis('image');
  xlabel('Chosen network label');
  ylabel('Batch network label');
  xticks(1:3);
  yticks(1:3);
  xticklabels(obdata.data.temp.names);
  yticklabels(obdata.data.temp.names);
  title('Training 90%, Mapseq');
  subplot(2, 3, 6);
  imagesc(mean(cat(3, resper2(:, C).barseq_corr), 3, 'omitnan'));
  colorbar;
  axis('image');
  xlabel('Chosen network label');
  ylabel('Batch network label');
  xticks(1:3);
  yticks(1:3);
  xticklabels(obdata.data.temp.names);
  yticklabels(obdata.data.temp.names);
  title('Training 90%, Barseq');
end
sgtitle('Comparison wrt the chosen network; max probability (Batch avg. of Label Correlation)');


FSIZE = 20;
NET = 3;

figure6(23) = figure;
set(figure6(14), 'name', 'Simplified fig 1');
subplot(1, 3, 1);
sto_conf = round(mean(cat(3, nn_sets(:, NET).train_conf), 3));
confusionchart(sto_conf, obdata.data.temp.names);
title('No training split');
set(gca, 'Fontsize', FSIZE);
subplot(1, 3, 3);
sto_conf = round(mean(cat(3, nn_sets_test(:, NET).test_conf), 3));
confusionchart(sto_conf, obdata.data.temp.names);
title('.1/.9 Split; testing');
set(gca, 'Fontsize', FSIZE);
subplot(1, 3, 2);
sto_conf = round(mean(cat(3, nn_sets_test(:, NET).train_conf), 3));
confusionchart(sto_conf, obdata.data.temp.names);
title('.1/.9 Split; training');
set(gca, 'Fontsize', FSIZE);
sgtitle('Comparison of classifier performance; maximum probability classification');

figure6(24) = figure;
set(figure6(15), 'name', 'Simplified fig 2');
subplot(1, 3, 1);
sto_conf = round(mean(cat(3, nn_sets(:, NET).mitral_train_conf), 3));
confusionchart(sto_conf, {'Negative', 'Positive'});
title('No training split');
set(gca, 'Fontsize', FSIZE);
subplot(1, 3, 3);
sto_conf = round(mean(cat(3, nn_sets_test(:, NET).mitral_test_conf), 3));
confusionchart(sto_conf, {'Negative', 'Positive'});
title('.1/.9 Split; testing');
set(gca, 'Fontsize', FSIZE);
subplot(1, 3, 2);
sto_conf = round(mean(cat(3, nn_sets_test(:, NET).mitral_train_conf), 3));
confusionchart(sto_conf, {'Negative', 'Positive'});
title('.1/.9 Split; training');
set(gca, 'Fontsize', FSIZE);
sgtitle('Comparison of classifier performance; probability threshold classification');

figure6(25) = figure;
set(figure6(16), 'name', 'Simplified fig 3');
subplot(1, 3, 1);
sto_conf = round(mean(cat(3, nn_sets(:, NET).mitral_train_conf), 3));
confusionchart(sto_conf, {'Negative', 'Positive'});
title('No training split');
set(gca, 'Fontsize', FSIZE);
subplot(1, 3, 3);
sto_conf = round(mean(cat(3, resper2(:, NET).mitral_test_conf), 3));
confusionchart(sto_conf, {'Negative', 'Positive'});
title('.2/.8 Split; testing');
set(gca, 'Fontsize', FSIZE);
subplot(1, 3, 2);
sto_conf = round(mean(cat(3, resper2(:, NET).mitral_train_conf), 3));
confusionchart(sto_conf, {'Negative', 'Positive'});
title('.2/.8 Split; training');
set(gca, 'Fontsize', FSIZE);
sgtitle('Comparison of classifier performance; probability threshold classification');

savefig(figure6(6:25), 'data/figures/classifier_nnDistribution.fig');