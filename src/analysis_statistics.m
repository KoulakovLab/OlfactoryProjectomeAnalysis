% Analysis plan
clc;

% Change the cutoff to 85%
% Also print out the values for the slopes
% Retrieve from dropbox the classifier xio used
% Bootstrap functions should take a downsampling parameter;
%   - To show that the size of our dataset is significant
%   - Mitral cells only

%---PARAMETERS---%
% Bootstrap sample size
BOOT = 1000;
NUMS = 1000;
PC_SOMA = 1;
IPR_THRESH = 8;
ELIM_AMNT = .85;
STARTTIME = tic;

% Load all the neccessary data
load('data/current.mat', ...
  'aliOB', 'data_ob', ...
  'aliPC', 'aPCAll', 'data_pc', 'data_pcuf');
load('data/classifier.mat', 'aliOB_classifier', ...
  'aliOB_classifier_vers', 'REG_COR', 'nn_sets');
load('data/raw/pcidxall.mat', 'pcidxall');

% Some figures
figure0 = gobjects(1);

%% Classification: Cell-type classificaton
aliOB.data.classifier = aliOB_classifier;
aliOB.data.classifier_vers = aliOB_classifier_vers;

% Run classification on the nn output of the templates
aliOB.data.temp.prob = aliOB.data.classifier.net(aliOB.data.temp.prjReg')';
[~, aliOB.data.temp.classId] = max(aliOB.data.temp.prob, [], 2);

% 2D plotting stats
aliOB.data.classifier.center = mean(REG_COR, 1);
aliOB.data.classifier.pts = REG_COR;

% Get condition counts
realBool = aliOB.data.temp.realId == (1:aliOB.data.temp.nIds);
aliOB.data.classifier.conPos = sum( realBool, 1);
aliOB.data.classifier.conNeg = sum(~realBool, 1);

% Get measured counts
measBool = aliOB.data.temp.classId == (1:aliOB.data.temp.nIds);
aliOB.data.classifier.claPos = sum( measBool, 1);
aliOB.data.classifier.claNeg = sum(~measBool, 1);

% Get number of correct and type I/II errors
aliOB.data.classifier.tp = sum(( realBool) & ( measBool), 1);
aliOB.data.classifier.fp = sum((~realBool) & ( measBool), 1);
aliOB.data.classifier.tn = sum((~realBool) & (~measBool), 1);
aliOB.data.classifier.fn = sum(( realBool) & (~measBool), 1);

% Get info for ROC curves for template data
[ aliOB.data.temp.fpr, ...
  aliOB.data.temp.tpr, ...
  aliOB.data.temp.auc, ...
  aliOB.data.temp.threshold] = aux.getROC( ...
  aliOB.data.temp.prob, ...
  aliOB.data.temp.realId == (1:aliOB.data.temp.nIds));

% Get the confusion matrix for template data/
aliOB.data.temp.confusionMatrix = confusionmat( ...
  aliOB.data.temp.realId, ...
  aliOB.data.temp.classId);

% Run this classifier on all the data points as well.
aliOB.brcName{2} = aliOB.data.temp.names;
aliOB.data.brc_nnetProb = cell(aliOB.nData, 1);
for d = 1:aliOB.nData
  % Get network results
  aliOB.data.brc_nnetProb{d} = ...
    aliOB.data.classifier.net(aliOB.prjRegSum{d}')';
  % Do the classifier, and add the result to barcode id
  [~, aliOB.brcId{d}(:, 2)] = max(aliOB.data.brc_nnetProb{d}, [], 2);
end

%% Organize: OB data

%-------------------------%
%---Olfactory Bulb Data---%
%-------------------------%

% OB injection data; the entire dataset merged
obdata = aliOB.getDataset('Merged');
obdata.name = 'OB injection';
obdata.brcName = aliOB.brcName;
[sto_max, sto_prj] = max(obdata.prjImg, [], 2);
[~, sto_ord] = sortrows([sto_prj, sto_max]);
sto_prj = sto_prj(sto_ord);
obdata.srcImg = obdata.srcImg(sto_ord, :);
obdata.prjImg = obdata.prjImg(sto_ord, :);
obdata.brcId = obdata.brcId(sto_ord, :);
obdata.brcCor = obdata.brcCor(sto_ord, :);
obdata.data.brcVis = sto_prj;
obdata.data.brcProb = aliOB.data.brc_nnetProb{end}(sto_ord, :);
obdata.data.brcDV = aliOB.data.brcDV{end}(sto_ord, :);
obdata.data.brcSoma = aliOB.data.brcSoma{end}(sto_ord, :);

% After this; place the pc indices here
obdata.brcName{3} = {'Broadly projecting'; 'Narrowly projecting'; ...
  'Non-mitral'; 'No PC projection'};
obdata.brcId = [obdata.brcId, zeros(obdata.nBrc, 1)];
obdata.brcId(pcidxall == 1, 3) = 1;
obdata.brcId(pcidxall == 2, 3) = 2;
obdata.brcId(pcidxall == 0, 3) = 3;
obdata.brcId(isnan(pcidxall), 3) = 4;
% Write this data back to aliOB and data_ob
aliOB.brcName = obdata.brcName;
for i = 1:aliOB.nData
  sto_ind = aux.findRowRemap(obdata.prjImg, aliOB.prjImg{i});
  aliOB.brcId{i} = [aliOB.brcId{i}, obdata.brcId(sto_ind, 3)];
end
for i = 1:length(data_ob)
  data_ob(i).brcName = obdata.brcName;
  sto_ind = aux.findRowRemap(obdata.prjImg, data_ob(i).prjImg);
  data_ob(i).brcId = [data_ob(i).brcId, obdata.brcId(sto_ind, 3)];
end

% An ordering;
% Sort according to cell type; then VIS: Very Important Slice
% Cell type ordering should be deep, tufted then mitral
obdata_type = aliOB.getDataset('Merged');
obdata_type.name = 'OB injection';
obdata_type.brcName = aliOB.brcName;
[~, sto_prj] = max(obdata_type.prjImg, [], 2);
sto_this = zeros(obdata_type.nBrc, 1);
sto_t = find(strcmp(obdata_type.data.temp.names, 'Tufted'), 1);
sto_d = find(strcmp(obdata_type.data.temp.names, 'Deep'  ), 1);
sto_m = find(strcmp(obdata_type.data.temp.names, 'Mitral'), 1);
sto_this(obdata_type.brcId(:, 2) == sto_d) = 1;
sto_this(obdata_type.brcId(:, 2) == sto_t) = 2;
sto_this(obdata_type.brcId(:, 2) == sto_m) = 3;
[~, sto_ord] = sortrows([sto_this, sto_prj]);
obdata_type.data.brcVis = sto_prj(sto_ord);
obdata_type.srcImg = obdata_type.srcImg(sto_ord, :);
obdata_type.prjImg = obdata_type.prjImg(sto_ord, :);
obdata_type.brcId = obdata_type.brcId(sto_ord, :);
obdata_type.brcCor = obdata_type.brcCor(sto_ord, :);
obdata_type.data.brcDV = aliOB.data.brcDV{end}(sto_ord, :);
obdata_type.data.brcSoma = aliOB.data.brcSoma{end}(sto_ord, :);

% An ordering; that is ordered with respect to distance
% Sort according to cell type; then to average projection distance
obdata_dist = aliOB.getDataset('Merged');
obdata_dist.name = 'OB injection';
obdata_dist.brcName = aliOB.brcName;
sto_prj = sum(obdata_dist.prjImg .* (1:obdata_dist.nPrjSli), 2) ...
  ./ sum(obdata_dist.prjImg, 2);
sto_t = find(strcmp(obdata_dist.data.temp.names, 'Tufted'), 1);
sto_d = find(strcmp(obdata_dist.data.temp.names, 'Deep'  ), 1);
sto_m = find(strcmp(obdata_dist.data.temp.names, 'Mitral'), 1);
sto_this = zeros(obdata_dist.nBrc, 1);
sto_this(obdata_dist.brcId(:, 2) == sto_d) = 1;
sto_this(obdata_dist.brcId(:, 2) == sto_t) = 2;
sto_this(obdata_dist.brcId(:, 2) == sto_m) = 3;
[~, sto_ord] = sortrows([sto_this, sto_prj]);
sto_prj = sto_prj(sto_ord);
obdata_dist.srcImg = obdata_dist.srcImg(sto_ord, :);
obdata_dist.prjImg = obdata_dist.prjImg(sto_ord, :);
obdata_dist.brcId = obdata_dist.brcId(sto_ord, :);
obdata_dist.brcCor = obdata_dist.brcCor(sto_ord, :);
obdata_dist.data.brcDis = sto_prj;
obdata_dist.data.brcDV = aliOB.data.brcDV{end}(sto_ord, :);
obdata_dist.data.brcSoma = aliOB.data.brcSoma{end}(sto_ord, :);

% OB injection with the specific partitionings
obdata_part(3, 1) = mapseqData;
% Get all measurement sets (seperate into just mitral cells)
obdata_part(1) = aliOB.getDataset('yc61' ) + aliOB.getDataset('yc65' );
obdata_part(2) = aliOB.getDataset('yc86' ) + aliOB.getDataset('yc92' );
obdata_part(3) = aliOB.getDataset('yc109') + aliOB.getDataset('yc111');
% Last partition is for the tiling figures
% Clean non-mitral barcodes
for i = length(obdata_part)
  obdata_part(i).name = [obdata_part(i).name, ' (mitral)'];
  % Get mitral class
  sto_m = find(strcmp(obdata_part(i).brcName{2}, 'Mitral'), 1);
  sto_rm = obdata_part(i).brcId(:, 2) ~= sto_m;
  obdata_part(i).srcImg(sto_rm, :) = [];
  obdata_part(i).prjImg(sto_rm, :) = [];
  obdata_part(i).brcId(sto_rm, :) = [];
end

% Create mitral/tufted/deep cell datasets
sto_t = find(strcmp(obdata.data.temp.names, 'Tufted'), 1);
sto_d = find(strcmp(obdata.data.temp.names, 'Deep'  ), 1);
sto_m = find(strcmp(obdata.data.temp.names, 'Mitral'), 1);
% Mitral with lower end removed
sto = ...
  (obdata.brcId(:, 2) == sto_m) & ...
  (obdata.data.brcProb(:, sto_m) >= ELIM_AMNT);
% Mitral
obdata_mitral = mapseqData;
obdata_mitral.name = [obdata.name, ' (mitral)'];
obdata_mitral.srcImg = obdata.srcImg(sto, :);
obdata_mitral.srcRegName = obdata.srcRegName;
obdata_mitral.nSrcRegSli = obdata.nSrcRegSli;
obdata_mitral.prjImg = obdata.prjImg(sto, :);
obdata_mitral.prjRegName = obdata.prjRegName;
obdata_mitral.nPrjRegSli = obdata.nPrjRegSli;
obdata_mitral.brcId = obdata.brcId(sto, :);
obdata_mitral.brcName = obdata.brcName;
obdata_mitral.brcCor = obdata.brcCor(sto, :);
obdata_mitral.data.temp = obdata.data.temp;
obdata_mitral.data.classifier = obdata.data.classifier;
% Tufted
sto = ...
  (obdata.brcId(:, 2) == sto_t) & ...
  (obdata.data.brcProb(:, sto_t) >= ELIM_AMNT);
obdata_tufted = mapseqData;
obdata_tufted.name = [obdata.name, ' (tufted)'];
obdata_tufted.srcImg = obdata.srcImg(sto, :);
obdata_tufted.srcRegName = obdata.srcRegName;
obdata_tufted.nSrcRegSli = obdata.nSrcRegSli;
obdata_tufted.prjImg = obdata.prjImg(sto, :);
obdata_tufted.prjRegName = obdata.prjRegName;
obdata_tufted.nPrjRegSli = obdata.nPrjRegSli;
obdata_tufted.brcId = obdata.brcId(sto, :);
obdata_tufted.brcName = obdata.brcName;
obdata_tufted.brcCor = obdata.brcCor(sto, :);
obdata_tufted.data.temp = obdata.data.temp;
obdata_tufted.data.classifier = obdata.data.classifier;
% Deep
sto = ...
  (obdata.brcId(:, 2) == sto_d) & ...
  (obdata.data.brcProb(:, sto_d) >= ELIM_AMNT);
obdata_deep = mapseqData;
obdata_deep.name = [obdata.name, ' (deep)'];
obdata_deep.srcImg = obdata.srcImg(sto, :);
obdata_deep.srcRegName = obdata.srcRegName;
obdata_deep.nSrcRegSli = obdata.nSrcRegSli;
obdata_deep.prjImg = obdata.prjImg(sto, :);
obdata_deep.prjRegName = obdata.prjRegName;
obdata_deep.nPrjRegSli = obdata.nPrjRegSli;
obdata_deep.brcId = obdata.brcId(sto, :);
obdata_deep.brcName = obdata.brcName;
obdata_deep.brcCor = obdata.brcCor(sto, :);
obdata_deep.data.temp = obdata.data.temp;
obdata_deep.data.classifier = obdata.data.classifier;

% Ipr Seperation
figure0(1) = figure;
set(figure0(1), 'name', 'IPR cutoff');
histogram(obdata_mitral.brcPrjIpr);
hold('on');
plot(ones(1, 2) * IPR_THRESH, ylim, '-k', 'LineWidth', 3);
hold('off');
savefig(figure0, 'data/figures/figure_0.fig');
savefig(figure0, '/home/sbp/Maestral/OBmapseq/batufigs/figure_0.fig');

% Low-ipr
% lowipr = obdata_mitral.brcPrjIpr <= IPR_THRESH;
lowipr = obdata_mitral.brcId(:, 3) == 2;
obdata_lowipr = mapseqData;
obdata_lowipr.name = [obdata_mitral.name(1:(end-1)), ...
  ', (', obdata_mitral.brcName{3}{2}, ')'];
obdata_lowipr.srcImg = obdata_mitral.srcImg(lowipr, :);
obdata_lowipr.srcRegName = obdata_mitral.srcRegName;
obdata_lowipr.nSrcRegSli = obdata_mitral.nSrcRegSli;
obdata_lowipr.prjImg = obdata_mitral.prjImg(lowipr, :);
obdata_lowipr.prjRegName = obdata_mitral.prjRegName;
obdata_lowipr.nPrjRegSli = obdata_mitral.nPrjRegSli;
obdata_lowipr.brcId = obdata_mitral.brcId(lowipr, :);
obdata_lowipr.brcName = obdata_mitral.brcName;
obdata_lowipr.brcCor = obdata_mitral.brcCor(lowipr, :);
obdata_lowipr.data.temp = obdata_mitral.data.temp;
obdata_lowipr.data.classifier = obdata_mitral.data.classifier;
% High-ipr
% highipr = ~lowipr;
highipr = obdata_mitral.brcId(:, 3) == 1;
obdata_highipr = mapseqData;
obdata_highipr.name = [obdata_mitral.name(1:(end-1)), ...
  ', (', obdata_mitral.brcName{3}{1}, ')'];
obdata_highipr.srcImg = obdata_mitral.srcImg(highipr, :);
obdata_highipr.srcRegName = obdata_mitral.srcRegName;
obdata_highipr.nSrcRegSli = obdata_mitral.nSrcRegSli;
obdata_highipr.prjImg = obdata_mitral.prjImg(highipr, :);
obdata_highipr.prjRegName = obdata_mitral.prjRegName;
obdata_highipr.nPrjRegSli = obdata_mitral.nPrjRegSli;
obdata_highipr.brcId = obdata_mitral.brcId(highipr, :);
obdata_highipr.brcName = obdata_mitral.brcName;
obdata_highipr.brcCor = obdata_mitral.brcCor(highipr, :);
obdata_highipr.data.temp = obdata_mitral.data.temp;
obdata_highipr.data.classifier = obdata_mitral.data.classifier;

% All-but-one 'abone' dataset
obdata_abone = ...
  + aliOB.getDataset('yc61') ...
  + aliOB.getDataset('yc65') ...
  + aliOB.getDataset('yc86') ...
  + aliOB.getDataset('yc92') ...
  + aliOB.getDataset('yc109');
obdata_abone.data.temp = aliOB.data.temp;
obdata_abone.data.classifier = aliOB.data.classifier;
obdata_abone.brcName = aliOB.brcName;
obdata_abone.brcName{1}(end) = [];
obdata_abone.name = [obdata.name, ' -yc111'];
[sto_max, sto_prj] = max(obdata_abone.prjImg, [], 2);
[~, sto_ord] = sortrows([sto_prj, sto_max]);
obdata_abone.srcImg = obdata_abone.srcImg(sto_ord, :);
obdata_abone.prjImg = obdata_abone.prjImg(sto_ord, :);
obdata_abone.brcId = obdata_abone.brcId(sto_ord, :);
obdata_abone.brcCor = obdata_abone.brcCor(sto_ord, :);
obdata_abone.data.brcVis = sto_prj(sto_ord);
obdata_abone.data.brcProb = aliOB.data.classifier.net( ...
  obdata_abone.prjRegSum')';

% Create tufted/mitral/deep/low/high seperation for all-but-one
% Mitral cells
sto = ...
  (obdata_abone.brcId(:, 2) == sto_m) & ...
  (obdata_abone.data.brcProb(:, sto_m) >= ELIM_AMNT);
obdata_abone_mitral = mapseqData;
obdata_abone_mitral.name = [obdata_abone.name, ' (mitral)'];
obdata_abone_mitral.srcImg = obdata_abone.srcImg(sto, :);
obdata_abone_mitral.srcRegName = obdata_abone.srcRegName;
obdata_abone_mitral.nSrcRegSli = obdata_abone.nSrcRegSli;
obdata_abone_mitral.prjImg = obdata_abone.prjImg(sto, :);
obdata_abone_mitral.prjRegName = obdata_abone.prjRegName;
obdata_abone_mitral.nPrjRegSli = obdata_abone.nPrjRegSli;
obdata_abone_mitral.brcId = obdata_abone.brcId(sto, :);
obdata_abone_mitral.brcName = obdata_abone.brcName;
obdata_abone_mitral.brcCor = obdata_abone.brcCor(sto, :);
obdata_abone_mitral.data.temp = obdata_abone.data.temp;
obdata_abone_mitral.data.classifier = obdata_abone.data.classifier;
% Tufted
sto = ...
  (obdata_abone.brcId(:, 2) == sto_t) & ...
  (obdata_abone.data.brcProb(:, sto_t) >= ELIM_AMNT);
obdata_abone_tufted = mapseqData;
obdata_abone_tufted.name = [obdata_abone.name, ' (tufted)'];
obdata_abone_tufted.srcImg = obdata_abone.srcImg(sto, :);
obdata_abone_tufted.srcRegName = obdata_abone.srcRegName;
obdata_abone_tufted.nSrcRegSli = obdata_abone.nSrcRegSli;
obdata_abone_tufted.prjImg = obdata_abone.prjImg(sto, :);
obdata_abone_tufted.prjRegName = obdata_abone.prjRegName;
obdata_abone_tufted.nPrjRegSli = obdata_abone.nPrjRegSli;
obdata_abone_tufted.brcId = obdata_abone.brcId(sto, :);
obdata_abone_tufted.brcName = obdata_abone.brcName;
obdata_abone_tufted.brcCor = obdata_abone.brcCor(sto, :);
obdata_abone_tufted.data.temp = obdata_abone.data.temp;
obdata_abone_tufted.data.classifier = obdata_abone.data.classifier;
% Deep
sto = ...
  (obdata_abone.brcId(:, 2) == sto_d) & ...
  (obdata_abone.data.brcProb(:, sto_d) >= ELIM_AMNT);
obdata_abone_deep = mapseqData;
obdata_abone_deep.name = [obdata_abone.name, ' (tufted)'];
obdata_abone_deep.srcImg = obdata_abone.srcImg(sto, :);
obdata_abone_deep.srcRegName = obdata_abone.srcRegName;
obdata_abone_deep.nSrcRegSli = obdata_abone.nSrcRegSli;
obdata_abone_deep.prjImg = obdata_abone.prjImg(sto, :);
obdata_abone_deep.prjRegName = obdata_abone.prjRegName;
obdata_abone_deep.nPrjRegSli = obdata_abone.nPrjRegSli;
obdata_abone_deep.brcId = obdata_abone.brcId(sto, :);
obdata_abone_deep.brcName = obdata_abone.brcName;
obdata_abone_deep.brcCor = obdata_abone.brcCor(sto, :);
obdata_abone_deep.data.temp = obdata_abone.data.temp;
obdata_abone_deep.data.classifier = obdata_abone.data.classifier;
% Low IPR
sto = obdata_abone_mitral.brcPrjIpr <= IPR_THRESH;
obdata_abone_lowipr = mapseqData;
obdata_abone_lowipr.name = [obdata_abone_mitral.name(1:(end-1)), ...
  ', IPR <=', num2str(IPR_THRESH), ')'];
obdata_abone_lowipr.srcImg = obdata_abone_mitral.srcImg(sto, :);
obdata_abone_lowipr.srcRegName = obdata_abone_mitral.srcRegName;
obdata_abone_lowipr.nSrcRegSli = obdata_abone_mitral.nSrcRegSli;
obdata_abone_lowipr.prjImg = obdata_abone_mitral.prjImg(sto, :);
obdata_abone_lowipr.prjRegName = obdata_abone_mitral.prjRegName;
obdata_abone_lowipr.nPrjRegSli = obdata_abone_mitral.nPrjRegSli;
obdata_abone_lowipr.brcId = obdata_abone_mitral.brcId(sto, :);
obdata_abone_lowipr.brcName = obdata_abone_mitral.brcName;
obdata_abone_lowipr.brcCor = obdata_abone_mitral.brcCor(sto, :);
obdata_abone_lowipr.data.temp = obdata_abone_mitral.data.temp;
obdata_abone_lowipr.data.classifier = obdata_abone_mitral.data.classifier;
% High IPR
sto = obdata_abone_mitral.brcPrjIpr > IPR_THRESH;
obdata_abone_highipr = mapseqData;
obdata_abone_highipr.name = [obdata_abone_mitral.name(1:(end-1)), ...
  ', IPR >', num2str(IPR_THRESH), ')'];
obdata_abone_highipr.srcImg = obdata_abone_mitral.srcImg(sto, :);
obdata_abone_highipr.srcRegName = obdata_abone_mitral.srcRegName;
obdata_abone_highipr.nSrcRegSli = obdata_abone_mitral.nSrcRegSli;
obdata_abone_highipr.prjImg = obdata_abone_mitral.prjImg(sto, :);
obdata_abone_highipr.prjRegName = obdata_abone_mitral.prjRegName;
obdata_abone_highipr.nPrjRegSli = obdata_abone_mitral.nPrjRegSli;
obdata_abone_highipr.brcId = obdata_abone_mitral.brcId(sto, :);
obdata_abone_highipr.brcName = obdata_abone_mitral.brcName;
obdata_abone_highipr.brcCor = obdata_abone_mitral.brcCor(sto, :);
obdata_abone_highipr.data.temp = obdata_abone_mitral.data.temp;
obdata_abone_highipr.data.classifier = obdata_abone_mitral.data.classifier;

% Seperate data_ob into respective datasets
data_ob_mitral(length(data_ob), 1) = mapseqData;
data_ob_tufted(length(data_ob), 1) = mapseqData;
data_ob_deep(length(data_ob), 1) = mapseqData;
data_ob_lowipr(length(data_ob), 1) = mapseqData;
data_ob_highipr(length(data_ob), 1) = mapseqData;
for i = 1:length(data_ob)
  % Get classifier on each element in data_ob
  data_ob(i).data.temp = aliOB.data.temp;
  data_ob(i).data.classifier = aliOB.data.classifier;
  data_ob(i).brcCor = aliOB.brcCor{i};
  % Run the classifier
  data_ob(i).data.brcProb = ...
    data_ob(i).data.classifier.net(data_ob(i).prjRegSum')';
  [~, data_ob(i).brcId(:, 4)] = max( ...
    data_ob(i).data.brcProb, [], 2);
  data_ob(i).brcName{4} = data_ob(i).data.temp.names';
  % Mitral
  sto = ...
    (data_ob(i).brcId(:, 4) == sto_m) & ...
    (data_ob(i).data.brcProb(:, sto_m) >= ELIM_AMNT);
  data_ob_mitral(i).name = [data_ob(i).name, ' (mitral)'];
  data_ob_mitral(i).srcImg = data_ob(i).srcImg(sto, :);
  data_ob_mitral(i).prjImg = data_ob(i).prjImg(sto, :);
  data_ob_mitral(i).brcId  = data_ob(i).brcId( sto, :);
  data_ob_mitral(i).srcRegName = data_ob(i).srcRegName;
  data_ob_mitral(i).nSrcRegSli = data_ob(i).nSrcRegSli;
  data_ob_mitral(i).prjRegName = data_ob(i).prjRegName;
  data_ob_mitral(i).nPrjRegSli = data_ob(i).nPrjRegSli;
  data_ob_mitral(i).brcName = data_ob(i).brcName;
  data_ob_mitral(i).brcCor = data_ob(i).brcCor(sto, :);
  data_ob_mitral(i).data.temp = data_ob(i).data.temp;
  data_ob_mitral(i).data.classifier = data_ob(i).data.classifier;
  % Tufted
  sto = ...
    (data_ob(i).brcId(:, 4) == sto_t) & ...
    (data_ob(i).data.brcProb(:, sto_t) >= ELIM_AMNT);
  data_ob_tufted(i).name = [data_ob(i).name, ' (tufted)'];
  data_ob_tufted(i).srcImg = data_ob(i).srcImg(sto, :);
  data_ob_tufted(i).prjImg = data_ob(i).prjImg(sto, :);
  data_ob_tufted(i).brcId  = data_ob(i).brcId( sto, :);
  data_ob_tufted(i).srcRegName = data_ob(i).srcRegName;
  data_ob_tufted(i).nSrcRegSli = data_ob(i).nSrcRegSli;
  data_ob_tufted(i).prjRegName = data_ob(i).prjRegName;
  data_ob_tufted(i).nPrjRegSli = data_ob(i).nPrjRegSli;
  data_ob_tufted(i).brcName = data_ob(i).brcName;
  data_ob_tufted(i).brcCor = data_ob(i).brcCor(sto, :);
  data_ob_tufted(i).data.temp = data_ob(i).data.temp;
  data_ob_tufted(i).data.classifier = data_ob(i).data.classifier;
  % Deep
  sto = ...
    (data_ob(i).brcId(:, 4) == sto_d) & ...
    (data_ob(i).data.brcProb(:, sto_d) >= ELIM_AMNT);
  data_ob_tufted(i).name = [data_ob(i).name, ' (tufted)'];
  data_ob_tufted(i).srcImg = data_ob(i).srcImg(sto, :);
  data_ob_tufted(i).prjImg = data_ob(i).prjImg(sto, :);
  data_ob_tufted(i).brcId  = data_ob(i).brcId( sto, :);
  data_ob_tufted(i).srcRegName = data_ob(i).srcRegName;
  data_ob_tufted(i).nSrcRegSli = data_ob(i).nSrcRegSli;
  data_ob_tufted(i).prjRegName = data_ob(i).prjRegName;
  data_ob_tufted(i).nPrjRegSli = data_ob(i).nPrjRegSli;
  data_ob_tufted(i).brcName = data_ob(i).brcName;
  data_ob_tufted(i).brcCor = data_ob(i).brcCor(sto, :);
  data_ob_tufted(i).data.temp = data_ob(i).data.temp;
  data_ob_tufted(i).data.classifier = data_ob(i).data.classifier;
  % Low-IPR
  sto = data_ob_mitral(i).brcPrjIpr <= IPR_THRESH;
  data_ob_lowipr(i).name = [data_ob_mitral(i).name(1:(end-1)), ...
    ', IPR <= ', num2str(IPR_THRESH), ')'];
  data_ob_lowipr(i).srcImg = data_ob_mitral(i).srcImg(sto, :);
  data_ob_lowipr(i).prjImg = data_ob_mitral(i).prjImg(sto, :);
  data_ob_lowipr(i).brcId  = data_ob_mitral(i).brcId( sto, :);
  data_ob_lowipr(i).srcRegName = data_ob_mitral(i).srcRegName;
  data_ob_lowipr(i).nSrcRegSli = data_ob_mitral(i).nSrcRegSli;
  data_ob_lowipr(i).prjRegName = data_ob_mitral(i).prjRegName;
  data_ob_lowipr(i).nPrjRegSli = data_ob_mitral(i).nPrjRegSli;
  data_ob_lowipr(i).brcName = data_ob_mitral(i).brcName;
  data_ob_lowipr(i).brcCor = data_ob_mitral(i).brcCor(sto, :);
  data_ob_lowipr(i).data.temp = data_ob_mitral(i).data.temp;
  data_ob_lowipr(i).data.classifier = data_ob_mitral(i).data.classifier;
  % High-IPR
  sto = data_ob_mitral(i).brcPrjIpr > IPR_THRESH;
  data_ob_highipr(i).name = [data_ob_mitral(i).name(1:(end-1)), ...
    ', IPR > ', num2str(IPR_THRESH), ')'];
  data_ob_highipr(i).srcImg = data_ob_mitral(i).srcImg(sto, :);
  data_ob_highipr(i).prjImg = data_ob_mitral(i).prjImg(sto, :);
  data_ob_highipr(i).brcId  = data_ob_mitral(i).brcId( sto, :);
  data_ob_highipr(i).srcRegName = data_ob_mitral(i).srcRegName;
  data_ob_highipr(i).nSrcRegSli = data_ob_mitral(i).nSrcRegSli;
  data_ob_highipr(i).prjRegName = data_ob_mitral(i).prjRegName;
  data_ob_highipr(i).nPrjRegSli = data_ob_mitral(i).nPrjRegSli;
  data_ob_highipr(i).brcName = data_ob_mitral(i).brcName;
  data_ob_highipr(i).brcCor = data_ob_mitral(i).brcCor(sto, :);
  data_ob_highipr(i).data.temp = data_ob_mitral(i).data.temp;
  data_ob_highipr(i).data.classifier = data_ob_mitral(i).data.classifier;
end

%--------------------------%
%---Piriform Cortex data---%
%--------------------------%

% PC injection data; with the original set
pcdata = aliPC.getDataset('Merged');
pcdata.name = 'PC injection';
pcdata.brcName = aliPC.brcName;

% Sort according to VIS: Very Important Slice
[~, sto_src] = max(pcdata.srcImg, [], 2);
[pcdata.data.brcVis, sto_ord] = sort(sto_src);
pcdata.srcImg = pcdata.srcImg(sto_ord, :);
pcdata.prjImg = pcdata.prjImg(sto_ord, :);
pcdata.brcId  = pcdata.brcId( sto_ord, :);
pcdata.data.brcOrigin = pcdata.srcImg == max(pcdata.srcImg, [], 2);

% PC injection data; with the entirety of the data
pcdata_nofilt = aPCAll.getDataset('Merged');
pcdata_nofilt.name = 'PC injection; no projection filter';
pcdata_nofilt.brcName = aPCAll.brcName;

% Sort according to VIS: Very Important Slice
[~, sto_src] = max(pcdata_nofilt.srcImg, [], 2);
[pcdata_nofilt.data.brcVis, sto_ord] = sort(sto_src);
pcdata_nofilt.srcImg = pcdata_nofilt.srcImg(sto_ord, :);
pcdata_nofilt.prjImg = pcdata_nofilt.prjImg(sto_ord, :);
pcdata_nofilt.brcId  = pcdata_nofilt.brcId( sto_ord, :);
pcdata_nofilt.data.brcOrigin = ...
  pcdata_nofilt.srcImg == max(pcdata_nofilt.srcImg, [], 2);

%% Analysis: Tiling

%------------------------------%
%---Tiling: Exponential Fits---%
%------------------------------%

% The exponential function
FN_OPTNS = fitoptions('exp1');
FN_OPTNS.StartPoint = 1;
FN_FUNCT = @(a, x) aux.limexpn(a, x);

% Repeat the analysis on all the following datasets
sto_ana = [ obdata, ...       % MERGED ALL
  obdata_mitral,  ...         % Mitral cells of all
  obdata_tufted,  ...         % Tufted cells of all
  obdata_deep,  ...           % Deep cells of all
  obdata_highipr, ...         % High-ipr mitral cells of all
  obdata_lowipr, ...          % Low-ipr mitral cells of all
  obdata_abone, ...           % ALL BUT ONE (-YC111)
  obdata_abone_mitral, obdata_abone_tufted, obdata_abone_deep, ...
  obdata_abone_highipr, obdata_abone_lowipr, ...
  data_ob', ...               % BRAINS: SEPERATED
  data_ob_mitral', data_ob_tufted', data_ob_deep', ...
  data_ob_lowipr', data_ob_highipr'];
fprintf('Fitting exponentials to tiling regions');


for i = 1:length(sto_ana)
  aux.progressbar(i, length(sto_ana));
  % Going to do fits on all the regions
  sto_ana(i).data.expFit = cell(sto_ana(i).nPrjReg, 1);
  sto_ana(i).data.expVal = nan(sto_ana(i).nPrjReg, 1);
  % Get the VIS
  [~, sto_ana(i).data.brcVis] = max(sto_ana(i).prjImg, [], 2);
  % Run for all available regions
  for r = 1:sto_ana(i).nPrjReg
    % Get the list of all the barcodes with maxima in this region
    brcIds = find(any( ...
      sto_ana(i).data.brcVis == sto_ana(i).prjRegInd{r}, 2));
    % Write the VIS regions her
    % Skip if not tileable
    if isempty(brcIds)
      continue;
    end
    brcVis = sto_ana(i).data.brcVis(brcIds);
    % Get xy values appropriate for the function
    brcXvals = linspace(0, 1, length(brcIds))';
    brcYvals = ...
      (brcVis - sum(sto_ana(i).nPrjRegSli(1:(r-1))) - 1 ) ...
      / (sto_ana(i).nPrjRegSli(r) - 1);
    % Do the fits in xy space
    sto_fit = fit(brcXvals, brcYvals, fittype(FN_FUNCT), FN_OPTNS);
    sto_ana(i).data.expVal(r) = coeffvalues(sto_fit);
    % Create surrogate values to plot
    fitXvals = linspace(0, 1, 100)';
    fitYvals = FN_FUNCT(sto_ana(i).data.expVal(r), fitXvals);
    fitBrc = linspace(min(brcIds, [], 'all'), max(brcIds, [], 'all'), 100)';
    fitSli = fitYvals * (sto_ana(i).nPrjRegSli(r) - 1) ...
      + sum(sto_ana(i).nPrjRegSli(1:(r-1))) + 1;
    sto_ana(i).data.expFit{r} = [fitSli, fitBrc];
  end
end

% Also calculate cross-correlation between orderings;
obdata_mitral.data.brcOrd = ...
  zeros(obdata_mitral.nBrc, obdata_mitral.nPrjReg);
for pr = 1:obdata_mitral.nPrjReg
  sto_img = obdata_mitral.getImgPrjReg(pr);
  [~, obdata_mitral.data.brcOrd(:, pr)] = max(sto_img, [], 2);
end
[obdata_mitral.data.regOrdCorr, obdata_mitral.data.regOrdPval] = ...
  corr(obdata_mitral.data.brcOrd);

%---Tiling: Lambda function for different values of splittings
SPLIT_NUM = 20;
SPLIT_INT = true;
% Repeat the analysis on all the following datasets
sto_ana = [obdata_mitral,  obdata_abone_mitral, data_ob_mitral'];
fprintf('Fitting exponentials using a range of IPR thresholds...');
for i = 1:length(sto_ana)
  aux.progressbar(i, length(sto_ana));
  sto_sli = sto_ana(i).nPrjRegSli(2);
  sto_sin = sto_ana(i).prjRegInd{2};
  % Get the image to this region
  sto_img = sto_ana(i).prjImg(:, sto_sin) ./ max(sto_ana(i).prjImg, [], 2);
  % Find the VIS and IPR for barcodes in APC; and restrict to here
  sto_brc = any(sto_ana(i).data.brcVis == sto_sin, 2);
  sto_img = sto_img(sto_brc, :);
  sto_vis = sto_ana(i).data.brcVis(sto_brc, :);
  sto_ipr = sto_ana(i).brcPrjIpr(sto_brc, :);
  % Get a range of IPR's; leaving at least 10% of the barcodes in each set
  if SPLIT_INT
    sto_ana(i).data.exp_iprVal = ...
      floor(prctile(sto_ipr, 10)):ceil(prctile(sto_ipr, 90));
  else
    sto_ana(i).data.exp_iprVal = linspace( ...
      prctile(sto_ipr, 10), prctile(sto_ipr, 90), SPLIT_NUM)';
  end
  sto_snu = length(sto_ana(i).data.exp_iprVal);
  sto_ana(i).data.exp_nLow = zeros(sto_snu, 1);
  sto_ana(i).data.exp_nHigh = zeros(sto_snu, 1);
  sto_ana(i).data.exp_lowGamma = zeros(sto_snu, 1);
  sto_ana(i).data.exp_highGamma = zeros(sto_snu, 1);
  % Background image of varios partitionings, to be stacked horizontally
  sto_ana(i).data.exp_imgBkg = zeros(sum(sto_brc), sto_snu * sto_sli);
  % * Corresponding X and Y values of the appropriate fits; in a cell array
  sto_ana(i).data.exp_imgXval = cell(sto_snu, 2);
  sto_ana(i).data.exp_imgYval = cell(sto_snu, 2);
  for s = 1:length(sto_ana(i).data.exp_iprVal)
      % Split the dataset
      sto_ind = sto_ipr > sto_ana(i).data.exp_iprVal(s);
      % Get frequency of the barcodes (First high; then low ipr)
      sto_nHi = sum(sto_ind);
      sto_nLo = sum(~sto_ind);
      sto_ana(i).data.exp_nLow(s) = sto_nLo;
      sto_ana(i).data.exp_nHigh(s)= sto_nHi;
      % Draw the background image
      sto_ana(i).data.exp_imgBkg( ...
        (1:sto_nHi), ...
        ((s - 1) * sto_sli + (1:sto_sli))) = sto_img(sto_ind, :);
      sto_ana(i).data.exp_imgBkg((1:sto_nLo) + sto_nHi, ...
        ((s - 1) * sto_sli + (1:sto_sli))) = sto_img(~sto_ind, :);
      % Fit the high-ipr gamma value
      sto_xval = linspace(0, 1, sto_ana(i).data.exp_nHigh(s))';
      sto_yval = (sto_vis(sto_ind) - sto_ana(i).nPrjRegSli(1) - 1) ...
        / (sto_ana(i).nPrjRegSli(2) - 1);
      sto_fit = fit(sto_xval, sto_yval, fittype(FN_FUNCT), FN_OPTNS);
      sto_ana(i).data.exp_highGamma(s) = coeffvalues(sto_fit);
      % Get the best fit line here; for high-ipr
      sto_xval = linspace(0, 1, 100)';
      sto_yval = FN_FUNCT(sto_ana(i).data.exp_highGamma(s), sto_xval);
      sto_ana(i).data.exp_imgXval{s, 1} = (linspace(1, sto_nHi, 100)');
      sto_ana(i).data.exp_imgYval{s, 1} = ...
        (sto_yval * (sto_sli - 1)) + ((s - 1) * sto_sli) + 1;
      % Fit the low-ipr gamma value
      sto_xval = linspace(0, 1, sto_ana(i).data.exp_nLow(s))';
      sto_yval = (sto_vis(~sto_ind) - sto_ana(i).nPrjRegSli(1) - 1) ...
        / (sto_ana(i).nPrjRegSli(2) - 1);
      sto_fit = fit(sto_xval, sto_yval, fittype(FN_FUNCT), FN_OPTNS);
      sto_ana(i).data.exp_lowGamma(s) = coeffvalues(sto_fit);
      % Get the best fit line here; for low-ipr
      sto_xval = linspace(0, 1, 100)';
      sto_yval = FN_FUNCT(sto_ana(i).data.exp_lowGamma(s), sto_xval);
      sto_ana(i).data.exp_imgXval{s, 2} = (linspace(1, sto_nLo, 100)') ...
        + sto_nHi;
      sto_ana(i).data.exp_imgYval{s, 2} = ...
        (sto_yval * (sto_sli - 1)) + ((s - 1) * sto_sli) + 1;
  end
end

%% New analysis

% Run the new analysis
fprintf('New (2020-11) analysis algorithms\n');
fprintf('Analysis on mitral cells, OB injection, PC coennervation.\n');
aux.newstatOBPC(obdata_mitral, BOOT);
for i = 1:length(data_ob_mitral)
  fprintf(['Same analysis, brain #', data_ob_mitral(i).name, '.\n']);
  aux.newstatOBPC(data_ob_mitral(i), BOOT);
end
fprintf('Analysis on narrow mitral cells, OB injection, PC coennervation.\n');
aux.newstatOBPC(obdata_lowipr, BOOT);
fprintf('Analysis on broad mitral cells, OB injection, PC coennervation.\n');
aux.newstatOBPC(obdata_highipr, BOOT);
fprintf('Analysis on mitral cells, OB injection, OT coennervation.\n');
aux.newstatOBOT(obdata_mitral, BOOT);
fprintf('Analysis on narrow mitral cells, OB injection, OT coennervation.\n');
aux.newstatOBOT(obdata_lowipr, BOOT);
fprintf('Analysis on broad mitral cells, OB injection, OT coennervation.\n');
aux.newstatOBOT(obdata_highipr, BOOT);
fprintf('Downsampling on mitral cells, OB injection, PC coennervation.');
aux.newdsOBPC(obdata_mitral, 10, NUMS, BOOT);
fprintf('Downsampling on narrow mitral cells, OB injection, PC coennervation.');
aux.newdsOBPC(obdata_lowipr, 10, NUMS, BOOT);
fprintf('Downsampling on broad mitral cells, OB injection, PC coennervation.');
aux.newdsOBPC(obdata_highipr,10, NUMS, BOOT);
fprintf('Checking neural network performances\n');
aux.newstatNN(obdata, obdata_mitral, nn_sets, BOOT);

%% Analysis: PC Prob.

%--------------------------------%
%---PC Conditional Probability---%
%--------------------------------%

fprintf('Running bootstrap on PC injection\n');
aux.newstatPCout(pcdata, BOOT);
fprintf('Running bootstrap on inter-PC injection\n');
aux.newstatPCPC(pcdata, BOOT);

%% Save records
save('data/processed.mat', ...
  'BOOT', 'PC_SOMA', 'IPR_THRESH', ...
  'aliOB', 'aliPC', ...
  'pcdata', 'pcdata_nofilt', ...
  'obdata', 'obdata_dist', 'obdata_part', ...
  'obdata_mitral', 'obdata_tufted', 'obdata_deep', 'obdata_lowipr', ...
  'obdata_highipr', 'obdata_abone', 'obdata_abone_mitral', ...
  'obdata_abone_tufted', 'obdata_abone_highipr', 'obdata_abone_lowipr', ...
  'data_ob', 'data_ob_mitral', 'data_ob_tufted', ...
  'data_ob_lowipr', 'data_ob_highipr');

RUNTIME = toc(STARTTIME);
toc(STARTTIME)

%% Organized savings
table_mitral = obdata_mitral.data.OBPC.oTable;
table_narrow = obdata_lowipr.data.OBPC.oTable;
table_broad = obdata_highipr.data.OBPC.oTable;
table_pcInj = pcdata.data.PCout.oTable;

filenames = {'data/table'};
for f = 1:length(filenames)
  filename = filenames{f};
  save([filename, 's.mat'], 'table_mitral',  'table_narrow',  'table_broad');
  writetable(table_mitral, [filename, '_mitral.csv'], 'WriteRowNames', true);
  writetable(table_narrow, [filename, '_narrow.csv'], 'WriteRowNames', true);
  writetable(table_broad,  [filename, '__broad.csv'], 'WriteRowNames', true);
  writetable(table_pcInj,  [filename, '_PCinjections.csv'], 'WriteRowNames', true);
  writetable(table_mitral, [filename, 's.xlsx'], 'WriteRowNames', true, 'Sheet', 'Mitral',        'Range', 'A1');
  writetable(table_narrow, [filename, 's.xlsx'], 'WriteRowNames', true, 'Sheet', 'Narrow',        'Range', 'A1');
  writetable(table_broad,  [filename, 's.xlsx'], 'WriteRowNames', true, 'Sheet', 'Broad',         'Range', 'A1');
  writetable(table_pcInj,  [filename, 's.xlsx'], 'WriteRowNames', true, 'Sheet', 'PC Injection',  'Range', 'A1');
end
