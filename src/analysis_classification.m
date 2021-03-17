% Classify: Cell-type classificaton

% Network parameters
NUM_NN = 500;
NN_NETWORKSIZE = {16, 20, [10, 10]};
NN_START = tic;
ELIM_AMNT = .85;
TEST_METHOD = 'threshold';
% Seed the RNJesus
rng(19);

% Load the required variables
load('data/current.mat', 'aliOB');

% Preambulatory
[CNV_TO_2D, CNV_TO_3D, REG_COR] = aux.genSimplex(2);

% Get individual ids to refer to specific cell type later
sto_mId = find(strcmp(aliOB.data.temp.names, 'Mitral'), 1);
sto_tId = find(strcmp(aliOB.data.temp.names, 'Tufted'), 1);
sto_dId = find(strcmp(aliOB.data.temp.names, 'Deep'), 1);

% Train NN's and pick the best performing one among them
fprintf('Training Neural Network classifiers . . . ');
% Get binary classification matrix
sto_temCls = double(aliOB.data.temp.realId == (1:3));

% Negative log function; but don't shoot to infinity if sthng is -
NEG_LOG = @(x) -1 * log(x + (1 / aliOB.data.temp.nBrc) * ...
  (isnan(x) | isinf(x) | (x == 0)));

% Store all tried network architectures
res = struct( ...
  'net',        cell(NUM_NN, length(NN_NETWORKSIZE)), ...
  'trainRes',   cell(NUM_NN, length(NN_NETWORKSIZE)), ...
  'prob_bar',   cell(NUM_NN, length(NN_NETWORKSIZE)), ...
  'class_id',   cell(NUM_NN, length(NN_NETWORKSIZE)), ...
  'net_loss',   cell(NUM_NN, length(NN_NETWORKSIZE)), ...
  'prob_cor',   cell(NUM_NN, length(NN_NETWORKSIZE)), ...
  'mitral_tpr', cell(NUM_NN, length(NN_NETWORKSIZE)), ...
  'mitral_fpr', cell(NUM_NN, length(NN_NETWORKSIZE)));

% Main loop for training NN's
for n = NUM_NN:-1:1
  aux.progressbar(NUM_NN - n + 1, NUM_NN);
  for i = 1:length(NN_NETWORKSIZE)
    % Initialize a neural network
    sto_inet = patternnet(NN_NETWORKSIZE{i});
    sto_inet.trainParam.showWindow = false;
    sto_inet.divideParam.trainRatio = 1;
    sto_inet.divideParam.valRatio= 0;
    sto_inet.divideParam.testRatio= 0;
    % Train a network with weighting losses
    res(n, i).title = num2str(NN_NETWORKSIZE{i});
    [res(n, i).net, res(n, i).trainRes] = train( ...
      sto_inet, aliOB.data.temp.prjReg', sto_temCls', [], [], ...
      aliOB.data.temp.nBrc ./ sqrt(sum(sto_temCls, 1)'));
    res(n, i).prob_bar = res(n, i).net(aliOB.data.temp.prjReg')';
    res(n, i).net_loss = perform( ...
      res(n, i).net, sto_temCls', res(n, i).prob_bar', ...
      1 ./ sum(sto_temCls, 1)');
    res(n, i).prob_cor = ...
      aliOB.data.temp.conv_fromBaryo(res(n, i).prob_bar);
  end
end

% Select the lowest cross-entropy loss score among all the networks
netLoss = vertcat(res(:).net_loss);
netMinLoss = min(netLoss);
netId = find(netLoss == netMinLoss, 1);
aliOB_classifier = res(netId);
[sto1, sto2] = ind2sub([NUM_NN, length(NN_NETWORKSIZE)], netId);
aliOB_classifier.netId = [sto1, sto2];

% Also among each group; pick out alternate versions
for j = length(NN_NETWORKSIZE):-1:1
  netLoss = vertcat(res(:, j).net_loss);
  netMinLoss = min(netLoss);
  netId = find(netLoss == netMinLoss, 1);
  aliOB_classifier_vers(j) = res(netId, j);
end
% Assign a name
for j = length(NN_NETWORKSIZE):-1:1
  aliOB_classifier_vers(j).name = [ ...
    'Neural network: ', num2str(NN_NETWORKSIZE{j})];
end

% Release my seed
rng('shuffle');
toc(NN_START);

%% With testing
% tic(NN_START);

% Seed the RNJesus
rng(17);

% Store all tried network architectures
resper = struct( ...
  'net',        cell(NUM_NN, length(NN_NETWORKSIZE)), ...
  'trainRes',   cell(NUM_NN, length(NN_NETWORKSIZE)), ...
  'prob_bar',   cell(NUM_NN, length(NN_NETWORKSIZE)), ...
  'class_id',   cell(NUM_NN, length(NN_NETWORKSIZE)), ...
  'net_loss',   cell(NUM_NN, length(NN_NETWORKSIZE)), ...
  'prob_cor',   cell(NUM_NN, length(NN_NETWORKSIZE)), ...
  'mitral_tpr', cell(NUM_NN, length(NN_NETWORKSIZE)), ...
  'mitral_fpr', cell(NUM_NN, length(NN_NETWORKSIZE)), ...
  'mitral_tnr', cell(NUM_NN, length(NN_NETWORKSIZE)), ...
  'mitral_fnr', cell(NUM_NN, length(NN_NETWORKSIZE)));

% Main loop for training NN's
% Use internal sampling to divide datasets
for n = NUM_NN:-1:1
  aux.progressbar(NUM_NN - n + 1, NUM_NN);
  for i = 1:length(NN_NETWORKSIZE)
    % Initialize a neural network
    sto_inet = patternnet(NN_NETWORKSIZE{i});
    sto_inet.trainParam.showWindow = false;
    sto_inet.divideParam.trainRatio = .9;
    sto_inet.divideParam.valRatio= 0;
    sto_inet.divideParam.testRatio= .1;
    % Train a network with weighting losses
    resper(n, i).title = num2str(NN_NETWORKSIZE{i});
    [resper(n, i).net, resper(n, i).trainRes] = train( ...
      sto_inet, aliOB.data.temp.prjReg', sto_temCls', [], [], ...
      aliOB.data.temp.nBrc ./ sqrt(sum(sto_temCls, 1)'));
  end
end

% Store all tried network architectures
resper2 = struct( ...
  'net',        cell(NUM_NN, length(NN_NETWORKSIZE)), ...
  'trainRes',   cell(NUM_NN, length(NN_NETWORKSIZE)), ...
  'prob_bar',   cell(NUM_NN, length(NN_NETWORKSIZE)), ...
  'class_id',   cell(NUM_NN, length(NN_NETWORKSIZE)), ...
  'net_loss',   cell(NUM_NN, length(NN_NETWORKSIZE)), ...
  'prob_cor',   cell(NUM_NN, length(NN_NETWORKSIZE)), ...
  'mitral_tpr', cell(NUM_NN, length(NN_NETWORKSIZE)), ...
  'mitral_fpr', cell(NUM_NN, length(NN_NETWORKSIZE)), ...
  'mitral_tnr', cell(NUM_NN, length(NN_NETWORKSIZE)), ...
  'mitral_fnr', cell(NUM_NN, length(NN_NETWORKSIZE)));

% Main loop for training NN's
% Use internal sampling to divide datasets
for n = NUM_NN:-1:1
  aux.progressbar(NUM_NN - n + 1, NUM_NN);
  for i = 1:length(NN_NETWORKSIZE)
    % Initialize a neural network
    sto_inet = patternnet(NN_NETWORKSIZE{i});
    sto_inet.trainParam.showWindow = false;
    sto_inet.divideParam.trainRatio = .8;
    sto_inet.divideParam.valRatio= 0;
    sto_inet.divideParam.testRatio= .2;
    % Train a network with weighting losses
    resper2(n, i).title = num2str(NN_NETWORKSIZE{i});
    [resper2(n, i).net, resper2(n, i).trainRes] = train( ...
      sto_inet, aliOB.data.temp.prjReg', sto_temCls', [], [], ...
      aliOB.data.temp.nBrc ./ sqrt(sum(sto_temCls, 1)'));
  end
end

% Release my seed
rng('shuffle');

%% Extract statistics from sets
TT = {'test', 'train'};
% Get data for the chosen network as well
sto = aliOB.getDataset('Merged');
barseqMat = aliOB.data.temp.prjReg;
barseqProb = aliOB_classifier.net(barseqMat')';
[~, barseq_tIdAct] = max(barseqProb, [], 2);
barseq_thr_tIdAct = barseqProb(:, sto_mId) >= ELIM_AMNT;
mapseqMat = sto.prjRegSum;
mapseqProb = aliOB_classifier.net(mapseqMat')';
[~, mapseq_tIdAct] = max(mapseqProb, [], 2);
mapseq_thr_tIdAct = mapseqProb(:, sto_mId) >= ELIM_AMNT;

for n = 1:NUM_NN
  for i = 1:length(NN_NETWORKSIZE)
    for tt = 1:length(TT)
      % Do for both training and testing sets
      T = TT{tt};
      % TESTING RUN STATISTICS
      sto_thisInd = resper(n, i).trainRes.([T, 'Ind']);
      sto_tMat = aliOB.data.temp.prjReg(sto_thisInd, :);
      sto_tProb = resper(n, i).net(sto_tMat')';
      sto_tIdAct = aliOB.data.temp.realId(sto_thisInd);
      sto_thr_tIdAct = sto_tIdAct == sto_mId;
      % Do identification using both methods
      [~, sto_tIdCls] = max(sto_tProb, [], 2);
      sto_thr_tIdCls = sto_tProb(:, sto_mId) >= ELIM_AMNT;
      % THRESHOLDING STATISTICS
      sto_mtp = sum(( sto_thr_tIdCls) & ( sto_thr_tIdAct));
      sto_mfp = sum(( sto_thr_tIdCls) & (~sto_thr_tIdAct));
      sto_mtn = sum((~sto_thr_tIdCls) & (~sto_thr_tIdAct));
      sto_mfn = sum((~sto_thr_tIdCls) & ( sto_thr_tIdAct));
      sto_mp = sto_mtp + sto_mfn;
      sto_mn = sto_mtn + sto_mfp;
      if sto_mp == 0
        sto_mp = nan;
      end
      if sto_mn == 0
        sto_mn = nan;
      end
      sto_mtpr = sto_mtp / sto_mp;
      sto_mfnr = sto_mfn / sto_mp;
      sto_mtnr = sto_mtn / sto_mn;
      sto_mfpr = sto_mfp / sto_mn;
      resper(n, i).(['mitral_', T, '_tp']) =  sto_mtp;
      resper(n, i).(['mitral_', T, '_fn']) =  sto_mfn;
      resper(n, i).(['mitral_', T, '_tn']) =  sto_mtn;
      resper(n, i).(['mitral_', T, '_fp']) =  sto_mfp;
      resper(n, i).(['mitral_', T, '_tpr']) = sto_mtpr;
      resper(n, i).(['mitral_', T, '_fnr']) = sto_mfnr;
      resper(n, i).(['mitral_', T, '_tnr']) = sto_mtnr;
      resper(n, i).(['mitral_', T, '_fpr']) = sto_mfpr;
      resper(n, i).(['mitral_', T, '_conf']) = ...
        confusionmat(sto_thr_tIdAct, sto_thr_tIdCls);
      % Each class statistics STATISTICS
      sto_tp = sum((sto_tIdCls == 1:3) & (sto_tIdAct == 1:3), 1);
      sto_fp = sum((sto_tIdCls == 1:3) & (sto_tIdAct ~= 1:3), 1);
      sto_tn = sum((sto_tIdCls ~= 1:3) & (sto_tIdAct ~= 1:3), 1);
      sto_fn = sum((sto_tIdCls ~= 1:3) & (sto_tIdAct == 1:3), 1);
      sto_p = sto_tp + sto_fn;
      sto_n = sto_tn + sto_fp;
      sto_p(sto_p(:) == 0) = nan;
      sto_n(sto_n(:) == 0) = nan;
      sto_tpr = sto_tp ./ sto_p;
      sto_fnr = sto_fn ./ sto_p;
      sto_tnr = sto_tn ./ sto_n;
      sto_fpr = sto_fp ./ sto_n;
      resper(n, i).([T, '_tp']) =  sto_tp;
      resper(n, i).([T, '_fn']) =  sto_fn;
      resper(n, i).([T, '_tn']) =  sto_tn;
      resper(n, i).([T, '_fp']) =  sto_fp;
      resper(n, i).([T, '_tpr']) = sto_tpr;
      resper(n, i).([T, '_fnr']) = sto_fnr;
      resper(n, i).([T, '_tnr']) = sto_tnr;
      resper(n, i).([T, '_fpr']) = sto_fpr;
      resper(n, i).([T, '_conf']) = ...
        confusionmat(sto_tIdAct, sto_tIdCls, 'ORDER', 1:3);
    end
    % Compare with the chosen network;
    sto_tProb = resper(n, i).net(barseqMat')';
    [~, sto_barseq_tIdCls] = max(sto_tProb, [], 2);
    sto_barseq_thr_tIdCls = sto_tProb(:, sto_mId) >= ELIM_AMNT;
    resper(n, i).('barseq_score_corr') =  corr(barseqProb, sto_tProb);
    sto_tProb = resper(n, i).net(mapseqMat')';
    [~, sto_mapseq_tIdCls] = max(sto_tProb, [], 2);
    sto_mapseq_thr_tIdCls = sto_tProb(:, sto_mId) >= ELIM_AMNT;
    resper(n, i).('mapseq_score_corr') =  corr(mapseqProb, sto_tProb);
    % Barseq, thresholding
    sto_mtp = sum(( sto_barseq_thr_tIdCls) & ( barseq_thr_tIdAct));
    sto_mfp = sum(( sto_barseq_thr_tIdCls) & (~barseq_thr_tIdAct));
    sto_mtn = sum((~sto_barseq_thr_tIdCls) & (~barseq_thr_tIdAct));
    sto_mfn = sum((~sto_barseq_thr_tIdCls) & ( barseq_thr_tIdAct));
    sto_mp = sto_mtp + sto_mfn;
    sto_mn = sto_mtn + sto_mfp;
    if sto_mp == 0
      sto_mp = nan;
    end
    if sto_mn == 0
      sto_mn = nan;
    end
    sto_mtpr = sto_mtp / sto_mp;
    sto_mfnr = sto_mfn / sto_mp;
    sto_mtnr = sto_mtn / sto_mn;
    sto_mfpr = sto_mfp / sto_mn;
    resper(n, i).('barseq_mitral_tp') =  sto_mtp;
    resper(n, i).('barseq_mitral_fn') =  sto_mfn;
    resper(n, i).('barseq_mitral_tn') =  sto_mtn;
    resper(n, i).('barseq_mitral_fp') =  sto_mfp;
    resper(n, i).('barseq_mitral_tpr') = sto_mtpr;
    resper(n, i).('barseq_mitral_fnr') = sto_mfnr;
    resper(n, i).('barseq_mitral_tnr') = sto_mtnr;
    resper(n, i).('barseq_mitral_fpr') = sto_mfpr;
    resper(n, i).('barseq_mitral_conf') = ...
      confusionmat(barseq_thr_tIdAct, sto_barseq_thr_tIdCls);
    resper(n, i).('barseq_mitral_corr') = ...
      corr(barseq_thr_tIdAct, sto_barseq_thr_tIdCls);
    % Mapseq, thresholding
    sto_mtp = sum(( sto_mapseq_thr_tIdCls) & ( mapseq_thr_tIdAct));
    sto_mfp = sum(( sto_mapseq_thr_tIdCls) & (~mapseq_thr_tIdAct));
    sto_mtn = sum((~sto_mapseq_thr_tIdCls) & (~mapseq_thr_tIdAct));
    sto_mfn = sum((~sto_mapseq_thr_tIdCls) & ( mapseq_thr_tIdAct));
    sto_mp = sto_mtp + sto_mfn;
    sto_mn = sto_mtn + sto_mfp;
    if sto_mp == 0
      sto_mp = nan;
    end
    if sto_mn == 0
      sto_mn = nan;
    end
    sto_mtpr = sto_mtp / sto_mp;
    sto_mfnr = sto_mfn / sto_mp;
    sto_mtnr = sto_mtn / sto_mn;
    sto_mfpr = sto_mfp / sto_mn;
    resper(n, i).('mapseq_mitral_tp') =  sto_mtp;
    resper(n, i).('mapseq_mitral_fn') =  sto_mfn;
    resper(n, i).('mapseq_mitral_tn') =  sto_mtn;
    resper(n, i).('mapseq_mitral_fp') =  sto_mfp;
    resper(n, i).('mapseq_mitral_tpr') = sto_mtpr;
    resper(n, i).('mapseq_mitral_fnr') = sto_mfnr;
    resper(n, i).('mapseq_mitral_tnr') = sto_mtnr;
    resper(n, i).('mapseq_mitral_fpr') = sto_mfpr;
    resper(n, i).('mapseq_mitral_conf') = ...
      confusionmat(mapseq_thr_tIdAct, sto_mapseq_thr_tIdCls);
    resper(n, i).('mapseq_mitral_corr') = ...
      corr(mapseq_thr_tIdAct, sto_mapseq_thr_tIdCls);
    % Barseq, max
    sto_mtp = sum((sto_barseq_tIdCls == 1:3) & (barseq_tIdAct == 1:3), 1);
    sto_mfp = sum((sto_barseq_tIdCls == 1:3) & (barseq_tIdAct ~= 1:3), 1);
    sto_mtn = sum((sto_barseq_tIdCls ~= 1:3) & (barseq_tIdAct ~= 1:3), 1);
    sto_mfn = sum((sto_barseq_tIdCls ~= 1:3) & (barseq_tIdAct == 1:3), 1);
    sto_mp = sto_mtp + sto_mfn;
    sto_mn = sto_mtn + sto_mfp;
    if sto_mp == 0
      sto_mp = nan;
    end
    if sto_mn == 0
      sto_mn = nan;
    end
    sto_mtpr = sto_mtp ./ sto_mp;
    sto_mfnr = sto_mfn ./ sto_mp;
    sto_mtnr = sto_mtn ./ sto_mn;
    sto_mfpr = sto_mfp ./ sto_mn;
    resper(n, i).('barseq_tp') =  sto_mtp;
    resper(n, i).('barseq_fn') =  sto_mfn;
    resper(n, i).('barseq_tn') =  sto_mtn;
    resper(n, i).('barseq_fp') =  sto_mfp;
    resper(n, i).('barseq_tpr') = sto_mtpr;
    resper(n, i).('barseq_fnr') = sto_mfnr;
    resper(n, i).('barseq_tnr') = sto_mtnr;
    resper(n, i).('barseq_fpr') = sto_mfpr;
    resper(n, i).('barseq_conf') = ...
      confusionmat(barseq_tIdAct, sto_barseq_tIdCls);
    resper(n, i).('barseq_corr') = ...
      corr(barseq_tIdAct == 1:3, sto_barseq_tIdCls == 1:3);
    % Mapseq, max
    sto_mtp = sum((sto_mapseq_tIdCls == 1:3) & (mapseq_tIdAct == 1:3), 1);
    sto_mfp = sum((sto_mapseq_tIdCls == 1:3) & (mapseq_tIdAct ~= 1:3), 1);
    sto_mtn = sum((sto_mapseq_tIdCls ~= 1:3) & (mapseq_tIdAct ~= 1:3), 1);
    sto_mfn = sum((sto_mapseq_tIdCls ~= 1:3) & (mapseq_tIdAct == 1:3), 1);
    sto_mp = sto_mtp + sto_mfn;
    sto_mn = sto_mtn + sto_mfp;
    if sto_mp == 0
      sto_mp = nan;
    end
    if sto_mn == 0
      sto_mn = nan;
    end
    sto_mtpr = sto_mtp ./ sto_mp;
    sto_mfnr = sto_mfn ./ sto_mp;
    sto_mtnr = sto_mtn ./ sto_mn;
    sto_mfpr = sto_mfp ./ sto_mn;
    resper(n, i).('mapseq_tp') =  sto_mtp;
    resper(n, i).('mapseq_fn') =  sto_mfn;
    resper(n, i).('mapseq_tn') =  sto_mtn;
    resper(n, i).('mapseq_fp') =  sto_mfp;
    resper(n, i).('mapseq_tpr') = sto_mtpr;
    resper(n, i).('mapseq_fnr') = sto_mfnr;
    resper(n, i).('mapseq_tnr') = sto_mtnr;
    resper(n, i).('mapseq_fpr') = sto_mfpr;
    resper(n, i).('mapseq_conf') = ...
      confusionmat(mapseq_tIdAct, sto_mapseq_tIdCls);
    resper(n, i).('mapseq_corr') = ...
      corr(mapseq_tIdAct == 1:3, sto_mapseq_tIdCls == 1:3);
  end
end

for n = 1:NUM_NN
  for i = 1:length(NN_NETWORKSIZE)
    for tt = 1:length(TT)
      % Do for both training and testing sets
      T = TT{tt};
      % TESTING RUN STATISTICS
      sto_thisInd = resper2(n, i).trainRes.([T, 'Ind']);
      sto_tMat = aliOB.data.temp.prjReg(sto_thisInd, :);
      sto_tProb = resper2(n, i).net(sto_tMat')';
      sto_tIdAct = aliOB.data.temp.realId(sto_thisInd);
      sto_thr_tIdAct = sto_tIdAct == sto_mId;
      % Do identification using both methods
      [~, sto_tIdCls] = max(sto_tProb, [], 2);
      sto_thr_tIdCls = false(size(sto_tProb, 1), 1);
      sto_thr_tIdCls(sto_tProb(:, sto_mId) >= ELIM_AMNT) = true;
      % THRESHOLDING STATISTICS
      sto_mtp = sum(( sto_thr_tIdCls) & ( sto_thr_tIdAct));
      sto_mfp = sum(( sto_thr_tIdCls) & (~sto_thr_tIdAct));
      sto_mtn = sum((~sto_thr_tIdCls) & (~sto_thr_tIdAct));
      sto_mfn = sum((~sto_thr_tIdCls) & ( sto_thr_tIdAct));
      sto_mp = sto_mtp + sto_mfn;
      sto_mn = sto_mtn + sto_mfp;
      if sto_mp == 0
        sto_mp = nan;
      end
      if sto_mn == 0
        sto_mn = nan;
      end
      sto_mtpr = sto_mtp / sto_mp;
      sto_mfnr = sto_mfn / sto_mp;
      sto_mtnr = sto_mtn / sto_mn;
      sto_mfpr = sto_mfp / sto_mn;
      resper2(n, i).(['mitral_', T, '_tp']) =  sto_mtp;
      resper2(n, i).(['mitral_', T, '_fn']) =  sto_mfn;
      resper2(n, i).(['mitral_', T, '_tn']) =  sto_mtn;
      resper2(n, i).(['mitral_', T, '_fp']) =  sto_mfp;
      resper2(n, i).(['mitral_', T, '_tpr']) = sto_mtpr;
      resper2(n, i).(['mitral_', T, '_fnr']) = sto_mfnr;
      resper2(n, i).(['mitral_', T, '_tnr']) = sto_mtnr;
      resper2(n, i).(['mitral_', T, '_fpr']) = sto_mfpr;
      resper2(n, i).(['mitral_', T, '_conf']) = ...
        confusionmat(sto_thr_tIdAct, sto_thr_tIdCls);
      % Each class statistics STATISTICS
      sto_tp = sum(( (sto_tIdCls == 1:3)) & ( (sto_tIdAct == 1:3)), 1);
      sto_fp = sum(( (sto_tIdCls == 1:3)) & (~(sto_tIdAct == 1:3)), 1);
      sto_tn = sum((~(sto_tIdCls == 1:3)) & (~(sto_tIdAct == 1:3)), 1);
      sto_fn = sum((~(sto_tIdCls == 1:3)) & ( (sto_tIdAct == 1:3)), 1);
      sto_p = sto_tp + sto_fn;
      sto_n = sto_tn + sto_fp;
      sto_p(sto_p(:) == 0) = nan;
      sto_n(sto_n(:) == 0) = nan;
      sto_tpr = sto_tp ./ sto_p;
      sto_fnr = sto_fn ./ sto_p;
      sto_tnr = sto_tn ./ sto_n;
      sto_fpr = sto_fp ./ sto_n;
      resper2(n, i).([T, '_tp']) =  sto_tp;
      resper2(n, i).([T, '_fn']) =  sto_fn;
      resper2(n, i).([T, '_tn']) =  sto_tn;
      resper2(n, i).([T, '_fp']) =  sto_fp;
      resper2(n, i).([T, '_tpr']) = sto_tpr;
      resper2(n, i).([T, '_fnr']) = sto_fnr;
      resper2(n, i).([T, '_tnr']) = sto_tnr;
      resper2(n, i).([T, '_fpr']) = sto_fpr;
      resper2(n, i).([T, '_conf']) = ...
        confusionmat(sto_tIdAct, sto_tIdCls, 'ORDER', 1:3);
    end
    % Compare with the chosen network;
    sto_tProb = resper2(n, i).net(barseqMat')';
    [~, sto_barseq_tIdCls] = max(sto_tProb, [], 2);
    sto_barseq_thr_tIdCls = sto_tProb(:, sto_mId) >= ELIM_AMNT;
    resper2(n, i).('barseq_score_corr') =  corr(barseqProb, sto_tProb);
    sto_tProb = resper2(n, i).net(mapseqMat')';
    [~, sto_mapseq_tIdCls] = max(sto_tProb, [], 2);
    sto_mapseq_thr_tIdCls = sto_tProb(:, sto_mId) >= ELIM_AMNT;
    resper2(n, i).('mapseq_score_corr') =  corr(mapseqProb, sto_tProb);
    % Barseq, thresholding
    sto_mtp = sum(( sto_barseq_thr_tIdCls) & ( barseq_thr_tIdAct));
    sto_mfp = sum(( sto_barseq_thr_tIdCls) & (~barseq_thr_tIdAct));
    sto_mtn = sum((~sto_barseq_thr_tIdCls) & (~barseq_thr_tIdAct));
    sto_mfn = sum((~sto_barseq_thr_tIdCls) & ( barseq_thr_tIdAct));
    sto_mp = sto_mtp + sto_mfn;
    sto_mn = sto_mtn + sto_mfp;
    if sto_mp == 0
      sto_mp = nan;
    end
    if sto_mn == 0
      sto_mn = nan;
    end
    sto_mtpr = sto_mtp / sto_mp;
    sto_mfnr = sto_mfn / sto_mp;
    sto_mtnr = sto_mtn / sto_mn;
    sto_mfpr = sto_mfp / sto_mn;
    resper2(n, i).('barseq_mitral_tp') =  sto_mtp;
    resper2(n, i).('barseq_mitral_fn') =  sto_mfn;
    resper2(n, i).('barseq_mitral_tn') =  sto_mtn;
    resper2(n, i).('barseq_mitral_fp') =  sto_mfp;
    resper2(n, i).('barseq_mitral_tpr') = sto_mtpr;
    resper2(n, i).('barseq_mitral_fnr') = sto_mfnr;
    resper2(n, i).('barseq_mitral_tnr') = sto_mtnr;
    resper2(n, i).('barseq_mitral_fpr') = sto_mfpr;
    resper2(n, i).('barseq_mitral_conf') = ...
      confusionmat(barseq_thr_tIdAct, sto_barseq_thr_tIdCls);
    resper2(n, i).('barseq_mitral_corr') = ...
      corr(barseq_thr_tIdAct, sto_barseq_thr_tIdCls);
    % Mapseq, thresholding
    sto_mtp = sum(( sto_mapseq_thr_tIdCls) & ( mapseq_thr_tIdAct));
    sto_mfp = sum(( sto_mapseq_thr_tIdCls) & (~mapseq_thr_tIdAct));
    sto_mtn = sum((~sto_mapseq_thr_tIdCls) & (~mapseq_thr_tIdAct));
    sto_mfn = sum((~sto_mapseq_thr_tIdCls) & ( mapseq_thr_tIdAct));
    sto_mp = sto_mtp + sto_mfn;
    sto_mn = sto_mtn + sto_mfp;
    if sto_mp == 0
      sto_mp = nan;
    end
    if sto_mn == 0
      sto_mn = nan;
    end
    sto_mtpr = sto_mtp / sto_mp;
    sto_mfnr = sto_mfn / sto_mp;
    sto_mtnr = sto_mtn / sto_mn;
    sto_mfpr = sto_mfp / sto_mn;
    resper2(n, i).('mapseq_mitral_tp') =  sto_mtp;
    resper2(n, i).('mapseq_mitral_fn') =  sto_mfn;
    resper2(n, i).('mapseq_mitral_tn') =  sto_mtn;
    resper2(n, i).('mapseq_mitral_fp') =  sto_mfp;
    resper2(n, i).('mapseq_mitral_tpr') = sto_mtpr;
    resper2(n, i).('mapseq_mitral_fnr') = sto_mfnr;
    resper2(n, i).('mapseq_mitral_tnr') = sto_mtnr;
    resper2(n, i).('mapseq_mitral_fpr') = sto_mfpr;
    resper2(n, i).('mapseq_mitral_conf') = ...
      confusionmat(mapseq_thr_tIdAct, sto_mapseq_thr_tIdCls);
    resper2(n, i).('mapseq_mitral_corr') = ...
      corr(mapseq_thr_tIdAct, sto_mapseq_thr_tIdCls);
    % Barseq, max
    sto_mtp = sum((sto_barseq_tIdCls == 1:3) & (barseq_tIdAct == 1:3), 1);
    sto_mfp = sum((sto_barseq_tIdCls == 1:3) & (barseq_tIdAct ~= 1:3), 1);
    sto_mtn = sum((sto_barseq_tIdCls ~= 1:3) & (barseq_tIdAct ~= 1:3), 1);
    sto_mfn = sum((sto_barseq_tIdCls ~= 1:3) & (barseq_tIdAct == 1:3), 1);
    sto_mp = sto_mtp + sto_mfn;
    sto_mn = sto_mtn + sto_mfp;
    if sto_mp == 0
      sto_mp = nan;
    end
    if sto_mn == 0
      sto_mn = nan;
    end
    sto_mtpr = sto_mtp ./ sto_mp;
    sto_mfnr = sto_mfn ./ sto_mp;
    sto_mtnr = sto_mtn ./ sto_mn;
    sto_mfpr = sto_mfp ./ sto_mn;
    resper2(n, i).('barseq_tp') =  sto_mtp;
    resper2(n, i).('barseq_fn') =  sto_mfn;
    resper2(n, i).('barseq_tn') =  sto_mtn;
    resper2(n, i).('barseq_fp') =  sto_mfp;
    resper2(n, i).('barseq_tpr') = sto_mtpr;
    resper2(n, i).('barseq_fnr') = sto_mfnr;
    resper2(n, i).('barseq_tnr') = sto_mtnr;
    resper2(n, i).('barseq_fpr') = sto_mfpr;
    resper2(n, i).('barseq_conf') = ...
      confusionmat(barseq_tIdAct, sto_barseq_tIdCls);
    resper2(n, i).('barseq_corr') = ...
      corr(barseq_tIdAct == 1:3, sto_barseq_tIdCls == 1:3);
    % Mapseq, max
    sto_mtp = sum((sto_mapseq_tIdCls == 1:3) & (mapseq_tIdAct == 1:3), 1);
    sto_mfp = sum((sto_mapseq_tIdCls == 1:3) & (mapseq_tIdAct ~= 1:3), 1);
    sto_mtn = sum((sto_mapseq_tIdCls ~= 1:3) & (mapseq_tIdAct ~= 1:3), 1);
    sto_mfn = sum((sto_mapseq_tIdCls ~= 1:3) & (mapseq_tIdAct == 1:3), 1);
    sto_mp = sto_mtp + sto_mfn;
    sto_mn = sto_mtn + sto_mfp;
    if sto_mp == 0
      sto_mp = nan;
    end
    if sto_mn == 0
      sto_mn = nan;
    end
    sto_mtpr = sto_mtp ./ sto_mp;
    sto_mfnr = sto_mfn ./ sto_mp;
    sto_mtnr = sto_mtn ./ sto_mn;
    sto_mfpr = sto_mfp ./ sto_mn;
    resper2(n, i).('mapseq_tp') =  sto_mtp;
    resper2(n, i).('mapseq_fn') =  sto_mfn;
    resper2(n, i).('mapseq_tn') =  sto_mtn;
    resper2(n, i).('mapseq_fp') =  sto_mfp;
    resper2(n, i).('mapseq_tpr') = sto_mtpr;
    resper2(n, i).('mapseq_fnr') = sto_mfnr;
    resper2(n, i).('mapseq_tnr') = sto_mtnr;
    resper2(n, i).('mapseq_fpr') = sto_mfpr;
    resper2(n, i).('mapseq_conf') = ...
      confusionmat(mapseq_tIdAct, sto_mapseq_tIdCls);
    resper2(n, i).('mapseq_corr') = ...
      corr(mapseq_tIdAct == 1:3, sto_mapseq_tIdCls == 1:3);
  end
end

for n = 1:NUM_NN
  for i = 1:length(NN_NETWORKSIZE)
    % TESTING RUN STATISTICS
    sto_tMat = aliOB.data.temp.prjReg;
    sto_tProb = res(n, i).net(sto_tMat')';
    sto_tIdAct = aliOB.data.temp.realId;
    sto_thr_tIdAct = sto_tIdAct == sto_mId;
    % Do identification using both methods
    [~, sto_tIdCls] = max(sto_tProb, [], 2);
    sto_thr_tIdCls = false(size(sto_tProb, 1), 1);
    sto_thr_tIdCls(sto_tProb(:, sto_mId) >= ELIM_AMNT) = true;
    % THRESHOLDING STATISTICS
    sto_mtp = sum(( sto_thr_tIdCls) & ( sto_thr_tIdAct));
    sto_mfp = sum(( sto_thr_tIdCls) & (~sto_thr_tIdAct));
    sto_mtn = sum((~sto_thr_tIdCls) & (~sto_thr_tIdAct));
    sto_mfn = sum((~sto_thr_tIdCls) & ( sto_thr_tIdAct));
    sto_mp = sto_mtp + sto_mfn;
    sto_mn = sto_mtn + sto_mfp;
    if sto_mp == 0
      sto_mp = nan;
    end
    if sto_mn == 0
      sto_mn = nan;
    end
    sto_mtpr = sto_mtp / sto_mp;
    sto_mfnr = sto_mfn / sto_mp;
    sto_mtnr = sto_mtn / sto_mn;
    sto_mfpr = sto_mfp / sto_mn;
    res(n, i).(['mitral_', T, '_tp']) =  sto_mtp;
    res(n, i).(['mitral_', T, '_fn']) =  sto_mfn;
    res(n, i).(['mitral_', T, '_tn']) =  sto_mtn;
    res(n, i).(['mitral_', T, '_fp']) =  sto_mfp;
    res(n, i).(['mitral_', T, '_tpr']) = sto_mtpr;
    res(n, i).(['mitral_', T, '_fnr']) = sto_mfnr;
    res(n, i).(['mitral_', T, '_tnr']) = sto_mtnr;
    res(n, i).(['mitral_', T, '_fpr']) = sto_mfpr;
    res(n, i).(['mitral_', T, '_conf']) = ...
      confusionmat(sto_thr_tIdAct, sto_thr_tIdCls);
    % Each class statistics STATISTICS
    sto_tp = sum(( (sto_tIdCls == 1:3)) & ( (sto_tIdAct == 1:3)), 1);
    sto_fp = sum(( (sto_tIdCls == 1:3)) & (~(sto_tIdAct == 1:3)), 1);
    sto_tn = sum((~(sto_tIdCls == 1:3)) & (~(sto_tIdAct == 1:3)), 1);
    sto_fn = sum((~(sto_tIdCls == 1:3)) & ( (sto_tIdAct == 1:3)), 1);
    sto_p = sto_tp + sto_fn;
    sto_n = sto_tn + sto_fp;
    sto_p(sto_p(:) == 0) = nan;
    sto_n(sto_n(:) == 0) = nan;
    sto_tpr = sto_tp ./ sto_p;
    sto_fnr = sto_fn ./ sto_p;
    sto_tnr = sto_tn ./ sto_n;
    sto_fpr = sto_fp ./ sto_n;
    res(n, i).([T, '_tp']) =  sto_tp;
    res(n, i).([T, '_fn']) =  sto_fn;
    res(n, i).([T, '_tn']) =  sto_tn;
    res(n, i).([T, '_fp']) =  sto_fp;
    res(n, i).([T, '_tpr']) = sto_tpr;
    res(n, i).([T, '_fnr']) = sto_fnr;
    res(n, i).([T, '_tnr']) = sto_tnr;
    res(n, i).([T, '_fpr']) = sto_fpr;
    res(n, i).([T, '_conf']) = ...
      confusionmat(sto_tIdAct, sto_tIdCls, 'ORDER', 1:3);
    res(n, i).('mapseq_mitral_corr') = ...
      corr(mapseq_thr_tIdAct, sto_mapseq_thr_tIdCls);
    % Compare with the chosen network;
    sto_tProb = res(n, i).net(barseqMat')';
    [~, sto_barseq_tIdCls] = max(sto_tProb, [], 2);
    sto_barseq_thr_tIdCls = sto_tProb(:, sto_mId) >= ELIM_AMNT;
    res(n, i).('barseq_score_corr') =  corr(barseqProb, sto_tProb);
    sto_tProb = res(n, i).net(mapseqMat')';
    [~, sto_mapseq_tIdCls] = max(sto_tProb, [], 2);
    sto_mapseq_thr_tIdCls = sto_tProb(:, sto_mId) >= ELIM_AMNT;
    res(n, i).('mapseq_score_corr') =  corr(mapseqProb, sto_tProb);
    % Barseq, thresholding
    sto_mtp = sum(( sto_barseq_thr_tIdCls) & ( barseq_thr_tIdAct));
    sto_mfp = sum(( sto_barseq_thr_tIdCls) & (~barseq_thr_tIdAct));
    sto_mtn = sum((~sto_barseq_thr_tIdCls) & (~barseq_thr_tIdAct));
    sto_mfn = sum((~sto_barseq_thr_tIdCls) & ( barseq_thr_tIdAct));
    sto_mp = sto_mtp + sto_mfn;
    sto_mn = sto_mtn + sto_mfp;
    if sto_mp == 0
      sto_mp = nan;
    end
    if sto_mn == 0
      sto_mn = nan;
    end
    sto_mtpr = sto_mtp / sto_mp;
    sto_mfnr = sto_mfn / sto_mp;
    sto_mtnr = sto_mtn / sto_mn;
    sto_mfpr = sto_mfp / sto_mn;
    res(n, i).('barseq_mitral_tp') =  sto_mtp;
    res(n, i).('barseq_mitral_fn') =  sto_mfn;
    res(n, i).('barseq_mitral_tn') =  sto_mtn;
    res(n, i).('barseq_mitral_fp') =  sto_mfp;
    res(n, i).('barseq_mitral_tpr') = sto_mtpr;
    res(n, i).('barseq_mitral_fnr') = sto_mfnr;
    res(n, i).('barseq_mitral_tnr') = sto_mtnr;
    res(n, i).('barseq_mitral_fpr') = sto_mfpr;
    res(n, i).('barseq_mitral_conf') = ...
      confusionmat(barseq_thr_tIdAct, sto_barseq_thr_tIdCls);
    res(n, i).('barseq_mitral_corr') = ...
      corr(barseq_thr_tIdAct, sto_barseq_thr_tIdCls);
    % Mapseq, thresholding
    sto_mtp = sum(( sto_mapseq_thr_tIdCls) & ( mapseq_thr_tIdAct));
    sto_mfp = sum(( sto_mapseq_thr_tIdCls) & (~mapseq_thr_tIdAct));
    sto_mtn = sum((~sto_mapseq_thr_tIdCls) & (~mapseq_thr_tIdAct));
    sto_mfn = sum((~sto_mapseq_thr_tIdCls) & ( mapseq_thr_tIdAct));
    sto_mp = sto_mtp + sto_mfn;
    sto_mn = sto_mtn + sto_mfp;
    if sto_mp == 0
      sto_mp = nan;
    end
    if sto_mn == 0
      sto_mn = nan;
    end
    sto_mtpr = sto_mtp / sto_mp;
    sto_mfnr = sto_mfn / sto_mp;
    sto_mtnr = sto_mtn / sto_mn;
    sto_mfpr = sto_mfp / sto_mn;
    res(n, i).('mapseq_mitral_tp') =  sto_mtp;
    res(n, i).('mapseq_mitral_fn') =  sto_mfn;
    res(n, i).('mapseq_mitral_tn') =  sto_mtn;
    res(n, i).('mapseq_mitral_fp') =  sto_mfp;
    res(n, i).('mapseq_mitral_tpr') = sto_mtpr;
    res(n, i).('mapseq_mitral_fnr') = sto_mfnr;
    res(n, i).('mapseq_mitral_tnr') = sto_mtnr;
    res(n, i).('mapseq_mitral_fpr') = sto_mfpr;
    res(n, i).('mapseq_mitral_conf') = ...
      confusionmat(mapseq_thr_tIdAct, sto_mapseq_thr_tIdCls);
    res(n, i).('mapseq_mitral_corr') = ...
      corr(mapseq_thr_tIdAct, sto_mapseq_thr_tIdCls);
    % Barseq, max
    sto_mtp = sum((sto_barseq_tIdCls == 1:3) & (barseq_tIdAct == 1:3), 1);
    sto_mfp = sum((sto_barseq_tIdCls == 1:3) & (barseq_tIdAct ~= 1:3), 1);
    sto_mtn = sum((sto_barseq_tIdCls ~= 1:3) & (barseq_tIdAct ~= 1:3), 1);
    sto_mfn = sum((sto_barseq_tIdCls ~= 1:3) & (barseq_tIdAct == 1:3), 1);
    sto_mp = sto_mtp + sto_mfn;
    sto_mn = sto_mtn + sto_mfp;
    if sto_mp == 0
      sto_mp = nan;
    end
    if sto_mn == 0
      sto_mn = nan;
    end
    sto_mtpr = sto_mtp ./ sto_mp;
    sto_mfnr = sto_mfn ./ sto_mp;
    sto_mtnr = sto_mtn ./ sto_mn;
    sto_mfpr = sto_mfp ./ sto_mn;
    res(n, i).('barseq_tp') =  sto_mtp;
    res(n, i).('barseq_fn') =  sto_mfn;
    res(n, i).('barseq_tn') =  sto_mtn;
    res(n, i).('barseq_fp') =  sto_mfp;
    res(n, i).('barseq_tpr') = sto_mtpr;
    res(n, i).('barseq_fnr') = sto_mfnr;
    res(n, i).('barseq_tnr') = sto_mtnr;
    res(n, i).('barseq_fpr') = sto_mfpr;
    res(n, i).('barseq_conf') = ...
      confusionmat(barseq_tIdAct, sto_barseq_tIdCls);
    res(n, i).('barseq_corr') = ...
      corr(barseq_tIdAct == 1:3, sto_barseq_tIdCls == 1:3);
    % Mapseq, max
    sto_mtp = sum((sto_mapseq_tIdCls == 1:3) & (mapseq_tIdAct == 1:3), 1);
    sto_mfp = sum((sto_mapseq_tIdCls == 1:3) & (mapseq_tIdAct ~= 1:3), 1);
    sto_mtn = sum((sto_mapseq_tIdCls ~= 1:3) & (mapseq_tIdAct ~= 1:3), 1);
    sto_mfn = sum((sto_mapseq_tIdCls ~= 1:3) & (mapseq_tIdAct == 1:3), 1);
    sto_mp = sto_mtp + sto_mfn;
    sto_mn = sto_mtn + sto_mfp;
    if sto_mp == 0
      sto_mp = nan;
    end
    if sto_mn == 0
      sto_mn = nan;
    end
    sto_mtpr = sto_mtp ./ sto_mp;
    sto_mfnr = sto_mfn ./ sto_mp;
    sto_mtnr = sto_mtn ./ sto_mn;
    sto_mfpr = sto_mfp ./ sto_mn;
    res(n, i).('mapseq_tp') =  sto_mtp;
    res(n, i).('mapseq_fn') =  sto_mfn;
    res(n, i).('mapseq_tn') =  sto_mtn;
    res(n, i).('mapseq_fp') =  sto_mfp;
    res(n, i).('mapseq_tpr') = sto_mtpr;
    res(n, i).('mapseq_fnr') = sto_mfnr;
    res(n, i).('mapseq_tnr') = sto_mtnr;
    res(n, i).('mapseq_fpr') = sto_mfpr;
    res(n, i).('mapseq_conf') = ...
      confusionmat(mapseq_tIdAct, sto_mapseq_tIdCls);
    res(n, i).('mapseq_corr') = ...
      corr(mapseq_tIdAct == 1:3, sto_mapseq_tIdCls == 1:3);
  end
end

%% Save data
nn_sets = res;
nn_sets_test = resper;
save('data/classifier.mat', ...
  'aliOB_classifier', 'aliOB_classifier_vers', 'REG_COR', ...
  'nn_sets', 'nn_sets_test', 'resper2');