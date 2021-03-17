function outset = obGetType(thisset, classifier, elim_thr, cell_type)
%OBGETMITRAL Summary of this function goes here

if ~exist('classifier', 'var')
  classifier = thisset.data.classifier;
end
if ~exist('elim_thr', 'var')
  elim_thr = .85;
end
if ~exist('cell_type', 'var')
  cell_type = 'Mitral';
end

% Get the mitral classifier
sto_m = find(strcmp(thisset.data.temp.names, cell_type), 1);

% Run the classifier
sto_score = classifier.net(thisset.prjRegSum')';
sto_ind = sto_score(:, sto_m) >= elim_thr;

% Create output set
outset = mapseqData;
outset.name = [thisset.name, ' (', cell_type, ')'];
outset.srcImg = thisset.srcImg(sto_ind, :);
outset.srcRegName = thisset.srcRegName;
outset.nSrcRegSli = thisset.nSrcRegSli;
outset.prjImg = thisset.prjImg(sto_ind, :);
outset.prjRegName = thisset.prjRegName;
outset.nPrjRegSli = thisset.nPrjRegSli;
outset.brcId = thisset.brcId(sto_ind, :);
outset.brcName = thisset.brcName;
outset.brcCor = thisset.brcCor(sto_ind, :);
outset.data.temp = thisset.data.temp;

% Put in this classifier
outset.data.classifier = classifier;
end

