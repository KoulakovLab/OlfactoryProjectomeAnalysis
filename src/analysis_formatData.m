% Script that loads the data in to one object array,
%   currently loads only for data gathered until October 2019
D = 1;

RM_THR = 5;

%% Olfactory bulb injections - Oct 2019
% load('data/obinj_2019_10.mat', 'rawData2');
% data_ob(length(rawData2), 1) = mapseqData;
% for d = 1:length(rawData2)
%     data_ob(d).name = strrep(rawData2(d).name, 'ss', '');
%     % Get the 'olfactory bulb' region and non-empty regions
%     obFlag = contains(rawData2(d).regName, 'OB');
%     exFlag = cellfun(@(x) ~isempty(x), rawData2(d).region);
%     % Load the images for the source region
%     sto = obFlag & exFlag;
%     data_ob(d).srcRegName = rawData2(d).regName(sto);
%     data_ob(d).nSrcRegSli = cellfun(@length, rawData2(d).region(sto));
%     data_ob(d).srcImg = cell2mat(cellfun(@(x) rawData2(d).matrix(:, x), ...
%         rawData2(d).region(sto), 'UniformOutput', 0)');
%     % Load the images for the projection region
%     sto = (~obFlag) & exFlag;
%     data_ob(d).prjRegName = rawData2(d).regName(sto);
%     data_ob(d).nPrjRegSli = cellfun(@length, rawData2(d).region(sto));
%     data_ob(d).prjImg = cell2mat(cellfun(@(x) rawData2(d).matrix(:, x), ...
%         rawData2(d).region(sto), 'UniformOutput', 0)');
% end
% 
% % The B86, B92, B109 and B111 have dorso-ventral stuff; arrange accordingly
% DV_DATA = {'B86', 'B92', 'B109', 'B111'};
% for d = 1:length(data_ob)
%     if any(strcmp(data_ob(d).name, DV_DATA))
%         % Find the source region with label OB
%         sto = find(strcmp(data_ob(d).srcRegName, 'OB'), 1);
%         % Dorsal barcodes are odd, and ventral barcodes are even
%         % Reorder the source image
%         data_ob(d).srcImg(:, data_ob(d).srcRegInd{sto}) = [ ...
%             data_ob(d).srcImg(:, data_ob(d).srcRegInd{sto}(1:2:end)), ...
%             data_ob(d).srcImg(:, data_ob(d).srcRegInd{sto}(2:2:end))];
%         % Relabel indices
%         data_ob(d).srcRegName = [ ...
%             data_ob(d).srcRegName(1:(sto-1)); ...
%             {'OB dorsal'; 'OB ventral'}; ...
%             data_ob(d).srcRegName((sto+1):end) ];
%         data_ob(d).nSrcRegSli = [ ...
%             data_ob(d).nSrcRegSli(1:(sto-1)); ...
%             .5 * data_ob(d).nSrcRegSli(sto) * [1; 1]; ...
%             data_ob(d).nSrcRegSli((sto+1):end) ];
%         % Do barcode identifiers by magnitude
%         data_ob(d).brcId = D + ( ...
%             sum(data_ob(d).getImgSrcReg('OB dorsal'), 2) < ...
%             sum(data_ob(d).getImgSrcReg('OB ventral'), 2) );
%         D = D + 2;
%     else
%         % Just identify by itself
%         data_ob(d).brcId = D * ones(data_ob(d).nBrc, 1);
%         D = D + 1;
%     end
% end

%% Olfactory bulb injections - Oct 2019
load('data/raw/pooleddata.mat', ...
  'alldata', 'alldatass', 'brainidx', 'somaslice', 'dvidx');
load('data/raw/obinj_clus_2020_05.mat', 'clustidxs', 'clustidxslabels');
load('data/raw/obinj_class_2020_05.mat', 'classidx', 'classlabel');
brainlabel = {'yc61'; 'yc65'; 'yc86'; 'yc92'; 'yc109'; 'yc111'};
regionlabel= {'AON'; 'APC'; 'PPC'; 'OT'; 'CoA'; 'lENT'};
data_ob(length(brainlabel), 1) = mapseqData;
ALLDATASS = [alldatass{:}];
for d = 1:length(brainlabel)
    data_ob(d).name = brainlabel{d};
    IND = (brainidx == d);
    % Olfactory bulb is empty; just do a full 1
    data_ob(d).srcRegName = {'OB'};
    data_ob(d).nSrcRegSli = ones(1);
    data_ob(d).srcImg = ones(sum(IND), 1);
    % Load the images for the projection region
    data_ob(d).prjRegName = regionlabel';
    data_ob(d).nPrjRegSli = cellfun(@(x) size(x, 2), alldatass);
    data_ob(d).prjImg = ALLDATASS(IND, :);
    data_ob(d).prjRegSum = alldata(IND, :);
    % Load extrenous info
    data_ob(d).data.brcDV = double(dvidx(IND, :));
    data_ob(d).data.brcSoma = somaslice(IND, :);
    % Create a barcode id mask
    data_ob(d).brcId = brainidx(IND);
    data_ob(d).brcName = cell(1, 1);
    data_ob(d).brcName{1} = brainlabel;
end

%% Olfactory bulb injections with BARSEQ- Feb 2020
load('data/raw/obbarseq_2020_02.mat', 'cbrainidx', 'cfiltproj', ...
    'cfiltprojss', 'cfiltroi', 'cfiltseq', 'cfiltslicename');
bname = {'YC107'; 'YC113'};
rname = {'AON'; 'APC'; 'PPC'; 'OT'; 'lENT'};%, 'Control'};
map = 'gtac';
data_bs(length(bname), 1) = mapseqData;
for d = 1:length(bname)
    data_bs(d).name = bname{d};
    % Get the index of the barcodes in this brain region
    sto = contains(cfiltslicename, bname{d});
    % Projection region info
    data_bs(d).prjRegName = rname;
    data_bs(d).nPrjRegSli = ones(length(rname), 1);
    data_bs(d).prjImg = cfiltproj(sto, 1:length(rname));
    % Source region info!
    % Convert get max pos names.
    posname = extractAfter(cfiltslicename(sto), 'Pos');
    % Pad with 0's
    ml = max(cellfun(@length, posname), [], 'all');
    posname = cellfun(@(x) ['Pos_', char('0' * ones(1, ...
        max(cellfun(@length, posname), [], 'all') - length(x))), x], ...
        posname, 'UniformOutput', false);
    % Convert and find unique id's
    data_bs(d).srcRegName = unique(posname);
    data_bs(d).nSrcRegSli = ones(length(data_bs(d).srcRegName), 1);
    % Create one-shot matrix for positions per barcode
    data_bs(d).srcImg = zeros(sum(sto), length(data_bs(d).srcRegName));
    for s = 1:length(data_bs(d).srcRegName)
        data_bs(d).srcImg(strcmp(posname, data_bs(d).srcRegName{s}), s) ...
            = 1;
    end
    % Barcode info
    data_bs(d).brcId = D * ones(data_bs(d).nBrc, 1);
    D = D + 1;
    % Additional info per barcode
    seqn = map(cfiltseq(sto, :));
    data_bs(d).data.brcSeq = arrayfun(@(x) seqn(x, :), ...
        (1:size(seqn, 1))', 'UniformOutput', 0);
    data_bs(d).data.brcXcor = cfiltroi(sto, 1);
    data_bs(d).data.brcYcor = cfiltroi(sto, 2);
    % Sort per max
    [val, sto] = max(data_bs(d).prjImg, [], 2);
    [~, ind] = sortrows([sto, val], 1:2, 'descend');
    data_bs(d).reorderBarcodes(ind);
end

% Load template
load('data/raw/obbarseq_2020_07.mat', ...
  'alldata', 'alldepths', 'allid', 'alltemplate');
load('data/raw/obbarseq_2020_02.mat', ...
  'template');
allid = allid(alltemplate);
alldata = alldata(alltemplate, :);
alldepths = alldepths(alltemplate, :);
for i = 1:3
  sto_ind = allid == i;
  template(i).data = alldata(sto_ind, :);
  template(i).pos = alldepths(sto_ind, :);
end

%% tSNE embedding for this data

% Create a tSNE embedding of the template + recorded data
fprintf('Running TSNE to put barcodes in space . . .\n');
rng(0);
sto_brcM = {data_ob.prjRegSum}';
sto_brcT = {template.data}';
sto_brc = cell2mat([sto_brcM; sto_brcT]);
% Get index
sto_numM = cellfun(@(x) size(x, 1), sto_brcM);
sto_cumM = cumsum([0; sto_numM]);
sto_numT = cellfun(@(x) size(x, 1), sto_brcT);
sto_cumT = sto_cumM(end) + cumsum([0; sto_numT]);
sto_idM = arrayfun(@(x, y) (x+1):y, ...
  sto_cumM(1:(end - 1)), sto_cumM(2:end), 'UniformOutput', false);
sto_idT = arrayfun(@(x, y) (x+1):y, ...
  sto_cumT(1:(end - 1)), sto_cumT(2:end), 'UniformOutput', false);
% Do tsne
sto_tsne = tsne(sto_brc, 'Perplexity', 50);
for i = 1:length(data_ob)
  data_ob(i).brcCor = sto_tsne(sto_idM{i}, :);
end
for i = 1:length(template)
  template(i).cor = sto_tsne(sto_idT{i}, :);
end
fprintf('\b Done!\n');
rng('shuffle');

%% Piriform cortex injections - May 2020
% XC 113,119,120 IG, has extra IG
% YC 123 has extra HY(hyperthalamus)
% YC 123 and XC 127 has extra CTX(cerebelum)
% XC 127 has extra H2O control (noted twice)
% XC 113, 119 has extra RT Neg Ctrl
% XC 127 has an extra slice in piriform
% The common regions across the data, use this ordering
regname = {...
    'OB', ...           1
    'OB caudal', ...    2
    'AON', ...          3
    'ORB', ...          4
    'TTD', ...          5
    'OT', ...           6
    'ACB', ...          7
    'CP', ...           8
    'APC caudal', ...   9
    'PPC caudal', ...  10
    'CoA', ...         11
    'MD', ...          12
    'Hip', ...         13
    'ENT'}; %          14
brainname = {'XC113', 'XC119', 'XC120', 'YC123', 'XC127'};

data_pc(5, 1) = mapseqData;
data_pcuf(5, 1) = mapseqData;

% XC113
%---Thresholded version
load('data/raw/filtBCmat113.mat');
data_pc(1).name = brainname{1};
% Injection site is col 16-28(APC),29-41(PPC)
src_id = 16:41;
prj_id = [1:6, 8:15];
% Load PC data
data_pc(1).srcImg = Bnorm113(:, src_id);
data_pc(1).srcRegName = {'APC', 'PPC'};
data_pc(1).nSrcRegSli = [13, 13];
% Load projection data
data_pc(1).prjImg = Bnorm113(:, prj_id);
data_pc(1).nPrjRegSli = ones(1, 14);
data_pc(1).prjRegName = regname;
% Get sequence
data_pc(1).data.brcSeq = cellfun(@(x) char(x), ...
    mat2cell(Bseq113, ones(size(Bseq113, 1), 1), size(Bseq113, 2)), ...
    'UniformOutput', false);
% Remove erroneus barcodes
unfit = any(sum(B113(:, src_id), 2) < B113(:, prj_id), 2);
fprintf('%i barcodes removed from %s\n', sum(unfit), data_pc(1).name);
data_pc(1).srcImg(unfit, :) = [];
data_pc(1).prjImg(unfit, :) = [];
data_pc(1).data.brcSeq(unfit, :) = [];
%---Non-thresholded version
load('data/raw/filtBCmat113noprojthresh-1.mat');
data_pcuf(1).name = brainname{1};
% Load PC data
data_pcuf(1).srcImg = Bnorm113(:, src_id);
data_pcuf(1).srcRegName = {'APC', 'PPC'};
data_pcuf(1).nSrcRegSli = [13, 13];
% Load projection data
data_pcuf(1).prjImg = Bnorm113(:, prj_id);
data_pcuf(1).nPrjRegSli = ones(1, 14);
data_pcuf(1).prjRegName = regname;
% Get sequence
data_pcuf(1).data.brcSeq = cellfun(@(x) char(x), ...
    mat2cell(Bseq113, ones(size(Bseq113, 1), 1), size(Bseq113, 2)), ...
    'UniformOutput', false);
% Remove erroneus barcodes; in this case; more than 3 slices in PC
unfit = sum(B113(:, src_id) == 0, 2) <= RM_THR;
fprintf('%i barcodes removed from %s (unfiltered)\n', sum(unfit), data_pcuf(1).name);
data_pcuf(1).srcImg(unfit, :) = [];
data_pcuf(1).prjImg(unfit, :) = [];
data_pcuf(1).data.brcSeq(unfit, :) = [];

% XC119
%---Thresholded version
load('data/raw/filtBCmat119.mat');
data_pc(2).name = brainname{2};
% Injection site is col 16-28(APC),29-41(PPC)
src_id = 16:41;
prj_id = [1:6, 8:15];
% Load PC data
data_pc(2).srcImg = Bnorm119(:, src_id);
data_pc(2).nSrcRegSli = [13, 13];
data_pc(2).srcRegName = {'APC', 'PPC'};
% Load projection data
data_pc(2).prjImg = Bnorm119(:, prj_id);
data_pc(2).nPrjRegSli = ones(1, 14);
data_pc(2).prjRegName = regname;
% Get sequence
data_pc(2).data.brcSeq = cellfun(@(x) char(x), ...
    mat2cell(Bseq119, ones(size(Bseq119, 1), 1), size(Bseq119, 2)), ...
    'UniformOutput', false);
% Remove erroneus barcodes
unfit = any(sum(B119(:, src_id), 2) < B119(:, prj_id), 2);
fprintf('%i barcodes removed from %s\n', sum(unfit), data_pc(2).name);
data_pc(2).srcImg(unfit, :) = [];
data_pc(2).prjImg(unfit, :) = [];
data_pc(2).data.brcSeq(unfit, :) = [];
%---Non-thresholded version
load('data/raw/filtBCmat119noprojthresh-1.mat');
data_pcuf(2).name = brainname{2};
% Load PC data
data_pcuf(2).srcImg = Bnorm119(:, src_id);
data_pcuf(2).nSrcRegSli = [13, 13];
data_pcuf(2).srcRegName = {'APC', 'PPC'};
% Load projection data
data_pcuf(2).prjImg = Bnorm119(:, prj_id);
data_pcuf(2).nPrjRegSli = ones(1, 14);
data_pcuf(2).prjRegName = regname;
% Get sequence
data_pcuf(2).data.brcSeq = cellfun(@(x) char(x), ...
    mat2cell(Bseq119, ones(size(Bseq119, 1), 1), size(Bseq119, 2)), ...
    'UniformOutput', false);
% Remove erroneus barcodes; in this case; more than 3 slices in PC
unfit = sum(B119(:, src_id) == 0, 2) <= RM_THR;
fprintf('%i barcodes removed from %s (unfiltered)\n', sum(unfit), data_pcuf(2).name);
data_pcuf(2).srcImg(unfit, :) = [];
data_pcuf(2).prjImg(unfit, :) = [];
data_pcuf(2).data.brcSeq(unfit, :) = [];

% XC120
%---Thresholded version
load('data/raw/filtBCmat120.mat');
data_pc(3).name = brainname{3};
% Injection site is col 16-28(APC),29-41(PPC)
src_id = 16:41;
prj_id = [1:6, 8:15];
% Load PC data
data_pc(3).srcImg = Bnorm120(:, src_id);
data_pc(3).nSrcRegSli = [13, 13];
data_pc(3).srcRegName = {'APC', 'PPC'};
% Load projection data
data_pc(3).prjImg = Bnorm120(:, prj_id);
data_pc(3).nPrjRegSli = ones(1, 14);
data_pc(3).prjRegName = regname;
% Get sequence
data_pc(3).data.brcSeq = cellfun(@(x) char(x), ...
    mat2cell(Bseq120, ones(size(Bseq120, 1), 1), size(Bseq120, 2)), ...
    'UniformOutput', false);
% Remove erroneus barcodes
unfit = any(sum(B120(:, src_id), 2) < B120(:, prj_id), 2);
fprintf('%i barcodes removed from %s\n', sum(unfit), data_pc(3).name);
data_pc(3).srcImg(unfit, :) = [];
data_pc(3).prjImg(unfit, :) = [];
data_pc(3).data.brcSeq(unfit, :) = [];
%---Non-Thresholded version
load('data/raw/filtBCmat120noprojthresh-1.mat');
data_pcuf(3).name = brainname{3};
% Load PC data
data_pcuf(3).srcImg = Bnorm120(:, src_id);
data_pcuf(3).nSrcRegSli = [13, 13];
data_pcuf(3).srcRegName = {'APC', 'PPC'};
% Load projection data
data_pcuf(3).prjImg = Bnorm120(:, prj_id);
data_pcuf(3).nPrjRegSli = ones(1, 14);
data_pcuf(3).prjRegName = regname;
% Get sequence
data_pcuf(3).data.brcSeq = cellfun(@(x) char(x), ...
    mat2cell(Bseq120, ones(size(Bseq120, 1), 1), size(Bseq120, 2)), ...
    'UniformOutput', false);
% Remove erroneus barcodes; in this case; more than 3 slices in PC
unfit = sum(B120(:, src_id) == 0, 2) <= RM_THR;
fprintf('%i barcodes removed from %s (unfiltered)\n', sum(unfit), data_pcuf(3).name);
data_pcuf(3).srcImg(unfit, :) = [];
data_pcuf(3).prjImg(unfit, :) = [];
data_pcuf(3).data.brcSeq(unfit, :) = [];

% YC123
%---Thresholded version
load('data/raw/filtBCmat123.mat');
data_pc(4).name = brainname{4};
% Injection site is col 17-29(APC),30-42(PPC)
src_id = 17:42;
prj_id = 1:14;
% Load PC data
data_pc(4).srcImg = Bnorm123(:, src_id);
data_pc(4).nSrcRegSli = [13, 13];
data_pc(4).srcRegName = {'APC', 'PPC'};
% Load projection data
data_pc(4).prjImg = Bnorm123(:, prj_id);
data_pc(4).nPrjRegSli = ones(1, 14);
data_pc(4).prjRegName = regname;
% Get sequence
data_pc(4).data.brcSeq = cellfun(@(x) char(x), ...
    mat2cell(Bseq123, ones(size(Bseq123, 1), 1), size(Bseq123, 2)), ...
    'UniformOutput', false);
% Remove erroneus barcodes
unfit = any(sum(B123(:, src_id), 2) < B123(:, prj_id), 2);
fprintf('%i barcodes removed from %s\n', sum(unfit), data_pc(4).name);
data_pc(4).srcImg(unfit, :) = [];
data_pc(4).prjImg(unfit, :) = [];
data_pc(4).data.brcSeq(unfit, :) = [];
%---Non-Thresholded version
load('data/raw/filtBCmat123noprojthresh-1.mat');
data_pcuf(4).name = brainname{4};
% Load PC data
data_pcuf(4).srcImg = Bnorm123(:, src_id);
data_pcuf(4).nSrcRegSli = [13, 13];
data_pcuf(4).srcRegName = {'APC', 'PPC'};
% Load projection data
data_pcuf(4).prjImg = Bnorm123(:, prj_id);
data_pcuf(4).nPrjRegSli = ones(1, 14);
data_pcuf(4).prjRegName = regname;
% Get sequence
data_pcuf(4).data.brcSeq = cellfun(@(x) char(x), ...
    mat2cell(Bseq123, ones(size(Bseq123, 1), 1), size(Bseq123, 2)), ...
    'UniformOutput', false);
% Remove erroneus barcodes; in this case; more than 3 slices in PC
unfit = sum(B123(:, src_id) == 0, 2) <= RM_THR;
fprintf('%i barcodes removed from %s (unfiltered)\n', sum(unfit), data_pcuf(4).name);
data_pcuf(4).srcImg(unfit, :) = [];
data_pcuf(4).prjImg(unfit, :) = [];
data_pcuf(4).data.brcSeq(unfit, :) = [];

% XC127
%---Thresholded version
load('data/raw/filtBCmat127.mat');
data_pc(5).name = brainname{5};
% Injection site is col 1-6,13,7-12(APC),14-26(PPC) [+27]
src_id = [1:6, 13, 7:12, 14:26];
% src_id = 1:26;
prj_id = 29:42;
% Load PC data
data_pc(5).srcImg = Bnorm(:, src_id);
data_pc(5).nSrcRegSli = [13, 13];
data_pc(5).srcRegName = {'APC', 'PPC'};
% Load projection data
data_pc(5).prjImg = Bnorm(:, prj_id);
data_pc(5).nPrjRegSli = ones(1, 14);
data_pc(5).prjRegName = regname;
% Get sequence
data_pc(5).data.brcSeq = cellfun(@(x) char(x), ...
    mat2cell(Bseq, ones(size(Bseq, 1), 1), size(Bseq, 2)), ...
    'UniformOutput', false);
% Remove erroneus barcodes
unfit = any(sum(B(:, src_id), 2) < B(:, prj_id), 2);
fprintf('%i barcodes removed from %s\n', sum(unfit), data_pc(5).name);
data_pc(5).srcImg(unfit, :) = [];
data_pc(5).prjImg(unfit, :) = [];
data_pc(5).data.brcSeq(unfit, :) = [];
%---Non-Thresholded version
load('data/raw/filtBCmat127noprojthresh-1.mat');
data_pcuf(5).name = brainname{5};
% Load PC data
data_pcuf(5).srcImg = B127norm(:, src_id);
data_pcuf(5).nSrcRegSli = [13, 13];
data_pcuf(5).srcRegName = {'APC', 'PPC'};
% Load projection data
data_pcuf(5).prjImg = B127norm(:, prj_id);
data_pcuf(5).nPrjRegSli = ones(1, 14);
data_pcuf(5).prjRegName = regname;
% Get sequence
data_pcuf(5).data.brcSeq = cellfun(@(x) char(x), ...
    mat2cell(B127seq, ones(size(B127seq, 1), 1), size(B127seq, 2)), ...
    'UniformOutput', false);
% Remove erroneus barcodes; in this case; more than 3 slices in PC
unfit = sum(B127(:, src_id) == 0, 2) <= RM_THR;
fprintf('%i barcodes removed from %s (unfiltered)\n', sum(unfit), data_pcuf(5).name);
data_pcuf(5).srcImg(unfit, :) = [];
data_pcuf(5).prjImg(unfit, :) = [];
data_pcuf(5).data.brcSeq(unfit, :) = [];

% Provide barcode id for brain labels
for i = 1:5
    % Write brain id
    data_pc(i).brcName = brainname;
    data_pc(i).brcId = i * ones(data_pc(i).nBrc, 1);
    % Write brain id
    data_pcuf(i).brcName = brainname;
    data_pcuf(i).brcId = i * ones(data_pcuf(i).nBrc, 1);
end


%% Alignment and getting of OB data
% Data of interest
% % DAT = {'65', 'BC61', 'B86', 'B92', 'B109', 'B111'};
DAT = {'yc65', 'yc61', 'yc86', 'yc92', 'yc109', 'yc111'};

% Some parameters
FACTOR = 2;         % Blow up to this factor

% thisone the data that's wanted
thisone = data_ob(arrayfun(@(x) any(strcmp(x.name, DAT)), data_ob));

% Merge dorso-ventral thisones
% dvs = arrayfun(@(x) any(strcmp('OB dorsal', x.srcRegName)), thisone);
% barNames = cell(length(thisone) + sum(dvs), 1);
% D = 1;
% for s = 1:length(thisone)
%     if dvs(s)
%         % Find dorsoventral ID's
%         dor_id = find(strcmp(thisone(s).srcRegName, 'OB dorsal'));
%         ven_id = find(strcmp(thisone(s).srcRegName, 'OB ventral'));
%         % Add dorso-ventral images on top of each other
%         imsum = ...
%             thisone(s).getImgSrcReg('OB dorsal') + ...
%             thisone(s).getImgSrcReg('OB ventral');
%         % Replace dorsal thisone info
%         thisone(s).srcImg(:, thisone(s).srcRegInd{dor_id}) = imsum;
%         thisone(s).srcRegName{dor_id} = 'OB';
%         % Remove ventral thisones
%         thisone(s).srcImg(:, thisone(s).srcRegInd{ven_id}) = [];
%         thisone(s).nSrcRegSli(ven_id) = [];
%         thisone(s).srcRegName(ven_id) = [];
%         % Find the dorsoventral markers
%         ids = unique(thisone(s).brcId);
%         nid = zeros(thisone(s).nBrc, 1);
%         nid(thisone(s).brcId == ids(1)) = D;
%         nid(thisone(s).brcId == ids(2)) = D + 1;
%         thisone(s).brcId = nid;
%         barNames{D} = [thisone(s).name, ' dorsal'];
%         barNames{D + 1} = [thisone(s).name, ' ventral'];
%         D = D + 2;
%     else
%         % Relabel markers
%         thisone(s).brcId = (0 .* thisone(s).brcId) + D;
%         barNames{D} = thisone(s).name;
%         D = D + 1;
%     end
% end

% % Find common region to all
% prjRegName = unique(vertcat(thisone.prjRegName), 'stable');
% prjRegName = prjRegName(cellfun( ...
%     @(x) all(arrayfun(@(y) any(strcmp(x, y.prjRegName)), thisone)), ...
%     prjRegName));
% srcRegName = unique(vertcat(thisone.srcRegName), 'stable');
% srcRegName = srcRegName(cellfun( ...
%     @(x) all(arrayfun(@(y) any(strcmp(x, y.srcRegName)), thisone)), ...
%     srcRegName));
% % Trim extrenous regions
% for s = 1:length(thisone)
%     prjInd = cellfun( ...
%         @(x) find(strcmp(x, thisone(s).prjRegName), 1), prjRegName);
%     thisone(s).prjImg = thisone(s).getImgPrjReg(prjRegName);
%     thisone(s).nPrjRegSli = thisone(s).nPrjRegSli(prjInd);
%     thisone(s).prjRegName = prjRegName;
%     srcInd = cellfun( ...
%         @(x) find(strcmp(x, thisone(s).srcRegName), 1), srcRegName);
%     thisone(s).srcImg = thisone(s).getImgSrcReg(srcRegName);
%     thisone(s).nSrcRegSli = thisone(s).nSrcRegSli(srcInd);
%     thisone(s).srcRegName = srcRegName;
% end
% 
% % Find new region sizes
% newNSrcRegSli = FACTOR * max(horzcat(thisone.nSrcRegSli), [], 2);
% newNPrjRegSli = FACTOR * max(horzcat(thisone.nPrjRegSli), [], 2);
% 
% % Align the data images to the new sizes
% for s = 1:length(thisone)
%     % Align the source region; do it region by region
%     srcImg = zeros(thisone(s).nBrc, sum(newNSrcRegSli));
%     for r = 1:length(newNSrcRegSli)
%         srcImg(:, (1:newNSrcRegSli(r)) + sum(newNSrcRegSli(1:(r-1)))) = ...
%             interp1(...
%             linspace(0, 1, thisone(s).nSrcRegSli(r)), ...
%             thisone(s).getImgSrcReg(srcRegName{r})', ...
%             linspace(0, 1, newNSrcRegSli(r)), ...
%             'pchip')';
%     end
%     thisone(s).srcImg = srcImg;
%     thisone(s).nSrcRegSli = newNSrcRegSli;
%     % Align the projection region; do it region by region
%     prjImg = zeros(thisone(s).nBrc, sum(newNPrjRegSli));
%     for r = 1:length(newNPrjRegSli)
%         prjImg(:, (1:newNPrjRegSli(r)) + sum(newNPrjRegSli(1:(r-1)))) = ...
%             interp1(...
%             linspace(0, 1, thisone(s).nPrjRegSli(r)), ...
%             thisone(s).getImgPrjReg(prjRegName{r})', ...
%             linspace(0, 1, newNPrjRegSli(r)), ...
%             'pchip')';
%     end
%     thisone(s).prjImg = prjImg;
%     thisone(s).nPrjRegSli = newNPrjRegSli;
% end

% Create a merged instance of all the data
% new_thisone = thisone(1);
% for s = 2:length(thisone)
%     new_thisone = new_thisone + thisone(s);
% end
% thisone = [thisone; new_thisone];
% thisone(end).name = 'Merged';
sto_barInfo = [thisone.data];
thisone(end+1) = mapseqData;
thisone(end).name = 'Merged';
thisone(end).srcImg = vertcat(thisone(1:(end-1)).srcImg);
thisone(end).srcRegName = thisone(1).srcRegName;
thisone(end).nSrcRegSli = thisone(1).nSrcRegSli;
thisone(end).prjImg = vertcat(thisone(1:(end-1)).prjImg);
thisone(end).prjRegName = thisone(1).prjRegName;
thisone(end).nPrjRegSli = thisone(1).nPrjRegSli;
thisone(end).brcId = vertcat(thisone(1:(end-1)).brcId);
thisone(end).brcCor = vertcat(thisone(1:(end-1)).brcCor);
thisone(end).brcName = thisone(1).brcName;
thisone(end).prjRegSum = vertcat(thisone(1:(end-1)).prjRegSum);
thisone(end).data.brcDV = vertcat(sto_barInfo.brcDV);
thisone(end).data.brcSoma = vertcat(sto_barInfo.brcSoma);
% Reorder according to the brightest thisone in projection
for s = 1:length(thisone)
    % Get the reordering
    [~, sto] = max(thisone(s).prjImg, [], 2, 'omitnan');
    [~, bor] = sort(sto, 'ascend');
    % Rewrite the images
    thisone(s).srcImg = thisone(s).srcImg(bor, :);
    thisone(s).prjImg = thisone(s).prjImg(bor, :);
    thisone(s).brcId  = thisone(s).brcId( bor, :);
    thisone(s).prjRegSum = thisone(s).prjRegSum(bor, :);
    thisone(s).brcCor = thisone(s).brcCor(bor, :);
    thisone(s).brcName = thisone(1).brcName;
    thisone(s).data.brcDV = thisone(s).data.brcDV(bor, :);
    thisone(s).data.brcSoma = thisone(s).data.brcSoma(bor, :);
end

% Create the aligned object, and find barcode labels
aliOB = mapseqDataAligned(thisone);
aliOB.brcName = thisone(1).brcName;
aliOB.data.template = template;
sto_barInfo = [thisone.data];
aliOB.data.brcDV = {sto_barInfo.brcDV}';
aliOB.data.brcSoma = {sto_barInfo.brcSoma}';

%% Alignment of the PC data
% Data of interest
DAT = {'XC113', 'XC119', 'XC120', 'YC123', 'XC127'};

% The data should be aligned already; just take it in
thisone = data_pc(  arrayfun(@(x) any(strcmp(x.name, DAT)), data_pc));
altered = data_pcuf(arrayfun(@(x) any(strcmp(x.name, DAT)), data_pcuf));

% Fix barcode labels
D = 1;
for s = 1:length(thisone)
    thisone(s).brcId = (0 .* thisone(s).brcId) + D;
    altered(s).brcId = (0 .* altered(s).brcId) + D;
    D = D + 1;
end

% Create merged dataset
new_thisone = thisone(1) + thisone(2);
for s = 3:length(thisone)
new_thisone = new_thisone + thisone(s);
end
thisone = [thisone; new_thisone];
thisone(end).name = 'Merged';
new_altered = altered(1) + altered(2);
for s = 3:length(altered)
new_altered = new_altered + altered(s);
end
altered = [altered; new_altered];
altered(end).name = 'Merged';

% Reorder according to the brightest thisone in source
for s = 1:length(thisone)
    % Get the reordering
    [~, sto] = max(thisone(s).srcImg, [], 2, 'omitnan');
    [~, bor] = sort(sto, 'ascend');
    % Rewrite the barcodes
    thisone(s).reorderBarcodes(bor);
    % Get the reordering
    [~, sto] = max(altered(s).srcImg, [], 2, 'omitnan');
    [~, bor] = sort(sto, 'ascend');
    % Rewrite the barcodes
    altered(s).reorderBarcodes(bor);
end

% Create the aligned object, and find barcode labels
aliPC = mapseqDataAligned(thisone);
% Barcode name identifier
aliPC.brcName = thisone(1).brcName;
aPCAll = mapseqDataAligned(altered);
aPCAll.brcName = altered(1).brcName;

% Reformat the template for aliOB
%--Classifier neural network
sto_templates = aliOB.data.template;

% Reformat template class
aliOB.data.temp.names = {sto_templates.name}';
aliOB.data.temp.prjReg = vertcat(sto_templates.data);
aliOB.data.temp.location = vertcat(sto_templates.pos);
aliOB.data.temp.nIds = length(aliOB.data.temp.names);
aliOB.data.temp.nBrc = size(aliOB.data.temp.prjReg, 1);
aliOB.data.temp.realId = cell2mat(arrayfun( ...
  @(x, y) y * ones(size(x.data, 1), 1), ...
  sto_templates', (1:3)', 'UniformOutput', 0));
aliOB.data.temp.tsneCor = vertcat(sto_templates.cor);
% Conversion (from Baryo := 3D) to lower dimensional (2D) model
[ aliOB.data.temp.conv_fromBaryo, ...
  aliOB.data.temp.conv_toBaryo] = ...
  aux.genSimplex(aliOB.data.temp.nIds - 1);

%% Save data
data = vertcat(data_ob, data_pc);
barseq = data_bs;

% Concatanate and save the formatted data
save('data/current.mat', 'data', 'data_ob', 'data_pc', 'data_pcuf', ...
    'barseq', 'aliOB', 'aliPC', 'aPCAll');


