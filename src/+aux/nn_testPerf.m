function nn_testPerf(DAT, NHY, NET, THR)
  %NN_TESTPERF Runs statictical testing on results from NN
  % INPUTS
  %   DAT: This will probably be obdata; mapseqData object in question
  %   NHY: This will probably be obdata_mitral; the mapseqData object with
  %     the full analysis already run
  %   NET: List that contains all the neural networks used
  %   THR: Threshold value used to seperate if a cell is mitral
  % Results are written in the DAT variable; which is a handle class
  
  [Ns, Nn] = size(NET);
  null_conProb = NHY.data.conProb;
  sto_m = find(strcmp(DAT.data.temp.names, 'Mitral'), 1);
  
  
  % Will calculate, for now
  %   the loss on the template set
  DAT.data.nn_loss = reshape([NET.net_loss], Ns, Nn);
  %   p value of the corresponding prob matrix (wrt NHY)
  DAT.data.nn_conProb = cell(Ns, Nn);
  DAT.data.nn_conProb_pVal = zeros(Ns, Nn);
  %   # of neurons identified per class
  DAT.data.nn_nBrcPerClass = zeros(Ns, Nn, 3);
  %   The tsne coordinates of the mitral barcodes
  DAT.data.nn_brcMitCor = cell(Ns, Nn);
  for s = 1:Ns
    for n = 1:Nn
      sto_net = NET(s, n).net;
      sto_prob = sto_net(DAT.prjRegSum')';
      % Check how many neurons per class for this classifier
      DAT.data.nn_nBrcPerClass(s, n, :) = sum(sto_prob >= THR, 1);
      % Create a seperate mitral dataset
      sto = sto_prob(:, sto_m) >= THR;
      sto_mit = mapseqData;
      sto_mit.srcImg = DAT.srcImg(sto, :);
      sto_mit.srcRegName = DAT.srcRegName;
      sto_mit.nSrcRegSli = DAT.nSrcRegSli;
      sto_mit.prjImg = DAT.prjImg(sto, :);
      sto_mit.prjRegName = DAT.prjRegName;
      sto_mit.nPrjRegSli = DAT.nPrjRegSli;
      sto_mit.brcCor = DAT.brcCor(sto, :);
      % Get the tsne coordinates
      DAT.data.nn_brcMitCor{s, n} = sto_mit.brcCor;
      % Get the probability matrix for PC for this one
      sto_regmat = [ ...
        sum(sto_mit.prjImg(:, sto_mit.prjRegInd{1}), 2), ...
        sum(sto_mit.prjImg(:, sto_mit.prjRegInd{4}), 2), ...
        sum(sto_mit.prjImg(:, sto_mit.prjRegInd{5}), 2), ...
        sum(sto_mit.prjImg(:, sto_mit.prjRegInd{6}), 2)];
      sto_pcmat = sto_mit.prjImg(:, [sto_mit.prjRegInd{2:3}]);
      DAT.data.nn_conProb{s, n} = aux.conProb(sto_regmat, sto_pcmat);
      % Record the p value on this run wrt ground truth
      [~, DAT.data.nn_conProb_pVal(s, n), ~] = ttest( ...
        DAT.data.nn_conProb{s, n}(:) - null_conProb(:));
    end
  end
  
end

