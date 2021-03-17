function newstatNN(SRC, MAIN, CLASS, BOOT)
%NEWSTATNN Get data; accross neural network types on how the analysis
%changes
% For each (t)ype of neural network; per one (n)etwork;
% 	* Check the TPR (TP/P) and FPR (FP/N) on the annotated dataset
%   * Check the network loss as a classifier
%   * Number of classified mitral cells
%   * Get the ob fits; and the p value wrt the main result using bootstrap

%#ok<*UNRCH>
rng('default');
TEMP = MAIN.data.temp;
NAME = {'16', '20', '10-10'};

% Main dataset information
regmat_main = [ ...
  sum(MAIN.prjImg(:, MAIN.prjRegInd{1}), 2), ...
  sum(MAIN.prjImg(:, MAIN.prjRegInd{4}), 2), ...
  sum(MAIN.prjImg(:, MAIN.prjRegInd{5}), 2), ...
  sum(MAIN.prjImg(:, MAIN.prjRegInd{6}), 2)];
apcmat_main = MAIN.prjImg(:, MAIN.prjRegInd{2});
ppcmat_main = MAIN.prjImg(:, MAIN.prjRegInd{3});
pcmat_main = [apcmat_main, ppcmat_main];

% Size information
[n_net, n_type] = size(CLASS);
n_reg = size(regmat_main, 2);
n_apc = size(apcmat_main, 2);
n_ppc = size(ppcmat_main, 2);
n_sli = size(pcmat_main, 2);
% Limits for fitting
l_spl = [1 - 1e-5; n_apc + .5; n_sli];
l_lin = [1 - 1e-5; n_sli];
l_apc = [1 - 1e-5; n_apc];
l_ppc = [1 - 1e-5; n_ppc];

% Store the output of the networks here
outstr = struct();
for t = n_type:-1:1
  for n = n_net:-1:1
    % Get the dataset for this neural network
    SET = aux.obGetType(SRC, CLASS(n, t));
    outstr(n, t).name = NAME{t};
    % Do TPR and FPR on the template data
    outstr(n, t).fpr = 1 - CLASS(n, t).mitral_tnr;
    outstr(n, t).tnr = CLASS(n, t).mitral_tnr;
    outstr(n, t).tpr = CLASS(n, t).mitral_tpr;
    outstr(n, t).loss = CLASS(n, t).mitral_tpr;
    % Dataset information
    regmat = [ ...
      sum(SET.prjImg(:, SET.prjRegInd{1}), 2), ...
      sum(SET.prjImg(:, SET.prjRegInd{4}), 2), ...
      sum(SET.prjImg(:, SET.prjRegInd{5}), 2), ...
      sum(SET.prjImg(:, SET.prjRegInd{6}), 2)];
    apcmat = SET.prjImg(:, SET.prjRegInd{2});
    ppcmat = SET.prjImg(:, SET.prjRegInd{3});
    pcmat = [apcmat, ppcmat];
    % Calculate the full OB matrix
    n_brc = size(regmat, 1);
    outstr(n, t).nBrc = n_brc;
    outstr(n, t).conProb_pc = aux.conProb(regmat, pcmat);
    outstr(n, t).conProb_apc = outstr(n, t).conProb_pc(:, 1:n_apc);
    outstr(n, t).conProb_ppc = outstr(n, t).conProb_pc(:, n_apc + (1:n_ppc));
    % Fit several lines using bootstrap
    s_bootApc = ceil(n_apc * rand(BOOT, n_apc));
    s_bootPpc = ceil(n_ppc * rand(BOOT, n_ppc));
    s_bootPc = [s_bootApc, n_apc + s_bootPpc];
    cpFit_b_spl = zeros(3, n_reg, BOOT);
    cpFit_b_lin = zeros(2, n_reg, BOOT);
    cpFit_b_apc = zeros(2, n_reg, BOOT);
    cpFit_b_ppc = zeros(2, n_reg, BOOT);
    cpFit_m_spl = zeros(2, n_reg, BOOT);
    cpFit_m_lin = zeros(1, n_reg, BOOT);
    cpFit_m_apc = zeros(1, n_reg, BOOT);
    cpFit_m_ppc = zeros(1, n_reg, BOOT);
    for r = 1:n_reg
      for i = 1:BOOT
        [cpFit_m_spl(:, r, i), cpFit_b_spl(:, r, i)] = aux.splineFit( ...
          s_bootPc(i, :), outstr(n, t).conProb_pc(r, s_bootPc(i, :)), l_spl);
        [cpFit_m_lin(:, r, i), cpFit_b_lin(:, r, i)] = aux.splineFit( ...
          s_bootPc(i, :), outstr(n, t).conProb_pc(r, s_bootPc(i, :)), l_lin);
        [cpFit_m_apc(:, r, i), cpFit_b_apc(:, r, i)] = aux.splineFit( ...
          s_bootApc(i, :), outstr(n, t).conProb_apc(r, s_bootApc(i, :)), l_apc);
        [cpFit_m_ppc(:, r, i), cpFit_b_ppc(:, r, i)] = aux.splineFit( ...
          s_bootPpc(i, :), outstr(n, t).conProb_ppc(r, s_bootPpc(i, :)), l_ppc);
      end
    end
    outstr(n, t).cpFit_b_spline = mean(cpFit_b_spl, 3);
    outstr(n, t).cpFit_m_spline = mean(cpFit_m_spl, 3);
    outstr(n, t).cpFit_b_linear = mean(cpFit_b_lin, 3);
    outstr(n, t).cpFit_m_linear = mean(cpFit_m_lin, 3);
    outstr(n, t).cpFit_b_apc = mean(cpFit_b_apc, 3);
    outstr(n, t).cpFit_m_apc = mean(cpFit_m_apc, 3);
    outstr(n, t).cpFit_b_ppc = mean(cpFit_b_ppc, 3);
    outstr(n, t).cpFit_m_ppc = mean(cpFit_m_ppc, 3);
    
    outstr(n, t).cpFit_b_spline_std = std(cpFit_b_spl, 0, 3);
    outstr(n, t).cpFit_m_spline_std = std(cpFit_m_spl, 0, 3);
    outstr(n, t).cpFit_b_linear_std = std(cpFit_b_lin, 0, 3);
    outstr(n, t).cpFit_m_linear_std = std(cpFit_m_lin, 0, 3);
    outstr(n, t).cpFit_b_apc_std = std(cpFit_b_apc, 0, 3);
    outstr(n, t).cpFit_m_apc_std = std(cpFit_m_apc, 0, 3);
    outstr(n, t).cpFit_b_ppc_std = std(cpFit_b_ppc, 0, 3);
    outstr(n, t).cpFit_m_ppc_std = std(cpFit_m_ppc, 0, 3);
    
    outstr(n, t).cpFit_b_spline_pval = zeros(3, n_reg);
    outstr(n, t).cpFit_m_spline_pval = zeros(2, n_reg);
    outstr(n, t).cpFit_b_linear_pval = zeros(2, n_reg);
    outstr(n, t).cpFit_m_linear_pval = zeros(1, n_reg);
    outstr(n, t).cpFit_b_apc_pval = zeros(2, n_reg);
    outstr(n, t).cpFit_m_apc_pval = zeros(1, n_reg);
    outstr(n, t).cpFit_b_ppc_pval = zeros(2, n_reg);
    outstr(n, t).cpFit_m_ppc_pval = zeros(1, n_reg);
    
    thisCalc = @(dat, nh) 1 - erf(abs(mean(( ...
      dat - nh), 3) ./ ...
      std(dat, 0, 3) ./ sqrt(2)));
    outstr(n, t).cpFit_b_spline_pval(:, :) = thisCalc(cpFit_b_spl, ...
      MAIN.data.OBPC.cpFit_b_spline);
    outstr(n, t).cpFit_m_spline_pval(:, :) = thisCalc(cpFit_m_spl, ...
      MAIN.data.OBPC.cpFit_m_spline);
    outstr(n, t).cpFit_b_linear_pval(:, :) = thisCalc(cpFit_b_lin, ...
      MAIN.data.OBPC.cpFit_b_linear);
    outstr(n, t).cpFit_m_linear_pval(:, :) = thisCalc(cpFit_m_lin, ...
      MAIN.data.OBPC.cpFit_m_linear);
    outstr(n, t).cpFit_b_apc_pval(:, :) = thisCalc(cpFit_b_apc, ...
      MAIN.data.OBPC.cpFit_b_apc);
    outstr(n, t).cpFit_m_apc_pval(:, :) = thisCalc(cpFit_m_apc, ...
      MAIN.data.OBPC.cpFit_m_apc);
    outstr(n, t).cpFit_b_ppc_pval(:, :) = thisCalc(cpFit_b_ppc, ...
      MAIN.data.OBPC.cpFit_b_ppc);
    outstr(n, t).cpFit_m_ppc_pval(:, :) = thisCalc(cpFit_m_ppc, ...
      MAIN.data.OBPC.cpFit_m_ppc);
  end
end

% Analyse based on network type
% 	* Check the TPR (TP/P) and FPR (FP/N) on the annotated dataset
%   * Check the network loss as a classifier
%   * Number of classified mitral cells
%   * Get the ob fits; and the p value wrt the main result using bootstrap
datstr = struct();
for t = 1:n_type
  datstr(t).name = outstr(1, t).name;
  % Get tpr, fpr, loss and barcode number
  sto_tpr = [outstr(:, t).tpr];
  datstr(t).tpr     = mean(sto_tpr);
  datstr(t).tpr_std = std(sto_tpr);
  sto_fpr = [outstr(:, t).fpr];
  datstr(t).fpr     = mean(sto_fpr);
  datstr(t).fpr_std = std(sto_fpr);
  sto_brc = [outstr(:, t).nBrc];
  datstr(t).nBrc     = mean(sto_brc);
  datstr(t).nBrc_std = std(sto_brc);
  sto_loss = [outstr(:, t).loss];
  datstr(t).loss     = mean(sto_loss);
  datstr(t).loss_std = std(sto_loss);
  % Get fit info
  conProb = zeros(n_reg, n_sli, n_net);
  cpFit_b_spl = zeros(3, n_reg, n_net);
  cpFit_b_lin = zeros(2, n_reg, n_net);
  cpFit_b_apc = zeros(2, n_reg, n_net);
  cpFit_b_ppc = zeros(2, n_reg, n_net);
  cpFit_m_spl = zeros(2, n_reg, n_net);
  cpFit_m_lin = zeros(1, n_reg, n_net);
  cpFit_m_apc = zeros(1, n_reg, n_net);
  cpFit_m_ppc = zeros(1, n_reg, n_net);
  cpFit_b_spl_pval = zeros(3, n_reg, n_net);
  cpFit_b_lin_pval = zeros(2, n_reg, n_net);
  cpFit_b_apc_pval = zeros(2, n_reg, n_net);
  cpFit_b_ppc_pval = zeros(2, n_reg, n_net);
  cpFit_m_spl_pval = zeros(2, n_reg, n_net);
  cpFit_m_lin_pval = zeros(1, n_reg, n_net);
  cpFit_m_apc_pval = zeros(1, n_reg, n_net);
  cpFit_m_ppc_pval = zeros(1, n_reg, n_net);
  for n = 1:n_net
    conProb(:, :, n) = outstr(n, t).conProb_pc;
    cpFit_b_spl(:, :, n) = outstr(n, t).cpFit_b_spline;
    cpFit_b_lin(:, :, n) = outstr(n, t).cpFit_b_linear;
    cpFit_b_apc(:, :, n) = outstr(n, t).cpFit_b_apc;
    cpFit_b_ppc(:, :, n) = outstr(n, t).cpFit_b_ppc;
    cpFit_m_spl(:, :, n) = outstr(n, t).cpFit_m_spline;
    cpFit_m_lin(:, :, n) = outstr(n, t).cpFit_m_linear;
    cpFit_m_apc(:, :, n) = outstr(n, t).cpFit_m_apc;
    cpFit_m_ppc(:, :, n) = outstr(n, t).cpFit_m_ppc;
    cpFit_b_spl_pval(:, :, n) = outstr(n, t).cpFit_b_spline_pval;
    cpFit_b_lin_pval(:, :, n) = outstr(n, t).cpFit_b_linear_pval;
    cpFit_b_apc_pval(:, :, n) = outstr(n, t).cpFit_b_apc_pval;
    cpFit_b_ppc_pval(:, :, n) = outstr(n, t).cpFit_b_ppc_pval;
    cpFit_m_spl_pval(:, :, n) = outstr(n, t).cpFit_m_spline_pval;
    cpFit_m_lin_pval(:, :, n) = outstr(n, t).cpFit_m_linear_pval;
    cpFit_m_apc_pval(:, :, n) = outstr(n, t).cpFit_m_apc_pval;
    cpFit_m_ppc_pval(:, :, n) = outstr(n, t).cpFit_m_ppc_pval;
  end
  datstr(t).conProb = mean(conProb, 3);
  datstr(t).conProb_std = std(conProb, 0, 3);
  % Intercepts
  datstr(t).cpFit_b_spline = mean(cpFit_b_spl, 3);
  datstr(t).cpFit_b_linear = mean(cpFit_b_lin, 3);
  datstr(t).cpFit_b_apc = mean(cpFit_b_apc, 3);
  datstr(t).cpFit_b_ppc = mean(cpFit_b_ppc, 3);
  datstr(t).cpFit_b_spline_std = std(cpFit_b_spl, 0, 3);
  datstr(t).cpFit_b_linear_std = std(cpFit_b_lin, 0, 3);
  datstr(t).cpFit_b_apc_std = std(cpFit_b_apc, 0, 3);
  datstr(t).cpFit_b_ppc_std = std(cpFit_b_ppc, 0, 3);
  datstr(t).cpFit_b_spline_pval = mean(cpFit_b_spl_pval, 3);
  datstr(t).cpFit_b_linear_pval = mean(cpFit_b_lin_pval, 3);
  datstr(t).cpFit_b_apc_pval = mean(cpFit_b_apc_pval, 3);
  datstr(t).cpFit_b_ppc_pval = mean(cpFit_b_ppc_pval, 3);
  datstr(t).cpFit_b_spline_pval_std = std(cpFit_b_spl_pval, 0, 3);
  datstr(t).cpFit_b_linear_pval_std = std(cpFit_b_lin_pval, 0, 3);
  datstr(t).cpFit_b_apc_pval_std = std(cpFit_b_apc_pval, 0, 3);
  datstr(t).cpFit_b_ppc_pval_std = std(cpFit_b_ppc_pval, 0, 3);
  % Slopes
  datstr(t).cpFit_m_spline = mean(cpFit_m_spl, 3);
  datstr(t).cpFit_m_linear = mean(cpFit_m_lin, 3);
  datstr(t).cpFit_m_apc = mean(cpFit_m_apc, 3);
  datstr(t).cpFit_m_ppc = mean(cpFit_m_ppc, 3);
  datstr(t).cpFit_m_spline_std = std(cpFit_m_spl, 0, 3);
  datstr(t).cpFit_m_linear_std = std(cpFit_m_lin, 0, 3);
  datstr(t).cpFit_m_apc_std = std(cpFit_m_apc, 0, 3);
  datstr(t).cpFit_m_ppc_std = std(cpFit_m_ppc, 0, 3);
  datstr(t).cpFit_m_spline_pval = mean(cpFit_m_spl_pval, 3);
  datstr(t).cpFit_m_linear_pval = mean(cpFit_m_lin_pval, 3);
  datstr(t).cpFit_m_apc_pval = mean(cpFit_m_apc_pval, 3);
  datstr(t).cpFit_m_ppc_pval = mean(cpFit_m_ppc_pval, 3);
  datstr(t).cpFit_m_spline_pval_std = std(cpFit_m_spl_pval, 0, 3);
  datstr(t).cpFit_m_linear_pval_std = std(cpFit_m_lin_pval, 0, 3);
  datstr(t).cpFit_m_apc_pval_std = std(cpFit_m_apc_pval, 0, 3);
  datstr(t).cpFit_m_ppc_pval_std = std(cpFit_m_ppc_pval, 0, 3);
end

MAIN.data.allnet = outstr;
MAIN.data.typenet = datstr;
rng('shuffle');

end

