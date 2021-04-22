function newstatPCout(set, BOOT)
%NEWSTATOB Three things that are, ugh . . .
% 1. Spearman correlation
% 2. Bootstrap slices for fit slope; check how often the slope changes
%   signs
% 3. Shuffle slices to calculate pvalue wrt shuffling
%#ok<*UNRCH>
rng('default');
USE_NORMALPDF = true;
CROSSVAL = 10;
BIN_SRC = true;
SKIP = 2;
PMIN = 1/BOOT;
pval_est = @(t, n) 2 * min(tcdf(t.*sqrt((n-2)./(1-t.^2)), n), tcdf(t.*sqrt((n-2)./(1-t.^2)), n, 'upper'));

% Store variables here
outstr = struct();

% Get the PC image to work with
if BIN_SRC
  pcmat = zeros(size(set.srcImg));
  pcmat(sub2ind(size(set.srcImg), (1:set.nBrc)', set.brcSrcVis)) = 1;
else
  pcmat = set.srcImg;
end
apcmat = pcmat(:, set.srcRegInd{1});
ppcmat = pcmat(:, set.srcRegInd{2});

% Collapse the regions into one; use only mitral cell data
regions = {'AON', 'OT', 'CoA', 'ENT'};
regmat = set.getImgPrjReg(regions);

% Useful
n_brc = size(regmat, 1);
n_reg = size(regmat, 2);
n_apc = size(apcmat, 2);
n_ppc = size(ppcmat, 2);
n_sli = size(pcmat, 2);
% Limits for fitting
l_spl = [1 + SKIP - 1e-5; n_apc + .5; n_sli - SKIP];
l_lin = [1 + SKIP - 1e-5; n_sli - SKIP];
l_apc = [1 + SKIP - 1e-5; n_apc + .5];
l_ppc = [n_apc + .5; n_sli - SKIP];
% Indices to get the skipped ones
x_fpc = ((1 + SKIP):(n_sli - SKIP))';
x_apc = ((1 + SKIP):n_apc)';
x_ppc = ((1 + n_apc):(n_sli - SKIP))';

% Calculate the full PC output matrix
outstr.conProb_pc = aux.conProb(regmat, pcmat);
outstr.conProb_apc = outstr.conProb_pc(:, 1:n_apc);
outstr.conProb_ppc = outstr.conProb_pc(:, n_apc + (1:n_ppc));

% 1. Calculate spearman correlation between regions and PC
[outstr.conProb_pc_sprCorr, outstr.conProb_pc_sprCorrPval] = ...
  corr(outstr.conProb_pc(:, x_fpc)', x_fpc, 'Type', 'Spearman');
[outstr.conProb_apc_sprCorr, outstr.conProb_apc_sprCorrPval] = ...
  corr(outstr.conProb_pc(:, x_apc)', x_apc, 'Type', 'Spearman');
[outstr.conProb_ppc_sprCorr, outstr.conProb_ppc_sprCorrPval] = ...
  corr(outstr.conProb_pc(:, x_ppc)', x_ppc, 'Type', 'Spearman');
outstr.conProb_pc_sprCorrPval2 = pval_est(outstr.conProb_pc_sprCorr, length(x_fpc));
outstr.conProb_apc_sprCorrPval2 = pval_est(outstr.conProb_apc_sprCorr, length(x_apc));
outstr.conProb_ppc_sprCorrPval2 = pval_est(outstr.conProb_ppc_sprCorr, length(x_ppc));

% 2. Fit several lines using bootstrap
s_bootApc = x_apc(ceil((n_apc - SKIP) * rand(BOOT, (n_apc - SKIP))));
s_bootPpc = x_ppc(ceil((n_ppc - SKIP) * rand(BOOT, (n_ppc - SKIP))));
s_bootPc = [s_bootApc, s_bootPpc];
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
      s_bootPc(i, :), outstr.conProb_pc(r, s_bootPc(i, :)), l_spl);
    [cpFit_m_lin(:, r, i), cpFit_b_lin(:, r, i)] = aux.splineFit( ...
      s_bootPc(i, :), outstr.conProb_pc(r, s_bootPc(i, :)), l_lin);
    [cpFit_m_apc(:, r, i), cpFit_b_apc(:, r, i)] = aux.splineFit( ...
      s_bootApc(i, :), outstr.conProb_pc(r, s_bootApc(i, :)), l_apc);
    [cpFit_m_ppc(:, r, i), cpFit_b_ppc(:, r, i)] = aux.splineFit( ...
      s_bootPpc(i, :), outstr.conProb_pc(r, s_bootPpc(i, :)), l_ppc);
  end
end
outstr.cpFit_b_spline = mean(cpFit_b_spl, 3);
outstr.cpFit_m_spline = mean(cpFit_m_spl, 3);
outstr.cpFit_b_linear = mean(cpFit_b_lin, 3);
outstr.cpFit_m_linear = mean(cpFit_m_lin, 3);
outstr.cpFit_b_apc = mean(cpFit_b_apc, 3);
outstr.cpFit_m_apc = mean(cpFit_m_apc, 3);
outstr.cpFit_b_ppc = mean(cpFit_b_ppc, 3);
outstr.cpFit_m_ppc = mean(cpFit_m_ppc, 3);

outstr.cpFit_b_spline_std = std(cpFit_b_spl, 0, 3);
outstr.cpFit_m_spline_std = std(cpFit_m_spl, 0, 3);
outstr.cpFit_b_linear_std = std(cpFit_b_lin, 0, 3);
outstr.cpFit_m_linear_std = std(cpFit_m_lin, 0, 3);
outstr.cpFit_b_apc_std = std(cpFit_b_apc, 0, 3);
outstr.cpFit_m_apc_std = std(cpFit_m_apc, 0, 3);
outstr.cpFit_b_ppc_std = std(cpFit_b_ppc, 0, 3);
outstr.cpFit_m_ppc_std = std(cpFit_m_ppc, 0, 3);

% Calculate p values of these quantities
outstr.cpFit_b_spline_pval = zeros(3, n_reg, 2);
outstr.cpFit_m_spline_pval = zeros(2, n_reg, 2);
outstr.cpFit_b_linear_pval = zeros(2, n_reg, 2);
outstr.cpFit_m_linear_pval = zeros(1, n_reg, 2);
outstr.cpFit_b_apc_pval = zeros(2, n_reg, 2);
outstr.cpFit_m_apc_pval = zeros(1, n_reg, 2);
outstr.cpFit_b_ppc_pval = zeros(2, n_reg, 2);
outstr.cpFit_m_ppc_pval = zeros(1, n_reg, 2);

if USE_NORMALPDF
  thisCalc = @(dat) 1 - erf(abs(mean(dat, 3)) ./ std(dat, 0, 3) ./ sqrt(2));
  outstr.cpFit_b_spline_pval(:, :, 1) = thisCalc(cpFit_b_spl);
  outstr.cpFit_m_spline_pval(:, :, 1) = thisCalc(cpFit_m_spl);
  outstr.cpFit_b_linear_pval(:, :, 1) = thisCalc(cpFit_b_lin);
  outstr.cpFit_m_linear_pval(:, :, 1) = thisCalc(cpFit_m_lin);
  outstr.cpFit_b_apc_pval(:, :, 1) = thisCalc(cpFit_b_apc);
  outstr.cpFit_m_apc_pval(:, :, 1) = thisCalc(cpFit_m_apc);
  outstr.cpFit_b_ppc_pval(:, :, 1) = thisCalc(cpFit_b_ppc);
  outstr.cpFit_m_ppc_pval(:, :, 1) = thisCalc(cpFit_m_ppc);
else
  thisCalc = @(dat) 2 * min(mean(dat > 0, 3), mean(dat < 0, 3));
  outstr.cpFit_b_spline_pval(:, :, 1) = thisCalc(cpFit_b_spl);
  outstr.cpFit_m_spline_pval(:, :, 1) = thisCalc(cpFit_m_spl);
  outstr.cpFit_b_linear_pval(:, :, 1) = thisCalc(cpFit_b_lin);
  outstr.cpFit_m_linear_pval(:, :, 1) = thisCalc(cpFit_m_lin);
  outstr.cpFit_b_apc_pval(:, :, 1) = thisCalc(cpFit_b_apc);
  outstr.cpFit_m_apc_pval(:, :, 1) = thisCalc(cpFit_m_apc);
  outstr.cpFit_b_ppc_pval(:, :, 1) = thisCalc(cpFit_b_ppc);
  outstr.cpFit_m_ppc_pval(:, :, 1) = thisCalc(cpFit_m_ppc);
end

% 3. Shuffling slices; then check percentage of what is
% different
cpFitS_b_spl = zeros(3, n_reg, BOOT);
cpFitS_b_lin = zeros(2, n_reg, BOOT);
cpFitS_b_apc = zeros(2, n_reg, BOOT);
cpFitS_b_ppc = zeros(2, n_reg, BOOT);
cpFitS_m_spl = zeros(2, n_reg, BOOT);
cpFitS_m_lin = zeros(1, n_reg, BOOT);
cpFitS_m_apc = zeros(1, n_reg, BOOT);
cpFitS_m_ppc = zeros(1, n_reg, BOOT);
for r = 1:n_reg
  for i = 1:BOOT
    s_shufApc = randperm(n_apc - SKIP);
    s_shufPpc = randperm(n_ppc - SKIP);
    s_shufPc = [s_shufApc, n_apc - SKIP + s_shufPpc];
    [cpFitS_m_spl(:, r, i), cpFitS_b_spl(:, r, i)] = aux.splineFit( ...
      s_bootPc(i, s_shufPc), outstr.conProb_pc(r, s_bootPc(i, :)), l_spl);
    [cpFitS_m_lin(:, r, i), cpFitS_b_lin(:, r, i)] = aux.splineFit( ...
      s_bootPc(i, s_shufPc), outstr.conProb_pc(r, s_bootPc(i, :)), l_lin);
    [cpFitS_m_apc(:, r, i), cpFitS_b_apc(:, r, i)] = aux.splineFit( ...
      s_bootApc(i, s_shufApc), outstr.conProb_pc(r, s_bootApc(i, :)), l_apc);
    [cpFitS_m_ppc(:, r, i), cpFitS_b_ppc(:, r, i)] = aux.splineFit( ...
      s_bootPpc(i, s_shufPpc), outstr.conProb_pc(r, s_bootPpc(i, :)), l_ppc);
  end
end
if USE_NORMALPDF
  thisCalc = @(nh, dat) 1 - erf(abs(mean((dat - nh), 3) ./ ...
    std(dat, 0, 3) ./ sqrt(2)));
  outstr.cpFit_b_spline_pval(:, :, 2) = thisCalc(outstr.cpFit_b_spline, cpFitS_b_spl);
  outstr.cpFit_m_spline_pval(:, :, 2) = thisCalc(outstr.cpFit_m_spline, cpFitS_m_spl);
  outstr.cpFit_b_linear_pval(:, :, 2) = thisCalc(outstr.cpFit_b_linear, cpFitS_b_lin);
  outstr.cpFit_m_linear_pval(:, :, 2) = thisCalc(outstr.cpFit_m_linear, cpFitS_m_lin);
  outstr.cpFit_b_apc_pval(:, :, 2) = thisCalc(outstr.cpFit_b_apc, cpFitS_b_apc);
  outstr.cpFit_m_apc_pval(:, :, 2) = thisCalc(outstr.cpFit_m_apc, cpFitS_m_apc);
  outstr.cpFit_b_ppc_pval(:, :, 2) = thisCalc(outstr.cpFit_b_ppc, cpFitS_b_ppc);
  outstr.cpFit_m_ppc_pval(:, :, 2) = thisCalc(outstr.cpFit_m_ppc, cpFitS_m_ppc);
else
  thisCalc = @(nh, dat) 2 * min( ...
    mean((dat - nh) > 0, 3), ...
    mean((dat - nh) < 0, 3));
  outstr.cpFit_b_spline_pval(:, :, 2) = thisCalc(outstr.cpFitS_b_spline, cpFitS_b_spl);
  outstr.cpFit_m_spline_pval(:, :, 2) = thisCalc(outstr.cpFit_m_spline, cpFitS_m_spl);
  outstr.cpFit_b_linear_pval(:, :, 2) = thisCalc(outstr.cpFit_b_linear, cpFitS_b_lin);
  outstr.cpFit_m_linear_pval(:, :, 2) = thisCalc(outstr.cpFit_m_linear, cpFitS_m_lin);
  outstr.cpFit_b_apc_pval(:, :, 2) = thisCalc(outstr.cpFit_b_apc, cpFitS_b_apc);
  outstr.cpFit_m_apc_pval(:, :, 2) = thisCalc(outstr.cpFit_m_apc, cpFitS_m_apc);
  outstr.cpFit_b_ppc_pval(:, :, 2) = thisCalc(outstr.cpFit_b_ppc, cpFitS_b_ppc);
  outstr.cpFit_m_ppc_pval(:, :, 2) = thisCalc(outstr.cpFit_m_ppc, cpFitS_m_ppc);
end

% Fix p-values to lower limits
outstr.cpFit_b_linear_pval(outstr.cpFit_b_linear_pval(:) < PMIN) = PMIN;
outstr.cpFit_b_spline_pval(outstr.cpFit_b_spline_pval(:) < PMIN) = PMIN;
outstr.cpFit_b_apc_pval(outstr.cpFit_b_apc_pval(:) < PMIN) = PMIN;
outstr.cpFit_b_ppc_pval(outstr.cpFit_b_ppc_pval(:) < PMIN) = PMIN;
outstr.cpFit_m_linear_pval(outstr.cpFit_m_linear_pval(:) < PMIN) = PMIN;
outstr.cpFit_m_spline_pval(outstr.cpFit_m_spline_pval(:) < PMIN) = PMIN;
outstr.cpFit_m_apc_pval(outstr.cpFit_m_apc_pval(:) < PMIN) = PMIN;
outstr.cpFit_m_ppc_pval(outstr.cpFit_m_ppc_pval(:) < PMIN) = PMIN;

% Get goodness of fit using K-fold cross validation
cpFit_linear_r2 = zeros(BOOT, CROSSVAL, n_reg);
cpFit_spline_r2 = zeros(BOOT, CROSSVAL, n_reg);
cpFit_region_r2 = zeros(BOOT, CROSSVAL, n_reg);
sto_matVald = zeros(n_reg, n_sli, CROSSVAL);
sto_matTest = zeros(n_reg, n_sli, CROSSVAL);
sto_matPredLin = zeros(n_reg, length(x_fpc), CROSSVAL);
sto_matPredSpl = zeros(n_reg, length(x_fpc), CROSSVAL);
sto_matPredReg = zeros(n_reg, length(x_fpc), CROSSVAL);
cpFit_b_spl = zeros(3, n_reg, CROSSVAL);
cpFit_b_lin = zeros(2, n_reg, CROSSVAL);
cpFit_b_apc = zeros(2, n_reg, CROSSVAL);
cpFit_b_ppc = zeros(2, n_reg, CROSSVAL);
cpFit_m_spl = zeros(2, n_reg, CROSSVAL);
cpFit_m_lin = zeros(1, n_reg, CROSSVAL);
cpFit_m_apc = zeros(1, n_reg, CROSSVAL);
cpFit_m_ppc = zeros(1, n_reg, CROSSVAL);
outstr.cpFit_boot_linR2 = zeros(BOOT, n_reg);
outstr.cpFit_boot_splR2 = zeros(BOOT, n_reg);
outstr.cpFit_boot_regR2 = zeros(BOOT, n_reg);
for b = 1:BOOT
  % Get crossvalidation indices
  sto_i = crossvalind('Kfold', n_brc, CROSSVAL);
  % Get the validation and testing set for each cross validation set
  for k = 1:CROSSVAL
    s_i = (sto_i == k);
    sto_matVald(:, :, k) = aux.conProb(regmat( s_i, :), pcmat( s_i, :));
    sto_matTest(:, :, k) = aux.conProb(regmat(~s_i, :), pcmat(~s_i, :));
  end
  % Run fits on the test sets for each dataset and each region
  for k = 1:CROSSVAL
    for r = 1:n_reg
      [cpFit_m_lin(:, r, k), cpFit_b_lin(:, r, k)] = aux.splineFit( ...
        x_fpc', sto_matTest(r, x_fpc, k), l_spl([1, 3]));
      [cpFit_m_spl(:, r, k), cpFit_b_spl(:, r, k)] = aux.splineFit( ...
        x_fpc', sto_matTest(r, x_fpc, k), l_spl(1:3));
      [cpFit_m_apc(:, r, k), cpFit_b_apc(:, r, k)] = aux.splineFit( ...
        x_apc', sto_matTest(r, x_apc, k), l_spl(1:2));
      [cpFit_m_ppc(:, r, k), cpFit_b_ppc(:, r, k)] = aux.splineFit( ...
        x_ppc', sto_matTest(r, x_ppc, k), l_spl(2:3));
      % Find m and b for fits
      sto_linmb = [l_spl(1), 1; l_spl(3), 1] \ [cpFit_b_lin(1, r, k); cpFit_b_lin(2, r, k)];
      sto_spl_a = [l_spl(1), 1; l_spl(2), 1] \ [cpFit_b_spl(1, r, k); cpFit_b_spl(2, r, k)];
      sto_spl_p = [l_spl(2), 1; l_spl(3), 1] \ [cpFit_b_spl(2, r, k); cpFit_b_spl(3, r, k)];
      sto_reg_a = [l_spl(1), 1; l_spl(2), 1] \ [cpFit_b_apc(1, r, k); cpFit_b_apc(2, r, k)];
      sto_reg_p = [l_spl(2), 1; l_spl(3), 1] \ [cpFit_b_ppc(1, r, k); cpFit_b_ppc(2, r, k)];
      % Using the m and b values; predict the matrix values
      sto_matPredLin(r, :, k) = (sto_linmb(1) * (x_fpc')) + sto_linmb(2);
      sto_matPredSpl(r, 1:length(x_apc), k) = ...
        (sto_spl_a(1) * (x_apc')) + sto_spl_a(2);
      sto_matPredSpl(r, length(x_apc) + (1:length(x_ppc)), k) = ...
        (sto_spl_p(1) * (x_ppc')) + sto_spl_p(2);
      sto_matPredReg(r, 1:length(x_apc), k) = ...
        (sto_reg_a(1) * (x_apc')) + sto_reg_a(2);
      sto_matPredReg(r, length(x_apc) + (1:length(x_ppc)), k) = ...
        (sto_reg_p(1) * (x_ppc')) + sto_reg_p(2);
      % Calculate the R2 value on the validation set
      cpFit_linear_r2(b, k, r) = sum( ...
        (sto_matPredLin(r, :, k) - sto_matVald(r, x_fpc, k)) .^ 2);
      cpFit_spline_r2(b, k, r) = sum( ...
        (sto_matPredSpl(r, :, k) - sto_matVald(r, x_fpc, k)) .^ 2);
      cpFit_region_r2(b, k, r) = sum( ...
        (sto_matPredReg(r, :, k) - sto_matVald(r, x_fpc, k)) .^ 2);
    end
  end
  % Record the R2 values for each fit
  for r = 1:n_reg
    outstr.cpFit_boot_linR2(b, r) = sum(cpFit_linear_r2(b, :, r), 'all');
    outstr.cpFit_boot_splR2(b, r) = sum(cpFit_spline_r2(b, :, r), 'all');
    outstr.cpFit_boot_regR2(b, r) = sum(cpFit_region_r2(b, :, r), 'all');
  end
end
% Get the average values for R2 for the fits
outstr.cpFit_linear_r2 = mean(outstr.cpFit_boot_linR2, 1)';
outstr.cpFit_spline_r2 = mean(outstr.cpFit_boot_splR2, 1)';
outstr.cpFit_region_r2 = mean(outstr.cpFit_boot_regR2, 1)';
% Get p value wrt the linear
s_s = std(outstr.cpFit_boot_linR2, 0, 1);
outstr.cpFit_spline_r2_pval = zeros(n_reg, 1);
outstr.cpFit_region_r2_pval = zeros(n_reg, 1);
for r = 1:n_reg
  [~, outstr.cpFit_spline_r2_pval(r)] = ztest( ...
    outstr.cpFit_boot_splR2(:, r), outstr.cpFit_linear_r2(r), s_s(r));
  [~, outstr.cpFit_region_r2_pval(r)] = ztest( ...
    outstr.cpFit_boot_regR2(:, r), outstr.cpFit_linear_r2(r), s_s(r));
end
outstr.cpFit_spline_r2_pval(outstr.cpFit_spline_r2_pval(:) < PMIN) = PMIN;
outstr.cpFit_region_r2_pval(outstr.cpFit_region_r2_pval(:) < PMIN) = PMIN;

% PCA
[~, outstr.conProb_sliPca3] = ...
    pca(outstr.conProb_pc(:, (1 + SKIP):(end - SKIP))', 'NumComponents', 3);
[~, outstr.conProb_sliPca2] = ...
    pca(outstr.conProb_pc(:, (1 + SKIP):(end - SKIP))', 'NumComponents', 2);

% Do table for easy reading
outstr.pTable = table( ...
  outstr.conProb_pc_sprCorrPval, ...
  outstr.conProb_apc_sprCorrPval, outstr.conProb_ppc_sprCorrPval, ...
  outstr.cpFit_m_spline_pval(1, :, 1)', outstr.cpFit_m_spline_pval(2, :, 1)', ...
  outstr.cpFit_m_linear_pval(1, :, 1)', ...
  outstr.cpFit_m_apc_pval(1, :, 1)', outstr.cpFit_m_ppc_pval(1, :, 1)', ...
  outstr.cpFit_m_spline_pval(1, :, 2)', outstr.cpFit_m_spline_pval(2, :, 2)', ...
  outstr.cpFit_m_linear_pval(1, :, 2)', ...
  outstr.cpFit_m_apc_pval(1, :, 2)', outstr.cpFit_m_ppc_pval(1, :, 2)', ...
  'RowNames', regions, 'VariableNames', { ...
  'Spearman Corr p', ...
  'Spearman Corr APC p', 'Spearman Corr PPC p', ...
  'BS: Piecewise APC p', 'BS: Piecewise PPC p', ...
  'BS: Linear p', ...
  'BS: APC p', 'BS: PPC p', ...
  'Shuf: Piecewise APC p', 'Shuf: Piecewise PPC p', ...
  'Shuf: Linear p', ...
  'Shuf: APC p', 'Shuf: PPC p'});

outstr.dTable = table( ...
  outstr.conProb_pc_sprCorr, ...
  outstr.conProb_apc_sprCorr, outstr.conProb_ppc_sprCorr, ...
  outstr.cpFit_m_spline(1, :)', outstr.cpFit_m_spline(2, :)', ...
  outstr.cpFit_m_linear', ...
  outstr.cpFit_m_apc', outstr.cpFit_m_ppc', ...
  'RowNames', regions, 'VariableNames', { ...
  'Spearman Corr', ...
  'Spearman Corr APC', 'Spearman Corr PPC', ...
  'Piecewise fit APC slope', 'Piecewise fit PPC slope', ...
  'Linear fit slope', ...
  'APC fit slope', 'PPC fit slope'});

outstr.oTable = table( ...
  outstr.conProb_pc_sprCorr, ...
  outstr.conProb_pc_sprCorrPval2, ...
  outstr.conProb_apc_sprCorr, ...
  outstr.conProb_apc_sprCorrPval2, ...
  outstr.conProb_ppc_sprCorr, ...
  outstr.conProb_ppc_sprCorrPval2, ...
  outstr.cpFit_m_linear', ...
  outstr.cpFit_m_linear_pval(1, :, 1)', ...
  outstr.cpFit_m_linear_pval(1, :, 2)', ...
  outstr.cpFit_linear_r2, ...
  outstr.cpFit_m_spline(1, :)', ...
  outstr.cpFit_m_spline_pval(1, :, 1)', ...
  outstr.cpFit_m_spline_pval(1, :, 2)', ...
  outstr.cpFit_m_spline(2, :)', ...
  outstr.cpFit_m_spline_pval(2, :, 1)', ...
  outstr.cpFit_m_spline_pval(2, :, 2)', ...
  outstr.cpFit_spline_r2, ...
  outstr.cpFit_spline_r2_pval, ...
  outstr.cpFit_m_apc', ...
  outstr.cpFit_m_apc_pval(1, :, 1)', ...
  outstr.cpFit_m_apc_pval(1, :, 2)', ...
  outstr.cpFit_m_ppc', ...
  outstr.cpFit_m_ppc_pval(1, :, 1)', ...
  outstr.cpFit_m_ppc_pval(1, :, 2)', ...
  outstr.cpFit_region_r2, ...
  outstr.cpFit_region_r2_pval, ...
  'RowNames', regions, 'VariableNames', { ...
  'PC Spearman Correlation', ...
  'PC Spearman Corr. p-value', ...
  'APC Spearman Correlation', ...
  'APC Spearman Corr. p-value', ...
  'PPC Spearman Correlation', ...
  'PPC Spearman Corr. p-value', ...
  'PC Linear Fit Slope', ...
  'PC Linear Fit Slope p-value (bootstrap)', ...
  'PC Linear Fit Slope p-value (shuffle)', ...
  'PC Linear Fit, mean r2 error', ...
  'Piecewise Fit, APC Slope', ...
  'Piecewise Fit, APC Slope p-value (bootstrap)', ...
  'Piecewise Fit, APC Slope p-value (shuffle)', ...
  'Piecewise Fit, PPC Slope', ...
  'Piecewise Fit, PPC Slope p-value (bootstrap)', ...
  'Piecewise Fit, PPC Slope p-value (shuffle)', ...
  'Piecewise Fit, mean r2 error', ...
  'Piecewise Fit, mean r2 error p-value', ...
  'APC Linear Fit Slope', ...
  'APC Linear Fit Slope p-value (bootstrap)', ...
  'APC Linear Fit Slope p-value (shuffle)', ...
  'PPC Linear Fit Slope', ...
  'PPC Linear Fit Slope p-value (bootstrap)', ...
  'PPC Linear Fit Slope p-value (shuffle)', ...
  'AP Linear Fits, mean r2 error', ...
  'AP Linear Fits, mean r2 error p-value'});

set.data.PCout = outstr;
rng('shuffle');

end

