function newstatOBOT(set, BOOT)
%NEWSTATOB Three things that are, ugh . . .
% 1. Spearman correlation
% 2. Bootstrap slices for fit slope; check how often the slope changes
%   signs
% 3. Shuffle slices to calculate pvalue wrt shuffling
%#ok<*UNRCH>
rng('default');
USE_NORMALPDF = true;

% Store variables here
outstr = struct();

% Collapse the regions into one; use only mitral cell data
regmat = [ ...
  sum(set.prjImg(:, set.prjRegInd{1}), 2), ...
  sum(set.prjImg(:, set.prjRegInd{2}), 2), ...
  sum(set.prjImg(:, set.prjRegInd{3}), 2), ...
  sum(set.prjImg(:, set.prjRegInd{5}), 2), ...
  sum(set.prjImg(:, set.prjRegInd{6}), 2)];
otmat = set.prjImg(:, [set.prjRegInd{4}]);

% Useful
n_brc = size(regmat, 1);
n_reg = size(regmat, 2);
n_sli = size(otmat, 2);
% Limits for fitting
l_lin = [1 - 1e-5; n_sli];

% Calculate the full OB matrix
outstr.conProb = aux.conProb(regmat, otmat);

% 1. Calculate spearman correlation between regions and PC
[outstr.conProb_sprCorr, outstr.conProb_sprCorrPval] = ...
  corr(outstr.conProb', (1:n_sli)', 'Type', 'Spearman');

% 2. Fit line using bootstrap
s_bootOt = ceil(n_sli * rand(BOOT, n_sli));
cpFit_b_lin = zeros(2, n_reg, BOOT);
cpFit_m_lin = zeros(1, n_reg, BOOT);
for r = 1:n_reg
  for i = 1:BOOT
    [cpFit_m_lin(:, r, i), cpFit_b_lin(:, r, i)] = aux.splineFit( ...
      s_bootOt(i, :), outstr.conProb(r, s_bootOt(i, :)), l_lin);
  end
end
outstr.cpFit_b_linear = mean(cpFit_b_lin, 3);
outstr.cpFit_m_linear = mean(cpFit_m_lin, 3);

outstr.cpFit_b_linear_std = std(cpFit_b_lin, 0, 3);
outstr.cpFit_m_linear_std = std(cpFit_m_lin, 0, 3);

outstr.cpFit_b_linear_pval = zeros(2, n_reg, 2);
outstr.cpFit_m_linear_pval = zeros(1, n_reg, 2);

if USE_NORMALPDF
  thisCalc = @(dat) 1 - erf(abs(mean(dat, 3)) ./ std(dat, 0, 3) ./ sqrt(2));
  outstr.cpFit_b_linear_pval(:, :, 1) = thisCalc(cpFit_b_lin);
  outstr.cpFit_m_linear_pval(:, :, 1) = thisCalc(cpFit_m_lin);
else
  thisCalc = @(dat) 2 * min(mean(dat > 0, 3), mean(dat < 0, 3));
  outstr.cpFit_b_linear_pval(:, :, 1) = thisCalc(cpFit_b_lin);
  outstr.cpFit_m_linear_pval(:, :, 1) = thisCalc(cpFit_m_lin);
end

% 3. Shuffling slices; then check percentage of what is
% different
cpFitS_b_lin = zeros(2, n_reg, BOOT);
cpFitS_m_lin = zeros(1, n_reg, BOOT);
for r = 1:n_reg
  for i = 1:BOOT
    s_shufOt = randperm(n_sli);
    [cpFitS_m_lin(:, r, i), cpFitS_b_lin(:, r, i)] = aux.splineFit( ...
      s_bootOt(i, s_shufOt), outstr.conProb(r, s_bootOt(i, :)), l_lin);
  end
end
if USE_NORMALPDF
  thisCalc = @(nh, dat) 1 - erf(abs(mean((dat - nh), 3) ./ ...
    std(dat, 0, 3) ./ sqrt(2)));
  outstr.cpFit_b_linear_pval(:, :, 2) = thisCalc(outstr.cpFit_b_linear, cpFitS_b_lin);
  outstr.cpFit_m_linear_pval(:, :, 2) = thisCalc(outstr.cpFit_m_linear, cpFitS_m_lin);
else
  thisCalc = @(nh, dat) 2 * min( ...
    mean((dat - nh) > 0, 3), ...
    mean((dat - nh) < 0, 3));
  outstr.cpFit_b_linear_pval(:, :, 2) = thisCalc(outstr.cpFit_b_linear, cpFitS_b_lin);
  outstr.cpFit_m_linear_pval(:, :, 2) = thisCalc(outstr.cpFit_m_linear, cpFitS_m_lin);
end

% Do table for easy reading
outstr.pTable = table( ...
  outstr.conProb_sprCorrPval, ...
  outstr.cpFit_m_linear_pval(1, :, 1)', ...
  outstr.cpFit_m_linear_pval(1, :, 2)', ...
  'RowNames', set.prjRegName([1, 2, 3, 5, 6]), 'VariableNames', { ...
  'Spearman Corr p', 'BS: Linear p', 'Shuffle: Linear p'});

outstr.dTable = table( ...
  outstr.conProb_sprCorr, ...
  outstr.cpFit_m_linear', ...
  'RowNames', set.prjRegName([1, 2, 3, 5, 6]), 'VariableNames', { ...
  'Spearman Corr', 'Linear fit slope'});

set.data.OBOT = outstr;
rng('shuffle');

end

