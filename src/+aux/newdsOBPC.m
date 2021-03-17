function newdsOB(set, N_D, NUMS, BOOT)
%NEWSTATOB Three things that are, ugh . . .
% 1. Spearman correlation
% 2. Just run bootstrap on the points for the fit certainty.
%#ok<*UNRCH>
rng('default');
DOWN = round(logspace(1, log10(set.nBrc), N_D));
USE_NORMALPDF = true;

% Store variables here
outstr = struct();
avgstr = struct();

% Collapse the regions into one; use only mitral cell data
regmat_base = [ ...
  sum(set.prjImg(:, set.prjRegInd{1}), 2), ...
  sum(set.prjImg(:, set.prjRegInd{4}), 2), ...
  sum(set.prjImg(:, set.prjRegInd{5}), 2), ...
  sum(set.prjImg(:, set.prjRegInd{6}), 2)];
apcmat_base = set.prjImg(:, set.prjRegInd{2});
ppcmat_base = set.prjImg(:, set.prjRegInd{3});
pcmat_base = [apcmat_base, ppcmat_base];

% Useful
n_reg = size(regmat_base, 2);
n_apc = size(apcmat_base, 2);
n_ppc = size(ppcmat_base, 2);
n_sli = size(pcmat_base, 2);
% Limits for fitting
l_spl = [1 - 1e-5; n_apc + .5; n_sli];
l_lin = [1 - 1e-5; n_sli];
l_apc = [1 - 1e-5; n_apc];
l_ppc = [1 - 1e-5; n_ppc];

% Do this for all downsampling sizes
for d = N_D:-1:1
  n_brc = DOWN(d);
  for t = NUMS:-1:1
    aux.progressbar(((N_D - d) * NUMS) + NUMS - t + 1, N_D * NUMS);
    s_this = randsample(set.nBrc, n_brc);
    regmat = regmat_base(s_this, :);
    apcmat = apcmat_base(s_this, :);
    ppcmat = ppcmat_base(s_this, :);
    pcmat =  pcmat_base( s_this, :);
    outstr(d, t).brcIds = s_this;
    
    % Calculate the full OB matrix
    outstr(d, t).conProb_pc = aux.conProb(regmat, pcmat);
    outstr(d, t).conProb_apc = outstr(d, t).conProb_pc(:, 1:n_apc);
    outstr(d, t).conProb_ppc = outstr(d, t).conProb_pc(:, n_apc + (1:n_ppc));
    
    % Convert 0's to NaN's
    conProb_pc = outstr(d, t).conProb_pc;
    conProb_pc(conProb_pc(:) == 0) = NaN;
    
    % 1. Calculate spearman correlation between regions and PC
    [outstr(d, t).conProb_pc_sprCorr, outstr(d, t).conProb_pc_sprCorrPval] = ...
      corr(conProb_pc', (1:n_sli)', 'Type', 'Spearman');
    [outstr(d, t).conProb_apc_sprCorr, outstr(d, t).conProb_apc_sprCorrPval] = ...
      corr(conProb_pc(:, 1:n_apc)', (1:n_apc)', 'Type', 'Spearman');
    [outstr(d, t).conProb_ppc_sprCorr, outstr(d, t).conProb_ppc_sprCorrPval] = ...
      corr(conProb_pc(:, n_apc + (1:n_ppc))', (1:n_ppc)', 'Type', 'Spearman');
    
    % 2. Fit several lines using bootstrap
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
          s_bootPc(i, :), outstr(d, t).conProb_pc(r, s_bootPc(i, :)), l_spl);
        [cpFit_m_lin(:, r, i), cpFit_b_lin(:, r, i)] = aux.splineFit( ...
          s_bootPc(i, :), outstr(d, t).conProb_pc(r, s_bootPc(i, :)), l_lin);
        [cpFit_m_apc(:, r, i), cpFit_b_apc(:, r, i)] = aux.splineFit( ...
          s_bootApc(i, :), outstr(d, t).conProb_apc(r, s_bootApc(i, :)), l_apc);
        [cpFit_m_ppc(:, r, i), cpFit_b_ppc(:, r, i)] = aux.splineFit( ...
          s_bootPpc(i, :), outstr(d, t).conProb_ppc(r, s_bootPpc(i, :)), l_ppc);
      end
    end
    outstr(d, t).cpFit_b_spline = mean(cpFit_b_spl, 3);
    outstr(d, t).cpFit_b_spline_std = std(cpFit_b_spl, 0, 3);
    outstr(d, t).cpFit_m_spline = mean(cpFit_m_spl, 3);
    outstr(d, t).cpFit_m_spline_std = std(cpFit_m_spl, 0, 3);
    outstr(d, t).cpFit_b_linear = mean(cpFit_b_lin, 3);
    outstr(d, t).cpFit_b_linear_std = std(cpFit_b_lin, 0, 3);
    outstr(d, t).cpFit_m_linear = mean(cpFit_m_lin, 3);
    outstr(d, t).cpFit_m_linear_std = std(cpFit_m_lin, 0, 3);
    outstr(d, t).cpFit_b_apc = mean(cpFit_b_apc, 3);
    outstr(d, t).cpFit_b_apc_std = std(cpFit_b_apc, 0, 3);
    outstr(d, t).cpFit_m_apc = mean(cpFit_m_apc, 3);
    outstr(d, t).cpFit_m_apc_std = std(cpFit_m_apc, 0, 3);
    outstr(d, t).cpFit_b_ppc = mean(cpFit_b_ppc, 3);
    outstr(d, t).cpFit_b_ppc_std = std(cpFit_b_ppc, 0, 3);
    outstr(d, t).cpFit_m_ppc = mean(cpFit_m_ppc, 3);
    outstr(d, t).cpFit_m_ppc_std = std(cpFit_m_ppc, 0, 3);
    
    % Initialize pvalue storage
    outstr(d, t).cpFit_b_spline_pval = zeros(3, n_reg, 2);
    outstr(d, t).cpFit_m_spline_pval = zeros(2, n_reg, 2);
    outstr(d, t).cpFit_b_linear_pval = zeros(2, n_reg, 2);
    outstr(d, t).cpFit_m_linear_pval = zeros(1, n_reg, 2);
    outstr(d, t).cpFit_b_apc_pval = zeros(2, n_reg, 2);
    outstr(d, t).cpFit_m_apc_pval = zeros(1, n_reg, 2);
    outstr(d, t).cpFit_b_ppc_pval = zeros(2, n_reg, 2);
    outstr(d, t).cpFit_m_ppc_pval = zeros(1, n_reg, 2);
    
    if USE_NORMALPDF
      thisCalc = @(dat) 1 - erf(abs(mean(dat, 3)) ...
        ./ std(dat, 0, 3) ./ sqrt(2));
      outstr(d, t).cpFit_b_spline_pval(:, :, 1) = thisCalc(cpFit_b_spl);
      outstr(d, t).cpFit_m_spline_pval(:, :, 1) = thisCalc(cpFit_m_spl);
      outstr(d, t).cpFit_b_linear_pval(:, :, 1) = thisCalc(cpFit_b_lin);
      outstr(d, t).cpFit_m_linear_pval(:, :, 1) = thisCalc(cpFit_m_lin);
      outstr(d, t).cpFit_b_apc_pval(:, :, 1) = thisCalc(cpFit_b_apc);
      outstr(d, t).cpFit_m_apc_pval(:, :, 1) = thisCalc(cpFit_m_apc);
      outstr(d, t).cpFit_b_ppc_pval(:, :, 1) = thisCalc(cpFit_b_ppc);
      outstr(d, t).cpFit_m_ppc_pval(:, :, 1) = thisCalc(cpFit_m_ppc);
    else
      thisCalc = @(dat) 2 * min(mean(dat > 0, 3), mean(dat < 0, 3));
      outstr(d, t).cpFit_b_spline_pval(:, :, 1) = thisCalc(cpFit_b_spl);
      outstr(d, t).cpFit_m_spline_pval(:, :, 1) = thisCalc(cpFit_m_spl);
      outstr(d, t).cpFit_b_linear_pval(:, :, 1) = thisCalc(cpFit_b_lin);
      outstr(d, t).cpFit_m_linear_pval(:, :, 1) = thisCalc(cpFit_m_lin);
      outstr(d, t).cpFit_b_apc_pval(:, :, 1) = thisCalc(cpFit_b_apc);
      outstr(d, t).cpFit_m_apc_pval(:, :, 1) = thisCalc(cpFit_m_apc);
      outstr(d, t).cpFit_b_ppc_pval(:, :, 1) = thisCalc(cpFit_b_ppc);
      outstr(d, t).cpFit_m_ppc_pval(:, :, 1) = thisCalc(cpFit_m_ppc);
    end
    
    % 3. Do the same shuffling slices; then check percentage of what is
    % different
    for r = 1:n_reg
      for i = 1:BOOT
        s_shufApc = randperm(n_apc);
        s_shufPpc = randperm(n_ppc);
        s_shufPc = randperm(n_apc + n_ppc);
        [cpFit_m_spl(:, r, i), cpFit_b_spl(:, r, i)] = aux.splineFit( ...
          s_bootPc(i, s_shufPc), outstr(d, t).conProb_pc(r, s_bootPc(i, :)), l_spl);
        [cpFit_m_lin(:, r, i), cpFit_b_lin(:, r, i)] = aux.splineFit( ...
          s_bootPc(i, s_shufPc), outstr(d, t).conProb_pc(r, s_bootPc(i, :)), l_lin);
        [cpFit_m_apc(:, r, i), cpFit_b_apc(:, r, i)] = aux.splineFit( ...
          s_bootApc(i, s_shufApc), outstr(d, t).conProb_apc(r, s_bootApc(i, :)), l_apc);
        [cpFit_m_ppc(:, r, i), cpFit_b_ppc(:, r, i)] = aux.splineFit( ...
          s_bootPpc(i, s_shufPpc), outstr(d, t).conProb_ppc(r, s_bootPpc(i, :)), l_ppc);
      end
    end
    
    % Calculate p-values again
    if USE_NORMALPDF
      thisCalc = @(nh, dat) 1 - erf(abs(mean((dat - nh), 3) ./ ...
        std(dat, 0, 3) ./ sqrt(2)));
      outstr(d, t).cpFit_b_spline_pval(:, :, 2) = thisCalc(outstr(d, t).cpFit_b_spline, cpFit_b_spl);
      outstr(d, t).cpFit_m_spline_pval(:, :, 2) = thisCalc(outstr(d, t).cpFit_m_spline, cpFit_m_spl);
      outstr(d, t).cpFit_b_linear_pval(:, :, 2) = thisCalc(outstr(d, t).cpFit_b_linear, cpFit_b_lin);
      outstr(d, t).cpFit_m_linear_pval(:, :, 2) = thisCalc(outstr(d, t).cpFit_m_linear, cpFit_m_lin);
      outstr(d, t).cpFit_b_apc_pval(:, :, 2) = thisCalc(outstr(d, t).cpFit_b_apc, cpFit_b_apc);
      outstr(d, t).cpFit_m_apc_pval(:, :, 2) = thisCalc(outstr(d, t).cpFit_m_apc, cpFit_m_apc);
      outstr(d, t).cpFit_b_ppc_pval(:, :, 2) = thisCalc(outstr(d, t).cpFit_b_ppc, cpFit_b_ppc);
      outstr(d, t).cpFit_m_ppc_pval(:, :, 2) = thisCalc(outstr(d, t).cpFit_m_ppc, cpFit_m_ppc);
    else
      thisCalc = @(nh, dat) 2 * min(mean((dat - nh) > 0, 3), ...
        mean((dat - nh) < 0, 3));
      outstr(d, t).cpFit_b_spline_pval(:, :, 2) = thisCalc(outstr(d, t).cpFit_b_spline, cpFit_b_spl);
      outstr(d, t).cpFit_m_spline_pval(:, :, 2) = thisCalc(outstr(d, t).cpFit_m_spline, cpFit_m_spl);
      outstr(d, t).cpFit_b_linear_pval(:, :, 2) = thisCalc(outstr(d, t).cpFit_b_linear, cpFit_b_lin);
      outstr(d, t).cpFit_m_linear_pval(:, :, 2) = thisCalc(outstr(d, t).cpFit_m_linear, cpFit_m_lin);
      outstr(d, t).cpFit_b_apc_pval(:, :, 2) = thisCalc(outstr(d, t).cpFit_b_apc, cpFit_b_apc);
      outstr(d, t).cpFit_m_apc_pval(:, :, 2) = thisCalc(outstr(d, t).cpFit_m_apc, cpFit_m_apc);
      outstr(d, t).cpFit_b_ppc_pval(:, :, 2) = thisCalc(outstr(d, t).cpFit_b_ppc, cpFit_b_ppc);
      outstr(d, t).cpFit_m_ppc_pval(:, :, 2) = thisCalc(outstr(d, t).cpFit_m_ppc, cpFit_m_ppc);
    end
    
    % Do table for easy reading
    outstr(d, t).pTable = table( ...
      outstr(d, t).conProb_pc_sprCorrPval, ...
      outstr(d, t).conProb_apc_sprCorrPval, outstr(d, t).conProb_ppc_sprCorrPval, ...
      outstr(d, t).cpFit_m_spline_pval(1, :, 1)', outstr(d, t).cpFit_m_spline_pval(2, :, 1)', ...
      outstr(d, t).cpFit_m_linear_pval(1, :, 1)', ...
      outstr(d, t).cpFit_m_apc_pval(1, :, 1)', outstr(d, t).cpFit_m_ppc_pval(1, :, 1)', ...
      outstr(d, t).cpFit_m_spline_pval(1, :, 2)', outstr(d, t).cpFit_m_spline_pval(2, :, 2)', ...
      outstr(d, t).cpFit_m_linear_pval(1, :, 2)', ...
      outstr(d, t).cpFit_m_apc_pval(1, :, 2)', outstr(d, t).cpFit_m_ppc_pval(1, :, 2)', ...
      'RowNames', set.prjRegName([1, 4, 5, 6]), 'VariableNames', { ...
      'Spearman Corr p', ...
      'Spearman Corr APC p', 'Spearman Corr PPC p', ...
      'BS: Piecewise APC p', 'BS: Piecewise PPC p', ...
      'BS: Linear p', ...
      'BS: APC p', 'BS: PPC p', ...
      'Shuf: Piecewise APC p', 'Shuf: Piecewise PPC p', ...
      'Shuf: Linear p', ...
      'Shuf: APC p', 'Shuf: PPC p'});
    
    outstr(d, t).dTable = table( ...
      outstr(d, t).conProb_pc_sprCorr, ...
      outstr(d, t).conProb_apc_sprCorr, outstr(d, t).conProb_ppc_sprCorr, ...
      outstr(d, t).cpFit_m_spline(1, :)', outstr(d, t).cpFit_m_spline(2, :)', ...
      outstr(d, t).cpFit_m_linear', ...
      outstr(d, t).cpFit_m_apc', outstr(d, t).cpFit_m_ppc', ...
      'RowNames', set.prjRegName([1,  4, 5, 6]), 'VariableNames', { ...
      'Spearman Corr', ...
      'Spearman Corr APC', 'Spearman Corr PPC', ...
      'Piecewise fit APC slope', 'Piecewise fit PPC slope', ...
      'Linear fit slope', ...
      'APC fit slope', 'PPC fit slope'});
  end
end

% Get stats
avgstr.sampleN = zeros(N_D, n_reg);
% Spearman pval
avgstr.spr_lin = zeros(N_D, NUMS, n_reg);
avgstr.spr_apc = zeros(N_D, NUMS, n_reg);
avgstr.spr_ppc = zeros(N_D, NUMS, n_reg);
% Bootstrap method
avgstr.bst_lin = zeros(N_D, NUMS, n_reg);
avgstr.bst_pwa = zeros(N_D, NUMS, n_reg);
avgstr.bst_pwp = zeros(N_D, NUMS, n_reg);
avgstr.bst_apc = zeros(N_D, NUMS, n_reg);
avgstr.bst_ppc = zeros(N_D, NUMS, n_reg);
% Shuffle method
avgstr.shf_lin = zeros(N_D, NUMS, n_reg);
avgstr.shf_pwa = zeros(N_D, NUMS, n_reg);
avgstr.shf_pwp = zeros(N_D, NUMS, n_reg);
avgstr.shf_apc = zeros(N_D, NUMS, n_reg);
avgstr.shf_ppc = zeros(N_D, NUMS, n_reg);
for d = 1:N_D
  for r = 1:n_reg
    avgstr.sampleN(d, r) = DOWN(d);
    for t = 1:NUMS
      % Spearman correlation p values
      avgstr.spr_lin(d, t, r) = outstr(d, t).conProb_pc_sprCorrPval(r);
      avgstr.spr_apc(d, t, r) = outstr(d, t).conProb_apc_sprCorrPval(r);
      avgstr.spr_ppc(d, t, r) = outstr(d, t).conProb_ppc_sprCorrPval(r);
      % Bootstrap p values
      avgstr.bst_lin(d, t, r) = outstr(d, t).cpFit_m_linear_pval(1, r, 1);
      avgstr.bst_pwa(d, t, r) = outstr(d, t).cpFit_m_spline_pval(1, r, 1);
      avgstr.bst_pwp(d, t, r) = outstr(d, t).cpFit_m_spline_pval(2, r, 1);
      avgstr.bst_apc(d, t, r) = outstr(d, t).cpFit_m_apc_pval(1, r, 1);
      avgstr.bst_ppc(d, t, r) = outstr(d, t).cpFit_m_ppc_pval(1, r, 1);
      % Piecewise fit p values, use bootstrap
      avgstr.shf_lin(d, t, r) = outstr(d, t).cpFit_m_linear_pval(1, r, 2);
      avgstr.shf_pwa(d, t, r) = outstr(d, t).cpFit_m_spline_pval(1, r, 2);
      avgstr.shf_pwp(d, t, r) = outstr(d, t).cpFit_m_spline_pval(2, r, 2);
      avgstr.shf_apc(d, t, r) = outstr(d, t).cpFit_m_apc_pval(1, r, 2);
      avgstr.shf_ppc(d, t, r) = outstr(d, t).cpFit_m_ppc_pval(1, r, 2);
    end
  end
end

set.data.OBPC_dsRaw = outstr;
set.data.OBPC_dsPvl = avgstr;
rng('shuffle');

end

