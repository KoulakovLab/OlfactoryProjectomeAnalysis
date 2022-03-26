% Figure 2 is me; OB output in slice-based
%   a)It includes The slice overview (two alternate views)
%   * The conditional probability in;
%       b) Matrix form
%       c) Plot form
%       d) Subspace form
%     The conditional probability is between the slices in AP axis of PC
%     and between target regions
% Supfig 2 is me
%   * Seperation of conditional probability into minibatches
%       - yc61  + yc65
%       - yc86  + yc92
%       - yc109 + yc111
%   * Conditional probability of all cells (compared to mitral only)
%   * Conditional probability; of tufted cells; but with OT instead of PC
figure2 = gobjects(10, 1);
CONFINT = .95;
STDAMNT = norminv(.5 * (1 + CONFINT));

sto_sn = obdata_mitral.nPrjRegSli(2) + obdata_mitral.nPrjRegSli(3);
sto_apc = 1:obdata_mitral.nPrjRegSli(2);
sto_ppc = obdata_mitral.nPrjRegSli(2) + (1:obdata_mitral.nPrjRegSli(3));
sto_regName = {'AON', 'OT', 'CoA', 'lENT'};

%--------------%
%---Overview---%
%--------------%

% Create the first overview plots
% Ordered by distance
figure2(1) = figure;
set(figure2(1), 'name', 'Overview');

% Distance
obdata.plotImg('prj');
title('Barcodes sorted by brightest projection');
colormap('parula');
colorbar;
set(gca, 'fontsize', 20);

% Ordered by brightests slice
figure2(2) = figure;
set(figure2(2), 'name', 'Overview');
subplot(1, 3, 1);
obdata_mitral.plotImg('prj');
colormap('parula');
subplot(1, 3, 2);
obdata_tufted.plotImg('prj');
colormap('parula');
subplot(1, 3, 3);
obdata_deep.plotImg('prj');
colormap('parula');
sgtitle('Barcodes sorted by brightest projection');
set(gca, 'fontsize', 20);

savefig(figure2(1:2), 'data/figures/obInjection_overview.fig');

%-----------------------------%
%---Conditional Probability---%
%-----------------------------%

% FIGURE 2B: The matrix
figure2(3) = figure;
set(figure2(3), 'name', 'Cond. Prob.');

% Image
subplot(1, 2, 1);
imagesc(obdata_mitral.data.OBPC.conProb_pc);
colormap(jet);
title('P_{OB \rightarrow PC}(OB \rightarrow region) (Mitral)');
xlabel('PC slice');
% Draw APC-PPC line
nSli = obdata_mitral.nPrjRegSli(2:3);
line((nSli(1) + .5) * ones(1, 2), ylim, 'LineWidth', 3, ...
  'Color', .5 * ones(1, 3));
xticks(cumsum([0, nSli(1:(end-1))]) + .5 * nSli);
xticklabels({'APC', 'PPC'});
ylabel('Region');
yticks(1:4);
yticklabels(obdata_mitral.prjRegName([1, 4, 5, 6]));
B = colorbar;
title(B, 'P_{OB \rightarrow PC}(OB \rightarrow reg.)');
set(gca, 'fontsize', 20);

% Plot
subplot(1, 2, 2);
aux.obmatPlotNEW(obdata_mitral);

% FIGURE 2C: Same with different populations
figure2(4) = figure;
set(figure2(4), 'name', 'Cond. Prob.');

% Image
subplot(2, 2, 1);
imagesc(obdata_lowipr.data.OBPC.conProb_pc);
colormap(jet);
xlabel('PC slice');
% Draw APC-PPC line
nSli = obdata_lowipr.nPrjRegSli(2:3);
line((nSli(1) + .5) * ones(1, 2), ylim, 'LineWidth', 3, ...
  'Color', .5 * ones(1, 3));
xticks(cumsum([0, nSli(1:(end-1))]) + .5 * nSli);
xticklabels({'APC', 'PPC'});
ylabel('Region');
yticks(1:4);
yticklabels(obdata_lowipr.prjRegName([1, 4, 5, 6]));
B = colorbar;
title('Narrow projecting');
title(B, 'P(reg|pc)');
set(gca, 'fontsize', 15);

subplot(2, 2, 3);
imagesc(obdata_highipr.data.OBPC.conProb_pc);
colormap(jet);
xlabel('PC slice');
% Draw APC-PPC line
nSli = obdata_highipr.nPrjRegSli(2:3);
line((nSli(1) + .5) * ones(1, 2), ylim, 'LineWidth', 3, ...
  'Color', .5 * ones(1, 3));
xticks(cumsum([0, nSli(1:(end-1))]) + .5 * nSli);
xticklabels({'APC', 'PPC'});
ylabel('Region');
yticks(1:4);
yticklabels(obdata_highipr.prjRegName([1, 4, 5, 6]));
B = colorbar;
title('Broad projecting');
title(B, 'P(reg|pc)');
set(gca, 'fontsize', 15);

% Plot
subplot(2, 2, 2);
aux.obmatPlotNEW(obdata_lowipr);
title('Narrow projecting');
subplot(2, 2, 4);
aux.obmatPlotNEW(obdata_highipr);
title('Broad projecting');

savefig(figure2(3:4), 'data/figures/obInjection_condProb.fig');

%------------------%
%---Downsampling---%
%------------------%
thisset = [obdata_mitral, obdata_lowipr, obdata_highipr];
thisreg = {'PC', 'APC', 'PPC'};
thistype = {'Mitral', 'Narrow', 'Broad'};
s_i = 4;
s_lw = .5;
s_ms = 7.5;
s_cols = lines(4);

% Downsampling; using the spearman correlation values
figure2(5) = figure;
set(figure2(5), 'name', 'DS: Spearman');

for s = 1:length(thisset)
  ss = thisset(s);
  x_lim = [min(ss.data.OBPC_dsPvl.sampleN, [], 'all'), ...
    max(ss.data.OBPC_dsPvl.sampleN, [], 'all')];
  % Do the PC plot
  subplot(3, 3, 1 + 3 * (s - 1));
  plot(x_lim, ones(1, 2) * .05, ':k', 'LineWidth', 2);
  hold('on');
  for r = 1:4
    yval = mean(ss.data.OBPC_dsPvl.spr_lin(:, :, r), 2, 'omitnan');
    yval(yval(:) <= 0) = min(yval(yval(:) > 0), [], 'all') / 10;
    yerr = std( ss.data.OBPC_dsPvl.spr_lin(:, :, r), 0, 2, 'omitnan');
    errorbar(ss.data.OBPC_dsPvl.sampleN(:, r), yval, yerr, ...
      '-d', 'Color', s_cols(r, :), 'CapSize', s_ms, 'MarkerSize', s_ms, ...
      'MarkerFaceColor', s_cols(r, :), 'MarkerEdgeColor', 'k', ...
      'LineWidth', s_lw);
  end
  hold('off');
  sto_leg = legend({'p=0.05', 'AON', 'OT', 'CoA', 'lENT'}, 'Location', 'southwest');
  title(sto_leg, 'Regions');
  set(gca,'xscale','log');
  set(gca,'yscale','log');
  xlim(x_lim);
  xlabel('Sample size');
  ylabel('P value');
  title([thisreg{1}, ' fit, ', thistype{s}]);
  % Do the APC plot
  subplot(3, 3, 2 + 3 * (s - 1));
  plot(x_lim, ones(1, 2) * .05, ':k', 'LineWidth', 2);
  hold('on');
  for r = 1:4
    yval = mean(ss.data.OBPC_dsPvl.spr_apc(:, :, r), 2, 'omitnan');
    yval(yval(:) <= 0) = min(yval(yval(:) > 0), [], 'all') / 10;
    yerr = std( ss.data.OBPC_dsPvl.spr_apc(:, :, r), 0, 2, 'omitnan');
    errorbar(ss.data.OBPC_dsPvl.sampleN(:, r), yval, yerr, ...
      '-d', 'Color', s_cols(r, :), 'CapSize', s_ms, 'MarkerSize', s_ms, ...
      'MarkerFaceColor', s_cols(r, :), 'MarkerEdgeColor', 'k', ...
      'LineWidth', s_lw);
  end
  hold('off');
  sto_leg = legend({'p=0.05', 'AON', 'OT', 'CoA', 'lENT'}, 'Location', 'southwest');
  title(sto_leg, 'Regions');
  set(gca,'xscale','log');
  set(gca,'yscale','log');
  xlim(x_lim);
  xlabel('Sample size');
  ylabel('P value');
  title([thisreg{2}, ' fit, ', thistype{s}]);
  % Do the PPC plot
  subplot(3, 3, 3 + 3 * (s - 1));
  plot(x_lim, ones(1, 2) * .05, ':k', 'LineWidth', 2);
  hold('on');
  for r = 1:4
    yval = mean(ss.data.OBPC_dsPvl.spr_ppc(:, :, r), 2, 'omitnan');
    yval(yval(:) <= 0) = min(yval(yval(:) > 0), [], 'all') / 10;
    yerr = std( ss.data.OBPC_dsPvl.spr_ppc(:, :, r), 0, 2, 'omitnan');
    errorbar(ss.data.OBPC_dsPvl.sampleN(:, r), yval, yerr, ...
      '-d', 'Color', s_cols(r, :), 'CapSize', s_ms, 'MarkerSize', s_ms, ...
      'MarkerFaceColor', s_cols(r, :), 'MarkerEdgeColor', 'k', ...
      'LineWidth', s_lw);
  end
  hold('off');
  sto_leg = legend({'p=0.05', 'AON', 'OT', 'CoA', 'lENT'}, 'Location', 'southwest');
  title(sto_leg, 'Regions');
  set(gca,'xscale','log');
  set(gca,'yscale','log');
  xlim(x_lim);
  xlabel('Sample size');
  ylabel('P value');
  title([thisreg{3}, ' fit, ', thistype{s}]);
end
sgtitle('Downsampling: Spearman Correlation');

% Downsampling; using the shuffle method
figure2(6) = figure;
set(figure2(6), 'name', 'DS: Shuffle');

for s = 1:length(thisset)
  ss = thisset(s);
  x_lim = [min(ss.data.OBPC_dsPvl.sampleN, [], 'all'), ...
    max(ss.data.OBPC_dsPvl.sampleN, [], 'all')];
  % Do the PC plot
  subplot(3, 3, 1 + 3 * (s - 1));
  plot(x_lim, ones(1, 2) * .05, ':k', 'LineWidth', 2);
  hold('on');
  for r = 1:4
    yval = mean(ss.data.OBPC_dsPvl.shf_lin(:, :, r), 2, 'omitnan');
    yval(yval(:) <= 0) = min(yval(yval(:) > 0), [], 'all') / 10;
    yerr = std( ss.data.OBPC_dsPvl.shf_lin(:, :, r), 0, 2, 'omitnan');
    errorbar(ss.data.OBPC_dsPvl.sampleN(:, r), yval, yerr, ...
      '-d', 'Color', s_cols(r, :), 'CapSize', s_ms, 'MarkerSize', s_ms, ...
      'MarkerFaceColor', s_cols(r, :), 'MarkerEdgeColor', 'k', ...
      'LineWidth', s_lw);
  end
  hold('off');
  sto_leg = legend({'p=0.05', 'AON', 'OT', 'CoA', 'lENT'}, 'Location', 'southwest');
  title(sto_leg, 'Regions');
  set(gca,'xscale','log');
  set(gca,'yscale','log');
  xlim(x_lim);
  xlabel('Sample size');
  ylabel('P value');
  title([thisreg{1}, ' fit, ', thistype{s}]);
  % Do the APC plot
  subplot(3, 3, 2 + 3 * (s - 1));
  plot(x_lim, ones(1, 2) * .05, ':k', 'LineWidth', 2);
  hold('on');
  for r = 1:4
    yval = mean(ss.data.OBPC_dsPvl.shf_apc(:, :, r), 2, 'omitnan');
    yval(yval(:) <= 0) = min(yval(yval(:) > 0), [], 'all') / 10;
    yerr = std( ss.data.OBPC_dsPvl.shf_apc(:, :, r), 0, 2, 'omitnan');
    errorbar(ss.data.OBPC_dsPvl.sampleN(:, r), yval, yerr, ...
      '-d', 'Color', s_cols(r, :), 'CapSize', s_ms, 'MarkerSize', s_ms, ...
      'MarkerFaceColor', s_cols(r, :), 'MarkerEdgeColor', 'k', ...
      'LineWidth', s_lw);
    yval = mean(ss.data.OBPC_dsPvl.shf_pwa(:, :, r), 2, 'omitnan');
    yval(yval(:) <= 0) = min(yval(yval(:) > 0), [], 'all') / 10;
    yerr = std( ss.data.OBPC_dsPvl.shf_pwa(:, :, r), 0, 2, 'omitnan');
    errorbar(ss.data.OBPC_dsPvl.sampleN(:, r), yval, yerr, ...
      ':s', 'Color', s_cols(r, :), 'CapSize', s_ms, 'MarkerSize', s_ms, ...
      'MarkerFaceColor', s_cols(r, :), 'MarkerEdgeColor', 'k', ...
      'LineWidth', s_lw);
  end
  hold('off');
  sto_leg = legend({'p=0.05', 'AON', 'AON (p.w.)', 'OT', 'OT (p.w.)', ...
    'CoA', 'CoA (p.w.)', 'lENT', 'lENT (p.w.)'}, ...
    'Location', 'southwest', 'NumColumns', 2);
  title(sto_leg, 'Regions');
  set(gca,'xscale','log');
  set(gca,'yscale','log');
  xlim(x_lim);
  xlabel('Sample size');
  ylabel('P value');
  title([thisreg{2}, ' fit, ', thistype{s}]);
  % Do the PPC plot
  subplot(3, 3, 3 + 3 * (s - 1));
  plot(x_lim, ones(1, 2) * .05, ':k', 'LineWidth', 2);
  hold('on');
  for r = 1:4
    yval = mean(ss.data.OBPC_dsPvl.shf_ppc(:, :, r), 2, 'omitnan');
    yval(yval(:) <= 0) = min(yval(yval(:) > 0), [], 'all') / 10;
    yerr = std( ss.data.OBPC_dsPvl.shf_ppc(:, :, r), 0, 2, 'omitnan');
    errorbar(ss.data.OBPC_dsPvl.sampleN(:, r), yval, yerr, ...
      '-d', 'Color', s_cols(r, :), 'CapSize', s_ms, 'MarkerSize', s_ms, ...
      'MarkerFaceColor', s_cols(r, :), 'MarkerEdgeColor', 'k', ...
      'LineWidth', s_lw);
    yval = mean(ss.data.OBPC_dsPvl.shf_pwp(:, :, r), 2, 'omitnan');
    yval(yval(:) <= 0) = min(yval(yval(:) > 0), [], 'all') / 10;
    yerr = std( ss.data.OBPC_dsPvl.shf_pwp(:, :, r), 0, 2, 'omitnan');
    errorbar(ss.data.OBPC_dsPvl.sampleN(:, r), yval, yerr, ...
      ':s', 'Color', s_cols(r, :), 'CapSize', s_ms, 'MarkerSize', s_ms, ...
      'MarkerFaceColor', s_cols(r, :), 'MarkerEdgeColor', 'k', ...
      'LineWidth', s_lw);
  end
  hold('off');
  sto_leg = legend({'p=0.05', 'AON', 'AON (p.w.)', 'OT', 'OT (p.w.)', ...
    'CoA', 'CoA (p.w.)', 'lENT', 'lENT (p.w.)'}, ...
    'Location', 'southwest', 'NumColumns', 2);
  title(sto_leg, 'Regions');
  set(gca,'xscale','log');
  set(gca,'yscale','log');
  xlim(x_lim);
  xlabel('Sample size');
  ylabel('P value');
  title([thisreg{3}, ' fit, ', thistype{s}]);
end
sgtitle('Downsampling: Shuffle p-value');

% Downsampling; using the bootstrap method
figure2(7) = figure;
set(figure2(7), 'name', 'DS: Bootstrap');

for s = 1:length(thisset)
  ss = thisset(s);
  x_lim = [min(ss.data.OBPC_dsPvl.sampleN, [], 'all'), ...
    max(ss.data.OBPC_dsPvl.sampleN, [], 'all')];
  % Do the PC plot
  subplot(3, 3, 1 + 3 * (s - 1));
  plot(x_lim, ones(1, 2) * .05, ':k', 'LineWidth', 2);
  hold('on');
  for r = 1:4
    yval = mean(ss.data.OBPC_dsPvl.bst_lin(:, :, r), 2, 'omitnan');
    yval(yval(:) <= 0) = min(yval(yval(:) > 0), [], 'all') / 10;
    yerr = std( ss.data.OBPC_dsPvl.bst_lin(:, :, r), 0, 2, 'omitnan');
    errorbar(ss.data.OBPC_dsPvl.sampleN(:, r), yval, yerr, ...
      '-d', 'Color', s_cols(r, :), 'CapSize', s_ms, 'MarkerSize', s_ms, ...
      'MarkerFaceColor', s_cols(r, :), 'MarkerEdgeColor', 'k', ...
      'LineWidth', s_lw);
  end
  hold('off');
  sto_leg = legend({'p=0.05', 'AON', 'OT', 'CoA', 'lENT'}, 'Location', 'southwest');
  title(sto_leg, 'Regions');
  set(gca,'xscale','log');
  set(gca,'yscale','log');
  xlim(x_lim);
  xlabel('Sample size');
  ylabel('P value');
  title([thisreg{1}, ' fit, ', thistype{s}]);
  % Do the APC plot
  subplot(3, 3, 2 + 3 * (s - 1));
  plot(x_lim, ones(1, 2) * .05, ':k', 'LineWidth', 2);
  hold('on');
  for r = 1:4
    yval = mean(ss.data.OBPC_dsPvl.bst_apc(:, :, r), 2, 'omitnan');
    yval(yval(:) <= 0) = min(yval(yval(:) > 0), [], 'all') / 10;
    yerr = std( ss.data.OBPC_dsPvl.bst_apc(:, :, r), 0, 2, 'omitnan');
    errorbar(ss.data.OBPC_dsPvl.sampleN(:, r), yval, yerr, ...
      '-d', 'Color', s_cols(r, :), 'CapSize', s_ms, 'MarkerSize', s_ms, ...
      'MarkerFaceColor', s_cols(r, :), 'MarkerEdgeColor', 'k', ...
      'LineWidth', s_lw);
    yval = mean(ss.data.OBPC_dsPvl.bst_pwa(:, :, r), 2, 'omitnan');
    yval(yval(:) <= 0) = min(yval(yval(:) > 0), [], 'all') / 10;
    yerr = std( ss.data.OBPC_dsPvl.bst_pwa(:, :, r), 0, 2, 'omitnan');
    errorbar(ss.data.OBPC_dsPvl.sampleN(:, r), yval, yerr, ...
      ':s', 'Color', s_cols(r, :), 'CapSize', s_ms, 'MarkerSize', s_ms, ...
      'MarkerFaceColor', s_cols(r, :), 'MarkerEdgeColor', 'k', ...
      'LineWidth', s_lw);
  end
  hold('off');
  sto_leg = legend({'p=0.05', 'AON', 'AON (p.w.)', 'OT', 'OT (p.w.)', ...
    'CoA', 'CoA (p.w.)', 'lENT', 'lENT (p.w.)'}, ...
    'Location', 'southwest', 'NumColumns', 2);
  title(sto_leg, 'Regions');
  set(gca,'xscale','log');
  set(gca,'yscale','log');
  xlim(x_lim);
  xlabel('Sample size');
  ylabel('P value');
  title([thisreg{2}, ' fit, ', thistype{s}]);
  % Do the PPC plot
  subplot(3, 3, 3 + 3 * (s - 1));
  plot(x_lim, ones(1, 2) * .05, ':k', 'LineWidth', 2);
  hold('on');
  for r = 1:4
    yval = mean(ss.data.OBPC_dsPvl.bst_ppc(:, :, r), 2, 'omitnan');
    yval(yval(:) <= 0) = min(yval(yval(:) > 0), [], 'all') / 10;
    yerr = std( ss.data.OBPC_dsPvl.bst_ppc(:, :, r), 0, 2, 'omitnan');
    errorbar(ss.data.OBPC_dsPvl.sampleN(:, r), yval, yerr, ...
      '-d', 'Color', s_cols(r, :), 'CapSize', s_ms, 'MarkerSize', s_ms, ...
      'MarkerFaceColor', s_cols(r, :), 'MarkerEdgeColor', 'k', ...
      'LineWidth', s_lw);
    yval = mean(ss.data.OBPC_dsPvl.bst_pwp(:, :, r), 2, 'omitnan');
    yval(yval(:) <= 0) = min(yval(yval(:) > 0), [], 'all') / 10;
    yerr = std( ss.data.OBPC_dsPvl.bst_pwp(:, :, r), 0, 2, 'omitnan');
    errorbar(ss.data.OBPC_dsPvl.sampleN(:, r), yval, yerr, ...
      ':s', 'Color', s_cols(r, :), 'CapSize', s_ms, 'MarkerSize', s_ms, ...
      'MarkerFaceColor', s_cols(r, :), 'MarkerEdgeColor', 'k', ...
      'LineWidth', s_lw);
  end
  hold('off');
  sto_leg = legend({'p=0.05', 'AON', 'AON (p.w.)', 'OT', 'OT (p.w.)', ...
    'CoA', 'CoA (p.w.)', 'lENT', 'lENT (p.w.)'}, ...
    'Location', 'southwest', 'NumColumns', 2);
  title(sto_leg, 'Regions');
  set(gca,'xscale','log');
  set(gca,'yscale','log');
  xlim(x_lim);
  xlabel('Sample size');
  ylabel('P value');
  title([thisreg{3}, ' fit, ', thistype{s}]);
end
sgtitle('Downsampling: Bootstrap p-value');

savefig(figure2(5:7), 'data/figures/obInjection_downsampling.fig');

% FIGURE 2B: The matrix
figure2(8) = figure;
set(figure2(8), 'name', 'OT Cond. Prob.');

% Image
subplot(1, 2, 1);
imagesc(obdata_mitral.data.OBOT.conProb);
colormap(jet);
title('P_{OB \rightarrow OT}(OB \rightarrow region) (Mitral)');
xlabel('OT slice');
xlabel('OT slice');
ylabel('Region');
yticks(1:5);
yticklabels(obdata_mitral.prjRegName([1, 2, 3, 5, 6]));
B = colorbar;
title(B, 'P_{OB \rightarrow OT}(OB \rightarrow reg.)');
set(gca, 'fontsize', 20);
% Plot
subplot(1, 2, 2);
aux.otmatPlotNEW(obdata_mitral);

% FIGURE 2C: Same with different populations
figure2(9) = figure;
set(figure2(9), 'name', 'OT Cond. Prob.');

% Image
subplot(2, 2, 1);
imagesc(obdata_lowipr.data.OBOT.conProb);
colormap(jet);
title('P_{OB \rightarrow OT}(OB \rightarrow region) (Mitral)');
xlabel('OT slice');
xlabel('OT slice');
ylabel('Region');
yticks(1:5);
yticklabels(obdata_lowipr.prjRegName([1, 2, 3, 5, 6]));
B = colorbar;
title(B, 'P_{OB \rightarrow OT}(OB \rightarrow reg.) narrow');
set(gca, 'fontsize', 20);
% Plot
subplot(2, 2, 3);
aux.otmatPlotNEW(obdata_lowipr);

% Image
subplot(2, 2, 2);
imagesc(obdata_highipr.data.OBOT.conProb);
colormap(jet);
title('P_{OB \rightarrow OT}(OB \rightarrow region) (Mitral)');
xlabel('OT slice');
xlabel('OT slice');
ylabel('Region');
yticks(1:5);
yticklabels(obdata_highipr.prjRegName([1, 2, 3, 5, 6]));
B = colorbar;
title(B, 'P_{OB \rightarrow OT}(OB \rightarrow reg.) broad');
set(gca, 'fontsize', 20);

% Plot
subplot(2, 2, 4);
aux.otmatPlotNEW(obdata_highipr);
savefig(figure2(8:9), 'data/figures/obInjection_OT.fig');

%% FIGURE ???: Histograms of different raw barcode counts
figure2(10) = figure;
clf(figure2(10));
set(figure2(10), 'name', 'Barcode counts');
sto_counts_all = sum(vertcat(data_ob.prjImg), 1);
sto_avgCounts = cell2mat(cellfun(@(x) mean(x, 1, 'omitnan'), ...
  {data_ob.prjImg}, 'UniformOutput', 0)');
% Aggregate plot
subplot(2, 1, 1);
plot(sto_counts_all);
xlabel('Slices');
ylabel('Total Count');
title('Total projection strength per slice');
xticks(cellfun(@mean, data_ob(1).prjRegInd));
xticklabels(data_ob(1).prjRegName);
xlim([1, data_ob(1).nPrjSli]);
this_y = ylim;
for r = 1:(data_ob(1).nPrjReg - 1)
  line([1, 1] * (data_ob(1).prjRegInd{r}(end) + .5), this_y, 'Color', 'k');
end
ylim(this_y);
% Different brains
subplot(2, 1, 2);
plot(sto_avgCounts');
legend({data_ob.name}, 'Location', 'southoutside', ...
  'NumColumns', length(data_ob), 'AutoUpdate', false);
xlabel('Slices');
ylabel('Average (per barcode) Count');
title('Total projection strength per slice');
xticks(cellfun(@mean, data_ob(1).prjRegInd));
xticklabels(data_ob(1).prjRegName);
xlim([1, data_ob(1).nPrjSli]);
this_y = ylim;
for r = 1:(data_ob(1).nPrjReg - 1)
  line([1, 1] * (data_ob(1).prjRegInd{r}(end) + .5), this_y, 'Color', 'k');
end
ylim(this_y);

