function nnmatPlotNEW(SET, IND, CONFINT)
%OBMATPLOT Plot olfactory bulb outputs to piriform and other region

if ~exist('CONFINT', 'var')
  CONFINT = .95;
end
STDAMNT = norminv(.5 * (1 + CONFINT));

% Color choices
sto_col = lines(4);

% Refer to each plot specifically
sto_plots = gobjects(12, 1);
sN = size(SET.data.OBPC.conProb_pc, 2);

% Some storage values
sto_xval = [(1:sN), (sN:-1:1)];
sto_apc = 1:SET.nPrjRegSli(2);
sto_ppc = SET.nPrjRegSli(2) + (1:SET.nPrjRegSli(3));
sto_regName = {'AON', 'OT', 'CoA', 'lENT'};
sto_fitX = [1; SET.nPrjRegSli(2) + .5; sN];
for i = 1:4
  xval_apc = sto_apc';
  yval_apc = SET.data.typenet(IND).conProb(i, sto_apc)';
  yerr_apc = SET.data.typenet(IND).conProb_std(i, sto_apc)';
  xval_ppc = sto_ppc';
  yval_ppc = SET.data.typenet(IND).conProb(i, sto_ppc)';
  yerr_ppc = SET.data.typenet(IND).conProb_std(i, sto_ppc)';
  xval = [xval_apc; xval_ppc];
  yval = [yval_apc; yval_ppc];
  yerr = [yerr_apc; yerr_ppc];
  % Draw the points
  sto_plots(i) = errorbar(xval, yval, yerr, 's', 'LineWidth', 2.5, ...
    'MarkerSize', 5, 'Color', sto_col(i, :));
  hold('on');
  % Draw the regional fits
  sto_plots(i + 4) = plot(sto_fitX, SET.data.typenet(IND).cpFit_b_spline(:, i), ...
    '-', 'Color', sto_col(i, :), 'LineWidth', 1.5);
  fill([sto_fitX; sto_fitX(end:-1:1)], [ ...
    SET.data.typenet(IND).cpFit_b_spline(:, i); flipud(...
    SET.data.typenet(IND).cpFit_b_spline(:, i))] + STDAMNT .* [ ...
    SET.data.typenet(IND).cpFit_b_spline_std(:, i); flipud(-...
    SET.data.typenet(IND).cpFit_b_spline_std(:, i))], ...
    sto_col(i, :), 'FaceAlpha', .3, 'EdgeAlpha', 0);
  % Draw the original fits
  sto_plots(i + 8) = plot(sto_fitX, SET.data.OBPC.cpFit_b_spline(:, i), ...
    '--', 'Color', sto_col(i, :), 'LineWidth', 1.5);
end
% Draw APC-PPC line
sto_yax = ylim;
plot((.5 + SET.nPrjRegSli(2)) * ones(1, 2), sto_yax, 'LineWidth', 3, ...
  'LineStyle', '--', 'Color', 'k');
hold('off');
ylim(sto_yax);
xlim([1, sN]);
grid('on');
xlabel('PC slice');
ylabel('P(OB \rightarrow region|OB \rightarrow PC)');
legend(sto_plots, [sto_regName, ...
  cellfun(@(x) [x, ' fit'], sto_regName, 'UniformOutput', false), ...
  cellfun(@(x) [x, ' actual'], sto_regName, 'UniformOutput', false)]);
title(['Conditional probability: P(OB \rightarrow region|OB \rightarrow PC) (', ...
  SET.data.typenet(IND).name, ')']);

end
