function obmatPlotNEW(set, CONFINT)
%OBMATPLOT Plot olfactory bulb outputs to piriform and other region

if ~exist('CONFINT', 'var')
  CONFINT = .95;
end
STDAMNT = norminv(.5 * (1 + CONFINT));

% Color choices
sto_col = lines(4);

% Refer to each plot specifically
sto_plots = gobjects(4, 1);
sN = size(set.data.OBPC.conProb_pc, 2);

% Some storage values
sto_xval = [(1:sN), (sN:-1:1)];
sto_apc = 1:set.nPrjRegSli(2);
sto_ppc = set.nPrjRegSli(2) + (1:set.nPrjRegSli(3));
sto_regName = {'AON', 'OT', 'CoA', 'lENT'};
sto_fitX = [1; set.nPrjRegSli(2) + .5; sN];
for i = 1:4
  xval_apc = sto_apc';
  yval_apc = set.data.OBPC.conProb_pc(i, sto_apc)';
  xyval_apc = [xval_apc, yval_apc];
  xval_ppc = sto_ppc';
  yval_ppc = set.data.OBPC.conProb_pc(i, sto_ppc)';
  xyval_ppc = [xval_ppc, yval_ppc];
  xval = [xval_apc; xval_ppc];
  yval = [yval_apc; yval_ppc];
  xyval = [xval, yval];
  % Draw the points as a line
  sto_plots(i) = plot(xval, yval, 'LineWidth', 2.5, ...
    'Color', sto_col(i, :));
  hold('on');
  % Draw the regional fits
  plot(sto_fitX, set.data.OBPC.cpFit_b_spline(:, i), ...
    '--', 'Color', sto_col(i, :), 'LineWidth', 1.5);
  fill([sto_fitX; sto_fitX(end:-1:1)], [ ...
    set.data.OBPC.cpFit_b_spline(:, i); flipud(...
    set.data.OBPC.cpFit_b_spline(:, i))] + STDAMNT .* [ ...
    set.data.OBPC.cpFit_b_spline_std(:, i); flipud(-...
    set.data.OBPC.cpFit_b_spline_std(:, i))], ...
    sto_col(i, :), 'FaceAlpha', .3, 'EdgeAlpha', 0);
end
% Draw APC-PPC line
sto_yax = ylim;
plot((.5 + set.nPrjRegSli(2)) * ones(1, 2), [0, 1], 'LineWidth', 3, ...
  'LineStyle', '--', 'Color', 'k');
hold('off');
ylim(sto_yax);
xlim([1, sN]);
grid('on');
xlabel('Piriform Cortex slice');
ylabel('P(region|Piriform Cortex)');
sto_leg = legend(sto_plots, {'AON', 'OT', 'CoA', 'lENT'});
title(sto_leg, 'Regions');
title('Conditional probability: P(OB \rightarrow region|OB \rightarrow PC) (Mitral)');
end
