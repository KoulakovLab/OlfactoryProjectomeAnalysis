function pcmatPlotNEW(SET, CONFINT)
%OBMATPLOT Plot olfactory bulb outputs to piriform and other region

if ~exist('CONFINT', 'var')
  CONFINT = .95;
end
STDAMNT = norminv(.5 * (1 + CONFINT));
SKIP = 2;

% Color choices
sto_col = lines(4);

% Refer to each plot specifically
sto_plots = gobjects(4, 1);
sN = size(SET.data.PCout.conProb_pc, 2);

% Some storage values
sto_xval = [((1 + SKIP):(sN - SKIP)), ((sN - SKIP):-1:(1 + SKIP))];
sto_apc = (1 + SKIP):(sN / 2);
sto_ppc = ((sN / 2) + 1):(sN - SKIP);
sto_regName = {'AON', 'OT', 'CoA', 'ENT'};
sto_fitX = [1 + SKIP; .5 * sN + .5; sN - SKIP];
for i = 1:4
  xval_apc = sto_apc';
  yval_apc = SET.data.PCout.conProb_pc(i, sto_apc)';
  xyval_apc = [xval_apc, yval_apc];
  xval_ppc = sto_ppc';
  yval_ppc = SET.data.PCout.conProb_pc(i, sto_ppc)';
  xyval_ppc = [xval_ppc, yval_ppc];
  xval = [xval_apc; xval_ppc];
  yval = [yval_apc; yval_ppc];
  xyval = [xval, yval];
  % Draw the points as a line
  sto_plots(i) = plot(xval, yval, 'LineWidth', 2.5, ...
    'Color', sto_col(i, :));
  hold('on');
  % Draw the regional fits
  plot(sto_fitX, SET.data.PCout.cpFit_b_spline(:, i), ...
    '--', 'Color', sto_col(i, :), 'LineWidth', 1.5);
  fill([sto_fitX; sto_fitX(end:-1:1)], [ ...
    SET.data.PCout.cpFit_b_spline(:, i); flipud(...
    SET.data.PCout.cpFit_b_spline(:, i))] + STDAMNT .* [ ...
    SET.data.PCout.cpFit_b_spline_std(:, i); flipud(-...
    SET.data.PCout.cpFit_b_spline_std(:, i))], ...
    sto_col(i, :), 'FaceAlpha', .3, 'EdgeAlpha', 0);
end
% Draw APC-PPC line
sto_yax = ylim;
plot((.5 + sN / 2) * ones(1, 2), [0, 1], 'LineWidth', 3, ...
  'LineStyle', '--', 'Color', 'k');
hold('off');
ylim(sto_yax);
xlim([1 + SKIP, sN - SKIP]);
grid('on');
xlabel('Piriform Cortex slice');
ylabel('P(region|Piriform Cortex)');
sto_leg = legend(sto_plots, {'AON', 'OT', 'CoA', 'lENT'});
title(sto_leg, 'Regions');
title('Conditional probability: P(PC \rightarrow region|PC)');
end
