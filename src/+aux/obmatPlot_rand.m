function obmatPlot_rand(set, CONFINT)
%OBMATPLOT Plot randomized olfactory bulb outputs to piriform and other region

if ~exist('CONFINT', 'var')
  CONFINT = .95;
end
STDAMNT = norminv(.5 * (1 + CONFINT));

% Color choices
sto_col = lines(4);

% Refer to each plot specifically
sto_plots = gobjects(4, 1);
sN = size(set.data.conProb, 2);

% Some storage values
sto_xval = [(1:sN), (sN:-1:1)];
sto_apc = 1:set.nPrjRegSli(2);
sto_ppc = set.nPrjRegSli(2) + (1:set.nPrjRegSli(3));
sto_textbox = cell((4 + 1), 2);
sto_textbox{1, 1} = 'Slope 1';
sto_textbox{1, 2} = 'Slope 2';
sto_regName = {'AON', 'OT', 'CoA', 'lENT'};
sto_fitX = [1; set.nPrjRegSli(2) + .5; sN];
for i = 1:4
  xval_apc = sto_apc';
  yval_apc = set.data.conProb_rand(i, sto_apc)';
  xyval_apc = [xval_apc, yval_apc];
  xval_ppc = sto_ppc';
  yval_ppc = set.data.conProb_rand(i, sto_ppc)';
  xyval_ppc = [xval_ppc, yval_ppc];
  xval = [xval_apc; xval_ppc];
  yval = [yval_apc; yval_ppc];
  xyval = [xval, yval];
  % Draw the points as a line
  sto_plots(i) = plot(xval, yval, 'LineWidth', 2.5, ...
    'Color', sto_col(i, :));
  hold('on');
  % Draw the regional fits
  plot(sto_fitX, set.data.conPFit_bAvg_rand(:, i), ...
    '--', 'Color', sto_col(i, :), 'LineWidth', 1.5);
  fill([sto_fitX; sto_fitX(end:-1:1)], [ ...
    set.data.conPFit_bAvg_rand(:, i); flipud(...
    set.data.conPFit_bAvg_rand(:, i))] + STDAMNT .* [ ...
    set.data.conPFit_bStd_rand(:, i); flipud(-...
    set.data.conPFit_bStd_rand(:, i))], ...
    sto_col(i, :), 'FaceAlpha', .3, 'EdgeAlpha', 0);
  % Fill the textbox annotation
  sto_textbox{1 + i, 1} = ['* m1=', ...
    num2str(set.data.conPFit_mAvg_rand(1, i), '%.3g'), ' (', ...
    sto_regName{i}, ')'];
  sto_textbox{1 + i, 2} = ['* m2=', ...
    num2str(set.data.conPFit_mAvg_rand(2, i), '%.3g'), ' (', ...
    sto_regName{i}, ')'];
end
% Draw APC-PPC line
sto_yax = ylim;
plot((.5 + set.nPrjRegSli(2)) * ones(1, 2), sto_yax, 'LineWidth', 3, ...
  'LineStyle', '--', 'Color', 'k');
hold('off');
ylim(sto_yax);
xlim([1, sN]);
grid('on');
xlabel('Piriform Cortex slice');
ylabel('P(region|Piriform Cortex)');
sto_leg = legend(sto_plots, {'AON', 'OT', 'CoA', 'lENT'});
title('Conditional probability: P(region|PC) (Mitral)');
title(sto_leg, 'Regions');
% set(gca, 'fontsize', 20);
% Write down the annotation
locn = get(gca, 'Position');
annotation('textbox', ...
  [max(0, locn(1) - .1), locn(2) + locn(4) / 4, min(locn(1), .1), locn(4) / 2], ...
  'String', sto_textbox(:), 'FitBoxToText', true);
end
