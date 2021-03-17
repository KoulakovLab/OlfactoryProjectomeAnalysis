function otmatPlotNEW(set, CONFINT)
%OBMATPLOT Plot olfactory bulb outputs to piriform and other region

if ~exist('CONFINT', 'var')
  CONFINT = .95;
end
STDAMNT = norminv(.5 * (1 + CONFINT));

% Color choices
sto_col = lines(5);

% Refer to each plot specifically
sto_plots = gobjects(5, 1);
sN = size(set.data.OBOT.conProb, 2);

% Some storage values
xval = (1:sN)';
sto_fitX = [1; sN];
sto_regName = {'AON', 'APC', 'PPC', 'CoA', 'lENT'};
for i = 1:5
  yval = set.data.OBOT.conProb(i, xval)';
  % Draw the points as a line
  sto_plots(i) = plot(xval, yval, 'LineWidth', 2.5, ...
    'Color', sto_col(i, :));
  hold('on');
  % Draw the regional fits
  plot(sto_fitX, set.data.OBOT.cpFit_b_linear(:, i), ...
    '--', 'Color', sto_col(i, :), 'LineWidth', 1.5);
  fill([sto_fitX; sto_fitX(end:-1:1)], [ ...
    set.data.OBOT.cpFit_b_linear(:, i); flipud(...
    set.data.OBOT.cpFit_b_linear(:, i))] + STDAMNT .* [ ...
    set.data.OBOT.cpFit_b_linear_std(:, i); flipud(-...
    set.data.OBOT.cpFit_b_linear_std(:, i))], ...
    sto_col(i, :), 'FaceAlpha', .3, 'EdgeAlpha', 0);
end
% Draw APC-PPC line
hold('off');
sto_yax = ylim;
ylim([max(0, sto_yax(1)), sto_yax(2)]);
xlim([1, sN]);
grid('on');
xlabel('Olfactory Tubercle slice');
ylabel('P(region|OT)');
sto_leg = legend(sto_plots, sto_regName);
title(sto_leg, 'Regions');
title('Conditional probability: P(OB \rightarrow region|OB \rightarrow OT)');
end
