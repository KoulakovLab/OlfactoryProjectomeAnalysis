% Figure 5 is me; PC output together wit OB output
STEP = 2;
SKIP_OB = 1;
SKIP_PC = 2;
CON_METHOD = 'bezier';
CON_NUMPTS = 100;

% CON_CMAP = @(x) aux.pmkmp(x, 'IsoL');
CON_CMAP = @(x) cool(x);
figure5 = gobjects(1, 1);
figure5(1) = figure;
set(figure5(1), 'name', 'Evolutions of Cond. Prob.');

sN_ob = size(obdata_mitral.data.OBPC.conProb_pc, 2);
sN_pc = size(pcdata.data.PCout.conProb_pc, 2);

sI_ob = (1 + SKIP_OB):(sN_ob - SKIP_OB);
sI_pc = (1 + SKIP_PC):(sN_pc - SKIP_PC);

if length(sI_ob) ~= length(sI_pc)
  sN = 2 * max(length(sI_ob), length(sI_pc));
  ob_vals = interp1(...
    linspace(0, 1, length(sI_ob)), ...
    obdata_mitral.data.OBPC.conProb_pc(:, sI_ob)', ...
    linspace(0, 1, sN));
  pc_vals = interp1(...
    linspace(0, 1, length(sI_pc)), ...
    pcdata.data.PCout.conProb_pc(:, sI_pc)', ...
    linspace(0, 1, sN));
else
  sN = length(sI_ob);
  ob_vals = obdata_mitral.data.OBPC.conProb_pc(:, sI_ob)';
  pc_vals = pcdata.data.PCout.conProb_pc(:, sI_pc)';
end

lc_vals = linspace(0, 1, sN);
sto_names = obdata_mitral.prjRegName([1, 4, 5, 6]);

for i = 1:4
  sto_plt = subplot(2, 2, i);
  scatter( ...
    ob_vals(1:end, i), ...
    pc_vals(1:end, i), 100, ...
    lc_vals(1:end), 'filled', 'MarkerEdgeColor', 'k');
  [sto_cor, ~, sto_cpv] = aux.spear( ...
    ob_vals(1:end, i), pc_vals(1:end, i))
  colormap(CON_CMAP(256));
  caxis([0 1]);
  if i == 1
    b = colorbar('Direction', 'reverse');
  else
    b = colorbar;
  end
  title(b, 'A \rightarrow P position');
  % Draw a connecting line
  hold('on');
  sto_con = aux.bezier([ ...
    ob_vals(1:end, i), ...
    pc_vals(1:end, i)], CON_NUMPTS);
  cl_vals = CON_CMAP(CON_NUMPTS);
  for j = 2:CON_NUMPTS
    plot(sto_con([j - 1; j], 1), sto_con([j - 1; j], 2), ...
      '-', 'LineWidth', 2, 'Color', cl_vals(j, :));
  end
  hold('off');
  % IF it's the first graph; reverse AP axis direction
  xlabel(['P(OB \rightarrow ', sto_names{i}, ' | OB \rightarrow PC)']);
  ylabel(['P(PC \rightarrow ', sto_names{i}, ' | PC)']);
  title(['Projection comparison: ', sto_names{i}]);
  set(gca, 'fontsize', 15);
  grid('on');

  % Annotate the p and rho value
  annotation('textbox', aux.pltbox(sto_plt) .* [1 1 0.5 0.5], 'String', ...
    {['\rho = ', num2str(sto_cor, 2)], ['p = ', num2str(sto_cpv, 2)]}, ...
    'FitBoxToText','on');
end

savefig(figure5, 'data/figures/mapping_overview.fig');