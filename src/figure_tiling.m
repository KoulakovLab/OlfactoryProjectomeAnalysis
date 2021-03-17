% Figure 3 is me; OB Tiling
%   Fig 3a)
%   * It's the tiling with the exponential fits
%   * Will do the following datasets
%      (1): All brains
%      (1): Only the first 5 brains
%      (6): All brains seperately
%   * Will have the following layouts
%      (1): All mitral, low-ipr, high-ipr
%      (1): All mitral, mitral-1, mitral-2
%      (1): All, mitral, tufted
%   Rest of the figures;
%   B) AON vs APC ordering figure
%       * Outside of their regions; they are either random, or some other
%       sorting-
%       * APC vs PPC sorting
%   Auxfigure 3a)
%   * Trend of APC exponential wrt different IPR thresholds
%   * Display each dataset seperately.


%% FIGURE 3A: Overview
ds_full = [ ...
  obdata, ...
  obdata_abone, ...
  data_ob'];
ds_mit0 = [ ...
  obdata_mitral, ...
  obdata_abone_mitral, ...
  data_ob_mitral'];
ds_tuft = [ ...
  obdata_tufted, ...
  obdata_abone_tufted, ...
  data_ob_tufted'];
ds_mitL = [ ...
  obdata_lowipr, ...
  obdata_abone_lowipr, ...
  data_ob_lowipr'];
ds_mitH = [ ...
  obdata_highipr, ...
  obdata_abone_highipr, ...
  data_ob_highipr'];
DS = [ds_full; ds_mit0; ds_tuft; ds_mitL; ds_mitH];
DS_CONF = [ ...
  2, 4, 5;    % Mitral; IPR
  1, 2, 3];   % Mitral vs Tufted
DS_CONF_NAME = {'IPR', 'Cell type'};

fig3a = gobjects(size(DS_CONF, 1), size(DS, 2));

for s = 1:size(DS, 2)               % For each dataset type
  for cc = 1:size(DS_CONF, 1)      % For each configuration
    fig3a(cc, s) = figure;
    set(fig3a(cc, s), 'name', ...
      ['Tiling: ', DS_CONF_NAME{cc}, ' ', DS(1, s).name]);
    sgtitle(['Tiling in ', DS(1, s).name, ', seperated by ', ...
      DS_CONF_NAME{cc}]);
    set(gca, 'FontSize', 20);
    for i = 1:size(DS_CONF, 2)
      % As a plot
      subplot(2, size(DS_CONF, 2), i);
      c = DS_CONF(cc, i);
      % Draw image
      DS(c, s).plotImg('prj', true, true);
      hold('on')
      % Seperate region with lines
      [~, max_id] = max(DS(c, s).prjImg, [] ,2);
      max_reg = sum(max_id > cumsum(DS(c, s).nPrjRegSli), 2) + 1;
      reg_brc = sum(max_reg == (1:DS(c, s).nPrjReg), 1);
      for r = 2:DS(c, s).nPrjReg
        plot(xlim, (.5 + sum(reg_brc(1:(r-1)))) * ones(1, 2), ...
          ':', 'Color', .5 * ones(1, 3), 'LineWidth', 2);
      end
      % Draw all the tiling fits
      for r = 1:4
        if isempty(DS(c, s).data.expFit{r})
          continue;
        end
        plot( ...
          DS(c, s).data.expFit{r}(:, 1), ...
          DS(c, s).data.expFit{r}(:, 2), ...
          'Color', [.9 * ones(1, 3), .75], 'LineWidth', 4);
      end
      hold('off');
      
      % As a histogram
      subplot(2, size(DS_CONF, 2), size(DS_CONF, 2) + i);
      % Linearize the indices
      [sto_ids, ~, sto_cls] = unique(DS(c, s).brcId(:, 1));
      sto_name = DS(c, s).brcName{1}(sto_ids);
      % Build the bar plot; each column will be one brain
      sto_cnt = zeros(DS(c, s).nPrjSli, length(sto_ids));
      for b = 1:length(sto_ids)
        sto_vis = ...
          DS(c, s).brcPrjVis(DS(c, s).brcId(:, 1) == sto_ids(b));
        sto_cnt(:, b) = sum(sto_vis == (1:DS(c, s).nPrjSli), 1)';
      end
      bar(sto_cnt, 1, 'stacked');
      % Draw seperators
      hold('on');
      for r = 2:DS(c, s).nPrjReg
        plot(ones(1, 2) * (.5 + ...
          max(DS(c, s).prjRegInd{r - 1}, [], 'all')), ...
          ylim, '-k', 'LineWidth', 1);
      end
      hold('off');
      % Labels
      legend(sto_name);
      ylabel('# of barcodes');
      xlabel('Slices');
      xticks(cellfun(@mean, DS(c, s).prjRegInd));
      xticklabels(DS(c, s).prjRegName);
    end
  end
end

% Save all 3A Figure plots
savefig(fig3a(:), 'data/figures/tiling_overview.fig');


%% FIGURE 3B-C: APC vs AON vs PPC orderings
SORT_METHOD = 'strength';
NOISE = 1/5 * randn(obdata_mitral.nBrc, 3);
figure3 = gobjects(3, 1);

% Create images with specific region orderings
sto_sli = obdata_mitral.prjRegInd(1:3);
[aon_str, aon_id] = max( ...
  obdata_mitral.prjImg(:, obdata_mitral.prjRegInd{1}), [], 2);
[apc_str, apc_id] = max( ...
  obdata_mitral.prjImg(:, obdata_mitral.prjRegInd{2}), [], 2);
[ppc_str, ppc_id] = max( ...
  obdata_mitral.prjImg(:, obdata_mitral.prjRegInd{3}), [], 2);
switch SORT_METHOD
  case 'random'
    [~, aon_ind] = sortrows([aon_id, rand(obdata_mitral.nBrc, 1)]);
    [~, apc_ind] = sortrows([apc_id, rand(obdata_mitral.nBrc, 1)]);
    [~, ppc_ind] = sortrows([ppc_id, rand(obdata_mitral.nBrc, 1)]);
  case 'strength'
    [~, aon_ind] = sortrows([aon_id, aon_str, rand(obdata_mitral.nBrc, 1)]);
    [~, apc_ind] = sortrows([apc_id, apc_str, rand(obdata_mitral.nBrc, 1)]);
    [~, ppc_ind] = sortrows([ppc_id, ppc_str, rand(obdata_mitral.nBrc, 1)]);
  otherwise
    [~, aon_ind] = sortrows(aon_id);
    [~, apc_ind] = sortrows(apc_id);
    [~, ppc_ind] = sortrows(ppc_id);
end
aon_img = obdata_mitral.prjImg(aon_ind, [sto_sli{:}]);
aon_img = aon_img ./ max(aon_img, [], 2);
apc_img = obdata_mitral.prjImg(apc_ind, [sto_sli{:}]);
apc_img = apc_img ./ max(apc_img, [], 2);
ppc_img = obdata_mitral.prjImg(ppc_ind, [sto_sli{:}]);
ppc_img = ppc_img ./ max(ppc_img, [], 2);

% FIGURE 3B: AON vs APC ordering
figure3(1) = figure;
set(figure3(1), 'name', 'Ordering, AON vs APC');
% Sorting image
subplot(2, 2, 1:2:3);
imagesc([aon_img; apc_img]);
colormap('parula');
% Draw y-axis label and lines
hold('on');
plot(xlim, (obdata_mitral.nBrc + .5) * ones(1, 2), ...
  ':', 'Color', .5 * ones(1, 3), 'LineWidth', 3);
yticks([0, obdata_mitral.nBrc] + ((obdata_mitral.nBrc + 1) / 2));
yticklabels({'AON sorted', 'APC sorted'});
ylabel('Barcodes (sorted by region maxima)');
% Draw x-axis label and lines
plot(cellfun( ...
  @(x, y) mean([max(x, [], 'all'), min(y, [], 'all')]), ...
  sto_sli(1:(end-1)), sto_sli(2:end)) .* ones(2, 1), ...
  [ylim; ylim]', ':', 'Color', .5 * ones(1, 3), 'LineWidth', 3);
hold('off');
xticks(cellfun(@mean, sto_sli));
xticklabels({'AON', 'APC', 'PPC'});
xlabel('Slices');
% Correlation and p-value
title(sprintf(...
  'Corr: %.2f P-Val: %.2f', ...
  obdata_mitral.data.regOrdCorr(1, 2), ...
  obdata_mitral.data.regOrdPval(1, 2)));
% Do in this weird space
subplot(2, 2, 2);
sto_img = zeros(obdata_mitral.nPrjRegSli(2), obdata_mitral.nPrjRegSli(1));
sto_ind = sub2ind(obdata_mitral.nPrjRegSli([2, 1]), ...
  obdata_mitral.data.brcOrd(:, 2), obdata_mitral.data.brcOrd(:, 1));
for b = 1:obdata_mitral.nBrc
  sto_img(sto_ind(b)) = sto_img(sto_ind(b)) + 1;
end
imagesc(sto_img);
xlabel('AON ordering');
ylabel('APC ordering');
title('Ordering counts');
colorbar;
subplot(2, 2, 4);
scatter( ...
  obdata_mitral.data.brcOrd(:, 1) + NOISE(:, 1), ...
  obdata_mitral.data.brcOrd(:, 2) + NOISE(:, 2), ...
  3, 'filled');
xticks(1:obdata_mitral.nPrjRegSli(1));
yticks(1:obdata_mitral.nPrjRegSli(2));
xlabel('AON ordering');
ylabel('APC ordering');
grid('on');
title('Ordering space');
sgtitle('Ordering; AON vs APC');
set(gca, 'FontSize', 20);

% FIGURE 3C: APC vs PPC ordering
figure3(2) = figure;
set(figure3(2), 'name', 'Ordering, APC vs PPC');
% Sorting image
subplot(2, 2, 1:2:3);
imagesc([apc_img; ppc_img]);
colormap('parula');
% Draw y-axis label and lines
hold('on');
plot(xlim, (obdata_mitral.nBrc + .5) * ones(1, 2), ...
  ':', 'Color', .5 * ones(1, 3), 'LineWidth', 3);
yticks([0, obdata_mitral.nBrc] + ((obdata_mitral.nBrc + 1) / 2));
yticklabels({'APC sorted', 'PPC sorted'});
ylabel('Barcodes (sorted by region maxima)');
% Draw x-axis label and lines
plot(cellfun( ...
  @(x, y) mean([max(x, [], 'all'), min(y, [], 'all')]), ...
  sto_sli(1:(end-1)), sto_sli(2:end)) .* ones(2, 1), ...
  [ylim; ylim]', ':', 'Color', .5 * ones(1, 3), 'LineWidth', 3);
hold('off');
xticks(cellfun(@mean, sto_sli));
xticklabels({'AON', 'APC', 'PPC'});
xlabel('Slices');
% Correlation and p-value
title(sprintf(...
  'Corr: %.2f P-Val: %.2f', ...
  obdata_mitral.data.regOrdCorr(2, 3), ...
  obdata_mitral.data.regOrdPval(2, 3)));
% Do in the weird space
subplot(2, 2, 2);
sto_img = zeros(obdata_mitral.nPrjRegSli(3), obdata_mitral.nPrjRegSli(2));
sto_ind = sub2ind(obdata_mitral.nPrjRegSli([3, 2]), ...
  obdata_mitral.data.brcOrd(:, 3), obdata_mitral.data.brcOrd(:, 2));
for b = 1:obdata_mitral.nBrc
  sto_img(sto_ind(b)) = sto_img(sto_ind(b)) + 1;
end
imagesc(sto_img);
xlabel('APC ordering');
ylabel('PPC ordering');
title('Ordering counts');
colorbar;
subplot(2, 2, 4);
scatter( ...
  obdata_mitral.data.brcOrd(:, 2) + NOISE(:, 1), ...
  obdata_mitral.data.brcOrd(:, 3) + NOISE(:, 2), ...
  3, 'filled');
xticks(1:obdata_mitral.nPrjRegSli(2));
yticks(1:obdata_mitral.nPrjRegSli(3));
xlabel('APC ordering');
ylabel('PPC ordering');
grid('on');
title('Ordering space');
sgtitle('Ordering; APC vs PPC');
set(gca, 'FontSize', 20);

% FIGURE 3D: Multi-ordering ordering
figure3(3) = figure;
set(figure3(3), 'name', 'Ordering, AON, APC & PPC');
% Sorting image
subplot(2, 2, 1:2:3);
imagesc([aon_img; apc_img; ppc_img]);
colormap('parula');
% Draw y-axis label and lines
hold('on');
plot([xlim; xlim]', ...
  ((obdata_mitral.nBrc * (1:2)) + .5) .* ones(2, 1), ...
  ':', 'Color', .5 * ones(1, 3), 'LineWidth', 3);
yticks(((0:2) * obdata_mitral.nBrc) + ((obdata_mitral.nBrc + 1) / 2));
yticklabels({'AON sorted', 'APC sorted', 'PPC sorted'});
ylabel('Barcodes (sorted by region maxima)');
% Draw x-axis label and lines
plot(cellfun( ...
  @(x, y) mean([max(x, [], 'all'), min(y, [], 'all')]), ...
  sto_sli(1:(end-1)), sto_sli(2:end)) .* ones(2, 1), ...
  [ylim; ylim]', ':', 'Color', .5 * ones(1, 3), 'LineWidth', 3);
hold('off');
xticks(cellfun(@mean, sto_sli));
xticklabels({'AON', 'APC', 'PPC'});
xlabel('Slices');
% Do in the weird space
subplot(2, 2, 2);
scatter3( ...
  obdata_mitral.data.brcOrd(:, 1) + NOISE(:, 1), ...
  obdata_mitral.data.brcOrd(:, 2) + NOISE(:, 2), ...
  obdata_mitral.data.brcOrd(:, 3) + NOISE(:, 3), ...
  3, 'filled');
xticks(1:obdata_mitral.nPrjRegSli(1));
yticks(1:obdata_mitral.nPrjRegSli(2));
zticks(1:obdata_mitral.nPrjRegSli(2));
xlabel('AON ordering');
ylabel('APC ordering');
zlabel('PPC ordering');
grid('on');
title('Ordering space');
% Show correlation plot
sto_plt = subplot(2, 2, 4);
sto_han = imagesc(obdata_mitral.data.regOrdCorr, ...
  [-1, 1] * (max(abs(obdata_mitral.data.regOrdCorr ...
  - eye(size(obdata_mitral.data.regOrdCorr))), [], 'all')));
set(gca, 'Color', zeros(1, 3));
sto_han.AlphaData = max(0, 1 - 2 * obdata_mitral.data.regOrdPval);
colormap(sto_plt, aux.red2cyan);
colorbar;
xticks(1:obdata_mitral.nPrjReg);
yticks(1:obdata_mitral.nPrjReg);
xticklabels(obdata_mitral.prjRegName);
yticklabels(obdata_mitral.prjRegName);
title('Correlation of ordering');
% Image title
sgtitle('Ordering; AON vs APC vs PPC');
set(gca, 'FontSize', 20);

% Save figures
savefig(figure3(:), 'data/figures/tiling_ordering.fig');

%% AUXILLARY FIGURE 3A: IPR trend with threshold
sto_ana = [obdata_mitral,  obdata_abone_mitral, data_ob_mitral'];
auxfig3 = gobjects(length(sto_ana), 2);
sto_col = lines(2);
sto_lcl = sto_col(2, :);
sto_hcl = sto_col(1, :);

for i = 1:length(sto_ana)
  sto_smp = size(sto_ana(i).data.exp_imgXval, 1);
  auxfig3(i, 1) = figure;
  set(auxfig3(i, 1), 'name', 'Bias trend vs. IPR threshold');
  set(gca, 'FontSize', 20);
  % Draw the abundance plot
  sto_axs = subplot(1, 2, 1);
  yyaxis('left');
  sto_plt = area( ...
    sto_ana(i).data.exp_iprVal, ...
    [sto_ana(i).data.exp_nLow, sto_ana(i).data.exp_nHigh], ...
    'FaceAlpha', .75);
  sto_plt(1).FaceColor = sto_lcl;
  sto_plt(2).FaceColor = sto_hcl;
  ylabel('Barcode # in partition');
  ylim([0, sto_ana(i).data.exp_nLow(1) + sto_ana(i).data.exp_nHigh(1)]);
  % Draw the gamma trend for both partitions
  yyaxis('right');
  plot(sto_ana(i).data.exp_iprVal, ...
    sto_ana(i).data.exp_lowGamma, ...
    '-v', 'Color', sto_lcl, 'LineWidth', 3, 'MarkerSize', 10, ...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
  hold('on');
  plot(sto_ana(i).data.exp_iprVal, ...
    sto_ana(i).data.exp_highGamma, ...
    '-^', 'Color', sto_hcl, 'LineWidth', 3, 'MarkerSize', 10, ...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
  ylabel('Exponential \gamma of fit');
  sto_vals = [sto_ana(i).data.exp_highGamma; sto_ana(i).data.exp_lowGamma];
  ylim([min(sto_vals), max(sto_vals)]);
  % Draw the full value as a line
  plot( ...
    [min(sto_ana(i).data.exp_iprVal), max(sto_ana(i).data.exp_iprVal)], ...
    ones(1, 2) * sto_ana(i).data.expVal(2), ...
    ':k', 'LineWidth', 3);
  hold('off');
  % Label the plot appropriately
  legend({ ...
    '# of low-IPR barcodes', '# of high-IPR barcodes', ...
    'Low-ipr \gamma', 'High-IPR \gamma', 'Whole set \gamma'}, ...
    'Location', 'southoutside');
  title('Partition size & \gamma evolution');
  % Modify plot a bit
  ax = gca;
  ax.FontSize = 15;
  ax.YAxis(1).Color = 'k';
  ax.YAxis(2).Color = 'k';
  xlim([min(sto_ana(i).data.exp_iprVal), max(sto_ana(i).data.exp_iprVal)]);
  % Draw the true images as well
  subplot(1, 2, 2);
  % Colorize the image to low and high ipr
  % Assign a master title
  sgtitle([sto_ana(i).name, '; IPR threshold vs Tiling (in APC)']);
  % Create a color coded image; orange is low while blue is high
  % High ipr on top
  sto_img = zeros([size(sto_ana(i).data.exp_imgBkg'), 3]);
  sto_sli = sto_ana(i).nPrjRegSli(2);
  for p = 1:size(sto_ana(i).data.exp_imgXval, 1)
    sto_sin = ((p - 1) * sto_sli) + (1:sto_sli);
    sto_lin = (1:sto_ana(i).data.exp_nLow(p)) + sto_ana(i).data.exp_nHigh(p);
    sto_hin = (1:sto_ana(i).data.exp_nHigh(p));
    for c = 1:3
      % The high-ipr barcodes in blue
      sto_img(sto_sin, sto_hin, c) = sto_hcl(c) * ...
        (sto_ana(i).data.exp_imgBkg(sto_hin, sto_sin)');
      % The low-ipr barcodes in orange
      sto_img(sto_sin, sto_lin, c) = sto_lcl(c) * ...
        (sto_ana(i).data.exp_imgBkg(sto_lin, sto_sin)');
    end
  end
  % Draw this colorized image
  image(sto_img);
  % Draw the fits on top
  hold('on');
  for p = 1:(2 * sto_smp)
    plot( ...
      sto_ana(i).data.exp_imgXval{p}, ...
      sto_ana(i).data.exp_imgYval{p}, ...
      'Color', [1.0 * ones(1, 3), .75], 'LineWidth', 2);
  end
  % Draw lines to seperate the samplings
  for p = 2:sto_smp
    plot(xlim, (sto_sli * (p - 1) + .5) * ones(1, 2), ...
      ':', 'Color', [.5 * ones(1, 3), .75], 'LineWidth', 1);
  end
  yticks(((1:sto_smp) - .5) * sto_sli);
  yticklabels(sprintfc('%.1f', sto_ana(i).data.exp_iprVal));
  ylabel('IPR Threshold');
  hold('off');
  xticks([]);
  title('APC partitions vs. Threshold');
  % Draw the correlation figures
  auxfig3(i, 2) = figure;
  set(auxfig3(i, 2), 'name', 'AUX.3B: IPR and proj. corr');
  set(gca, 'FontSize', 20);
  % IPR vs VIS
  subplot(1, 2, 1);
  sto_dat = [sto_ana(i).brcPrjIpr, sto_ana(i).brcPrjVis];
  sto_cov = cov(sto_dat);
  sto_avg = mean(sto_dat, 1);
  scatter(sto_dat(:, 1), sto_dat(:, 2), 50, 'filled');
  hold('on');
  aux.covEllipse(sto_avg, sto_cov);
  hold('off');
  xlabel('Inverse Participation Ratio (IPR)');
  ylabel('Brightest slice location');
  % IPR vs Dist
  subplot(1, 2, 2);
  sto_dat = [sto_ana(i).brcPrjIpr, sto_ana(i).brcPrjDist];
  sto_cov = cov(sto_dat);
  sto_avg = mean(sto_dat, 1);
  scatter(sto_dat(:, 1), sto_dat(:, 2), 50, 'filled');
  hold('on');
  aux.covEllipse(sto_avg, sto_cov);
  hold('off');
  xlabel('Inverse Participation Ratio (IPR)');
  ylabel('Average projection distance');
  sgtitle([sto_ana(i).name, ', IPR vs Spatial trends']);
end

savefig(auxfig3(:), 'data/figures/tiling_IPRbiases.fig');

