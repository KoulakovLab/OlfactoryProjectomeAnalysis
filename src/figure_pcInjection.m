% Figure 4 is me; PC output as slice based
%   a) Overview of the data
%   * The conditional probability in;
%       b) Matrix form
%       c) Plot form
%       d) Subspace form
%   * Inter piriform connections
%       e) Cond. prob Matrix & discrepency of the connections
%           * Both regional and slice based connectome
% Altfig 4 is the same; but with filtered barcodes (and without conProb)
% Auxillary figure is me seperating the dataset into square blocks
%   * Has a overview image
%   * Has both conditional probability and example discrepency images

sto_cmp = [...
    hsv2rgb([0.0 * ones(128, 1), linspace(1, 0, 128)', ones(128, 1)]); ...
    hsv2rgb([0.5 * ones(128, 1), linspace(0, 1, 128)', ones(128, 1)])];
figure4 = gobjects(3, 1);
SKIP = 2;

%--------------%
%---OVERVIEW---%
%--------------%

% FIGURE 4A: Overview
figure4(1) = figure;
set(figure4(1), 'name', 'Overview');
% Put these plots separately side by side
subplot(1, 2, 1);
pcdata.plotImg('src');
colormap('parula');
colorbar;
% Find the ant-pos seperation barcode
sto_apBrc = find(pcdata.data.brcVis > pcdata.nSrcRegSli(1), 1);
hold('on');
plot(xlim, (sto_apBrc + .5) * ones(1, 2), ...
    ':', 'LineWidth', 3, 'Color', .5 * ones(1, 3));
hold('off');
title('Injection Site');
xlabel('Piriform Cortex (Slices)');
set(gca, 'XTickMode', 'auto', 'XTickLabelMode', 'auto')
set(gca, 'fontsize', 15);
subplot(1, 2, 2);
pcdata.plotImg('prj');
hold('on');
plot(xlim, (sto_apBrc + .5) * ones(1, 2), ...
    ':', 'LineWidth', 3, 'Color', .5 * ones(1, 3));
hold('off');
colorbar;
title('Projection sites');
xlabel('Target region');
ylabel('');
yticklabels([]);
set(gca, 'fontsize', 13);

savefig(figure4(1), 'data/figures/pcInjection_overview.fig');

%-----------------------------%
%---CONDITIONAL PROBABILITY---%
%-----------------------------%

figure4(2) = figure;
set(figure4(2), 'name', 'Cond. Prob.');

% FIGURE 4B: Conprob matrix
subplot(2, 2, 1);
imagesc(pcdata.data.PCout.conProb_pc);
colormap(jet);
title('P_{PC slice}(PC \rightarrow region)');
xlabel('PC slice');
% Draw A-P line
nSli = pcdata.nSrcRegSli(1:2);
line((nSli(1) + .5) * ones(1, 2), ylim, 'LineWidth', 3, ...
    'Color', .5 * ones(1, 3));
xticks(cumsum([0, nSli(1:(end-1))]) + .5 * nSli);
xticklabels({'APC', 'PPC'});
ylabel('Region');
yticks(1:4);
regions = {'AON', 'OT', 'CoA', 'ENT'};
yticklabels(regions);
B = colorbar;
title(B, 'P_{PC}(PC \rightarrow reg.)');
set(gca, 'fontsize', 14);

% Conprob plot
subplot(2, 2, 2);
aux.pcmatPlotNEW(pcdata);

% Conprob PCA
sN = size(pcdata.data.PCPC.conProb, 2) - 2 * SKIP;
sto_apc = 1:(pcdata.nSrcRegSli(1) - SKIP);
sto_ppc = (pcdata.nSrcRegSli(1) - SKIP) + (1:(pcdata.nSrcRegSli(2) - SKIP));
for d = 1:2
    subplot(2, 2, 2 + d)
    if d == 1
        scatter3(...
            pcdata.data.PCout.conProb_sliPca3(:, 1), ...
            pcdata.data.PCout.conProb_sliPca3(:, 2), ...
            pcdata.data.PCout.conProb_sliPca3(:, 3), ...
            100, 1:sN, 'filled');
        % Do fits to APC & PPC
        hold('on');
        [v, x] = aux.fitLineMultiN(...
            pcdata.data.PCout.conProb_sliPca3(sto_apc, :));
        plot3(...
            [x(1)-v(1), x(1)+v(1)], ...
            [x(2)-v(2), x(2)+v(2)], ...
            [x(3)-v(3), x(3)+v(3)], ...
            'LineWidth', 3, 'LineStyle', ':', 'Color', 'k');
        [v, x] = aux.fitLineMultiN(...
            pcdata.data.PCout.conProb_sliPca3(sto_ppc, :));
        plot3(...
            [x(1)-v(1), x(1)+v(1)], ...
            [x(2)-v(2), x(2)+v(2)], ...
            [x(3)-v(3), x(3)+v(3)], ...
            'LineWidth', 3, 'LineStyle', ':', 'Color', 'k');
        zlabel('PCA_{3}');
    elseif d == 2
        scatter(...
            pcdata.data.PCout.conProb_sliPca2(:, 1), ...
            pcdata.data.PCout.conProb_sliPca2(:, 2), ...
            100, 1:sN, 'filled');
        % Do fits to APC & PPC
        hold('on');
        [v, x] = aux.fitLineMultiN(...
            pcdata.data.PCout.conProb_sliPca2(sto_apc, :));
        plot(...
            [x(1)-v(1), x(1)+v(1)], ...
            [x(2)-v(2), x(2)+v(2)], ...
            'LineWidth', 3, 'LineStyle', ':', 'Color', 'k');
        [v, x] = aux.fitLineMultiN(...
            pcdata.data.PCout.conProb_sliPca2(sto_ppc, :));
        plot(...
            [x(1)-v(1), x(1)+v(1)], ...
            [x(2)-v(2), x(2)+v(2)], ...
            'LineWidth', 3, 'LineStyle', ':', 'Color', 'k');
    end
    hold('off');
    xlabel('PCA_{1}');
    ylabel('PCA_{2}');
    colormap('cool');
    axis('image');
    grid('on');
    sto_bar = colorbar;
    set(sto_bar, 'YDir', 'reverse');
    title(sto_bar, 'Slice (A \rightarrow P)');
    title('Piriform Slices in PCA space of cond. prob.');
    set(gca, 'fontsize', 20);
end

savefig(figure4(2), 'data/figures/pcInjection_condProb.fig');

% FIGURE 4C: Inter-PC connectome
figure4(3) = figure;
set(figure4(3), 'name', 'Inter-PC Connectome');
sgtitle('Inter-PC connectome');

% The conditional prob image
sto_plt = subplot(2, 2, 1);
colormap(sto_plt, 'jet');
sto_img = imagesc(log10(pcdata.data.PCPC.conProb));
set(gca, 'Color', .5 * ones(1, 3));
mask = double( ...
    (1 - tril(triu(ones(pcdata.nSrcSli), -PC_SOMA-1), PC_SOMA)) ...
    & (pcdata.data.PCPC.conProb_pval < .05));
sto_img.AlphaData = mask;
xlabel('Location of Soma');
ylabel('Projected slice');
title('Conditional probability; PC projections given soma loc.');
hold('on');
plot((.5 + pcdata.nSrcRegSli(1)) * ones(1, 2), ylim, ...
    ':', 'Color', .5 * ones(1, 3), 'LineWidth', 3);
plot(xlim, (.5 + pcdata.nSrcRegSli(1)) * ones(1, 2), ...
    ':', 'Color', .5 * ones(1, 3), 'LineWidth', 3);
hold('off');
set(gca, 'fontsize', 15);
b = colorbar;
title(b, 'log_{10} Prob.');
axis('image');
% AP discrepency
sto_plt = subplot(2, 2, 2);
sto_mat = pcdata.data.PCPC.conProb - pcdata.data.PCPC.conProb';
% Make into log_10 values
sto_img = log10(abs(sto_mat));
sto_min = floor(min(sto_img(~isinf(sto_img(:))), [], 'all'));
sto_max = ceil(max(sto_img(~isinf(sto_img(:))), [], 'all'));
sto_img(isinf(sto_img(:))) = sto_min;
% Make it so that + values range from [0,1], and neg range from [-1,0]
sto_img = sign(sto_mat) .* ((sto_img - sto_min) / (sto_max - sto_min));
% Create labels for the colorbar running appropriately
sto_cbl = [...
    sprintfc('-10^{%d}', (sto_max:-1:(sto_min+1))'); ...
    {'0'}; ...
    sprintfc('+10^{%d}', ((sto_min+1):sto_max)')];
sto_cbt = linspace(-1, 1, length(sto_cbl))';
% Create a mask that removes the top triangle; and obscures the low-p vals
sto_mask = tril(ones(pcdata.nSrcSli), -PC_SOMA-1);
% Draw img
imagesc(sto_img, 'AlphaData', sto_mask);
axis('image');
set(gca, 'Color', .1 * ones(1, 3));
colormap(sto_plt, sto_cmp);
xlabel('Anterior slice');
ylabel('Posterior slice');
title('AP axis Discrepency: (A \rightarrow P) - (P \rightarrow A)');
% Colorbar
sto_cbr = colorbar;
title(sto_cbr, '\Delta log_{10}prob.');
sto_cbr.Ticks = sto_cbt;
sto_cbr.TickLabels = sto_cbl;
% Seperate APC and PPC
hold('on');
plot((.5 + pcdata.nSrcRegSli(1)) * ones(1, 2), ylim, ...
    ':', 'Color', .5 * ones(1, 3), 'LineWidth', 3);
plot(xlim, (.5 + pcdata.nSrcRegSli(1)) * ones(1, 2), ...
    ':', 'Color', .5 * ones(1, 3), 'LineWidth', 3);
hold('off');
set(gca, 'fontsize', 15);
% Regional connectome
sto_axs = subplot(2, 2, 3);
imagesc(pcdata.data.PCPC.conProbReg);
colormap(sto_axs, bone);
xlabel('From');
xticks(1:2);
xticklabels({'APC', 'PPC'});
ylabel('To');
yticks(1:2);
yticklabels({'APC', 'PPC'});
b = colorbar;
title(b, 'Prob.');
sgtitle('Regional connectome in Piriform Cortex: Cond. Prob.');
% Regional connectome as a directed graph
subplot(2, 2, 4);
dg = digraph((pcdata.data.PCPC.conProbReg'), ...
    {'APC', 'PPC'});
plot(dg, '-o', ...
    'LineWidth', 5, 'ArrowSize', 25, ...
    'EdgeLabel', dg.Edges.Weight);
  
savefig(figure4(3), 'data/figures/pcInjection_innerPC.fig');
