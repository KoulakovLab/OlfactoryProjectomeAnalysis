function plotSrcPrjProb(o, figid)
% PLOTSRCPRJPROB Displays the results of prob analysis
DIMRED_TYPE = 'pca';
STYLE = {'--', ':'};

% Open new figure if closed
if ~isfield(o.data, 'prob_figHandle')
    o.data.prob_figHandle = figure('name', ['Probability: ', o.name]);
elseif ~ishandle(o.data.prob_figHandle)
    o.data.prob_figHandle = figure('name', ['Probability: ', o.name]);
else
    % Don't steal focus
    set(0, 'CurrentFigure', o.data.prob_figHandle);
end

f = figure(figid);
colormap(aux.inferno);
set(f, 'name', ['Cond. prob. as matrix, ', o.name]);
imagesc(o.data.prob_srcPrj);
title(['P_{source}(target) for ', o.name]);
xlabel('Source site');
xticks(cellfun(@mean, o.srcRegInd));
xticklabels(o.srcRegName);
ylabel('Projection site');
yticks(cellfun(@mean, o.prjRegInd));
yticklabels(o.prjRegName);
B = colorbar;
title(B, 'Prob');

f = figure(figid + 1);
colormap(aux.pmkmp(256, 'IsoL'));
set(f, 'name', ['PCA space, ', o.name]);
switch DIMRED_TYPE
    case 'nnmf'
        aux.scattermat(o.data.prob_nnmf_sliScore, 25, (1:o.nSrcSli)', ...
            'filled');
        title(['NNMF for joint prob. for ', o.name]);
        if isfield(o.data, 'prob_nnmfFit')
            hold('on');
            for h = 1:size(o.data.prob_nnmfFit.vel, 2)
                vel = o.data.prob_nnmfFit.vel(:, h)';
                int = o.data.prob_nnmfFit.int(:, h)';
                par = o.data.prob_nnmfFit.par(:, h);
                plot3( ...
                    int(1) + vel(1) * par, ...
                    int(2) + vel(2) * par, ...
                    int(3) + vel(3) * par, ...
                    'Color', .5 * ones(1, 3), 'LineWidth', 2, ...
                    'LineStyle', STYLE{mod(h - 1, length(STYLE)) + 1});
            end
        end
    case 'pca'
        aux.scattermat(o.data.prob_pca_sliScore, 25, (1:o.nSrcSli)', ...
            'filled');
        title(['PCA for joint prob. for ', o.name]);
        if isfield(o.data, 'prob_pcaFit')
            hold('on');
            for h = 1:size(o.data.prob_pcaFit.vel, 2)
                vel = o.data.prob_pcaFit.vel(:, h)';
                int = o.data.prob_pcaFit.int(:, h)';
                par = o.data.prob_pcaFit.par(:, h);
                plot3( ...
                    int(1) + vel(1) * par, ...
                    int(2) + vel(2) * par, ...
                    int(3) + vel(3) * par, ...
                    'Color', .5 * ones(1, 3), 'LineWidth', 2, ...
                    'LineStyle', STYLE{mod(h - 1, length(STYLE)) + 1});
            end
        end
end
B = colorbar;
title(B, 'Slice no');
axis('image');

f = figure(figid + 2);
set(f, 'name', ['Cond. probabilites plot, ', o.name]);
P = gobjects(o.nPrjSli, 1);
COL_REG = lines(o.nPrjSli);
% Plot baselines and probability distributions
for s = 1:o.nPrjSli
    P(s) = plot(o.data.prob_srcPrj(s, :), ...
        '-', 'LineWidth', 2, 'Color', COL_REG(s, :));
    hold('on');
end

% Plot the patterns
LEG = cell(size(o.data.prob_pcaFit.vel, 2) * o.nPrjReg, 1);
for h = 1:size(o.data.prob_pcaFit.vel, 2)
    x = o.data.prob_fitSli(h, :);
    %FIXTHIS
    if strcmp(o.name, 'PC injection') && (h == 1)
        x = x(end:-1:1);
    end
    for s = 1:o.nPrjSli
        r = o.getPrjSliReg(s);
        P(h * o.nPrjSli + s) = plot( x, ...
            o.data.prob_pcaFit.prob{h}(s, :), ...
            STYLE{mod(h - 1, length(STYLE)) + 1}, 'LineWidth', 2, ...
            'Color', COL_REG(s, :));
        LEG{(h - 1)* o.nPrjSli + s} = [num2str(h), ' section, region ', ...
            o.prjRegName{r}, ' slice ', num2str(s)];
    end
end
hold('off');
title(['Conditional probability for ', o.name]);
xlabel('Source site');
xlim([1; o.nSrcSli]);
ylabel('Probability');
YLIM = ylim;
ylim([max(YLIM(1), 0), YLIM(2)]);
grid('on');
legend(P, [ ...
    cellfun(@(x) ['P(', x, '|source)'], o.prjRegName, 'UniformOutput', 0); ...
    LEG], 'Location', 'eastoutside');
%...
%     arrayfun(@(x) ['Fit region ', num2str(x)], ...
%     (1:size(o.data.prob_pcaFit.vel, 2))', 'UniformOutput', 0) ]);


end