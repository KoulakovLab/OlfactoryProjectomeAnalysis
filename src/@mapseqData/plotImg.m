function hand = plotImg(o, part, sep, dcl)
%PLOTIMG Plot the image of the current injection image

%-----DEFAULTS-----% BEGIN
% Default to drawing both source and projection images
if ~exist('part', 'var')
    part = 'all';
end
% Default to seperating different regions using lines
if ~exist('sep', 'var')
    sep = true;
end
% Default to no different color accross brains, an
if ~exist('dcl', 'var')
    dcl = false;
end
%-----DEFAULTS-----% END

%-----PREPROCESS-----% BEGIN
% Check which image is requested; and get region list with ids
switch part
    case {'prj', 'projection', 'target'}
        IMG = o.prjImg ./ max(o.prjImg, [], 2);
        REGIND = o.prjRegInd;
        REGNAME = o.prjRegName;
    case {'src', 'source', 'injection'}
        IMG = o.srcImg ./ max(o.srcImg, [], 2);
        REGIND = o.srcRegInd;
        REGNAME = o.srcRegName;
    case {'all', 'both', ''}
        IMG = [ ...
            o.srcImg ./ max(o.srcImg, [], 2), ...
            o.prjImg ./ max(o.prjImg, [], 2)];
        REGIND = [ ...
            o.srcRegInd, ...
            cellfun(@(x) x + o.nSrcSli, o.prjRegInd, 'UniformOutput', 0)];
        REGNAME = [ ...
            cellfun(@(x) ['Source: ', x], o.srcRegName, ...
            'UniformOutput', 0), ...
            cellfun(@(x) ['Target: ', x], o.prjRegName, ...
            'UniformOutput', 0)];
    otherwise
        error('Invalid part specification');
end
% Scale image so each barcode is visible throughout
IMG = IMG ./ (max(0, max(IMG, [], 2)) + (max(IMG, [], 2) <= 0));
%-----PREPROCESS-----% END

%-----DRAW-----% BEGIN
if dcl% Multicolor image
    [IDS, ~, CLASS] = unique(o.brcId(:, 1));
    if length(IDS) < 2
        hand = imagesc(IMG);
    else
        CMAP = aux.pmkmp(size(IDS, 1), 'IsoL');
        IMGLAB = o.brcName{1};
        % Recreate the images in 3 channels
        STO = zeros([size(IMG), 3]);
        for c = 1:3
            STO(:, :, c) = IMG .* CMAP(CLASS, c);
        end
        % Draw image
        hand = image(STO);
        % Create color bar for image
        hold('on');
        for b = 1:length(IDS)
            plot(NaN, NaN, 's', ...
                'MarkerSize', 10, ...
                'MarkerFaceColor', CMAP(b, :));
        end
        L = legend(IMGLAB, ...
            'AutoUpdate', 'off', ...
            'Location', 'northeast');
        title(L, 'Sample');
        hold('off');
    end
else                            % Single color image
    hand = imagesc(IMG);
end
%-----DRAW-----% END

%-----AXES-----% BEGIN
% Y axes: barcodes
ylabel('Barcodes');
% X axes: regional slices
switch part
    case {'prj', 'projection', 'target'}
        xlabel('Projection Regions');
    case {'src', 'source', 'injection'}
        xlabel('Source Regions');
    otherwise
        xlabel('Regions');
end
xticks(cellfun(@mean, REGIND));
xtickangle(45);
xticklabels(REGNAME);
if sep
    hold('on');
    for r = 2:length(REGIND)
        plot(ones(1, 2) * (.5 + max(REGIND{r - 1}, [], 'all')), ylim, ...
            ':', 'Color', .5 * ones(1,3), 'LineWidth', 3);
    end
    % Do also source-injection seperation
    if any(strcmp({'all', 'both', ''}, part), 'all')
        plot(ones(1, 2) * (.5 + max(REGIND{o.nSrcReg}, [], 'all')), ...
            ylim, '-', ...
            'Color', .5 * ones(1,3), 'LineWidth', 3);
    end
    hold('off');
end
% Title
switch part
    case {'prj', 'projection', 'target'}
        title([o.name, ', projection site']);
    case {'src', 'source', 'injection'}
        title([o.name, ', injection site']);
    case {'all', 'both', ''}
        title(o.name);
end
%-----AXES-----% END

end

