function flattenRegions(o, ins)
%FLATTENREGIONS Sums accross regions
%   Converts multi-slice information to integrated accross regions
%   USAGE
%       <instance>.flattenRegions(<specifier>)
%   INPUTS
%       instance:   An instance of the mapseqDataAligned class
%       specifier:  Either 'source' (or src), 'projection' (or prj) or 'all'

% Flatten the source slices
if any(strcmp(ins, {'all', 'source', 'src'}))
    for b = 1:o.nBrc
        newImg = zeros(o.nBrc(b), o.nSrcReg);
        for r = 1:o.nSrcReg
            newImg(:, r) = sum(o.srcImg{b}(:, o.srcRegInd{r}), 2);
        end
        o.srcImg{b} = newImg;
        o.nSrcRegSli = 1 + (0 .* o.nSrcRegSli);
    end
end

% Flatten the projection slices
if any(strcmp(ins, {'all', 'projection', 'prj'}))
    for b = 1:o.nBrc
        newImg = zeros(o.nBrc(b), o.nPrjReg);
        for r = 1:o.nPrjReg
            newImg(:, r) = sum(o.prjImg{b}(:, o.prjRegInd{r}), 2);
        end
        o.prjImg{b} = newImg;
        o.nPrjRegSli = 1 + (0 .* o.nPrjRegSli);
    end
end

end

