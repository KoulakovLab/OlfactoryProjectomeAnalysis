function trimPrjRegions(obj, trim)
%TRIMREGIONS Trim TO the specified regions in projection end

if ischar(trim)
    trim = {trim};
end
if iscell(trim)
    % Get correspondant indices
    trim = cellfun(@(x) find(strcmpi(x, obj.prjRegName), 1), trim);
elseif isvector(trim)
    if (~all((round(trim) == trim))) || any(trim <= 0) || any(trim > obj.nPrjReg)
        error('Invalid indices');
    end
end


obj.prjImg = obj.prjImg(:, horzcat(obj.prjRegInd{trim}));
if ~empty(obj.different_prjRegSum)
  obj.prjRegSum = obj.prjRegSum(:, trim);
end

obj.prjRegName = obj.prjRegName(trim);
obj.nPrjRegSli = obj.nPrjRegSli(trim);

end

