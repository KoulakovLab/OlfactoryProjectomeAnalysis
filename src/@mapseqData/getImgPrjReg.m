function prjImg = getImgPrjReg(o, reg)
%GETIMGPRJREG Gets the image from the specified region on the source array

% Convert to cell array if just char
if ischar(reg)
    reg = {reg};
end

% Find indices if string input(s) are given
if iscell(reg)
    ind = cellfun(@(x) find(strcmp(x, o.prjRegName), 1), reg);
else
    ind = reg;
end

% Get indices
lin = horzcat(o.prjRegInd{ind});

% Output the sliced array
prjImg = o.prjImg(:, lin);

end

