function cleanEmpty(o)
%CLEANEMPTY Cleans up empty/invalid barcodes and slices

% Empty barcodes
brcEmp = (sum(o.srcImg, 2) == 0) & (sum(o.prjImg, 2) == 0);
o.srcImg(brcEmp, :) = [];
o.prjImg(brcEmp, :) = [];

% Empty source slices
srcSliEmp = find(arrayfun(@(x) (x == 0) | isnan(x), sum(o.srcImg, 1)));
srcRegEmp = arrayfun(@(x) find(x >= cumsum(o.nSrcRegSli), 1), srcSliEmp);
o.nSrcRegSli(srcRegEmp) = o.nSrcRegSli(srcRegEmp) - 1;
o.srcImg(:, srcSliEmp) = [];

% Empty projection slices
prjSliEmp = find(arrayfun(@(x) (x == 0) | isnan(x), sum(o.prjImg, 1)));
prjRegEmp = arrayfun(@(x) find(x >= cumsum(o.nPrjRegSli), 1), prjSliEmp);
o.nPrjRegSli(prjRegEmp) = o.nPrjRegSli(prjRegEmp) - 1;
o.prjImg(:, prjSliEmp) = [];

% Clean empty slices
o.prjRegName(o.nPrjRegSli == 0) = [];
o.nSrcRegSli(o.nPrjRegSli == 0) = [];
o.srcRegName(o.nSrcRegSli == 0) = [];
o.nSrcRegSli(o.nSrcRegSli == 0) = [];

end

