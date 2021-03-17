function getIprSrc(o)
%CALCIPR Calculate the barcode IPR's of the data.
%   IPR stands for inverse participation ratio; which is the inverse of the
%   fractional sum of squares. (If one-hot-vector, it is 1. If uniformly
%   diffuse vector, it is D: vector dimension)

% The IPR function; do note that it gives 0 on a zero vector
IPR = @(x) (sum(x, 2) .^ 2) ./ (sum(x .^2, 2) + (sum(x, 2) == 0));

% Calculate the entire IPR of the images
o.brc_src_ipr = cellfun(IPR, o.src_img, 'UniformOutput', 0);
o.brc_prj_ipr = cellfun(IPR, o.prj_img, 'UniformOutput', 0);

% Calculate the barcode to region matrix of each individual region
o.brc_reg_src_ipr = cellfun(@(x) cell2mat( ...
    arrayfun(@(i) IPR(x(:, o.src_reg_ind{i})), ...
    1:o.num_src_reg, 'UniformOutput', 0)), o.src_img, 'UniformOutput', 0);
o.brc_reg_prj_ipr = cellfun(@(x) cell2mat( ...
    arrayfun(@(i) IPR(x(:, o.src_prj_ind{i})), ...
    1:o.num_prj_reg, 'UniformOutput', 0)), o.prj_img, 'UniformOutput', 0);

end

