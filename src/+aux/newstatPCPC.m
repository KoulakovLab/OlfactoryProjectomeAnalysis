function newstatPC(set, BOOT)
%NEWSTATOB Three things that are, ugh . . .
% 1. Spearman correlation
% 2. Bootstrap slices for fit slope; check how often the slope changes
%   signs
% 3. Shuffle slices to calculate pvalue wrt shuffling
%#ok<*UNRCH>
rng('default');
USE_NORMALPDF = true;
SKIP = 2;
PC_SOMA = 1;

% Store variables here
outstr = struct();

% Useful
n_brc = size(set.nBrc, 1);
n_sli = set.nSrcSli;
n_apc = set.nSrcRegSli(1);
n_ppc = set.nSrcRegSli(2);

% Get the PC image to work with
idmat = zeros(size(set.srcImg));
idmat(sub2ind(size(set.srcImg), (1:set.nBrc)', set.brcSrcVis)) = 1;
sto_mask = idmat;
for s = 1:PC_SOMA
  sto_mask = sto_mask ...
    + [idmat(:, (1+s):end), zeros(size(idmat, 1), s)] ...
    + [zeros(size(idmat, 1), s), idmat(:, 1:(end-s))];
end
conmat = set.srcImg .* (1 - sto_mask);

idmat_reg = [sum(idmat(:, SKIP:n_apc), 2), ...
  sum(idmat(:, (n_apc + 1):(end - SKIP)), 2)];
conmat_reg = [sum(conmat(:, SKIP:n_apc), 2), ...
  sum(conmat(:, (n_apc + 1):(end - SKIP)), 2)];

% Calculate the joint probability matrix
outstr.conProb = aux.conProb(conmat, idmat);
outstr.conProbReg = aux.conProb(conmat_reg, idmat_reg);

% Calculate p-values on the elements
cp_boot = zeros(n_sli, n_sli, BOOT);
cr_boot = zeros(2, 2, BOOT);
cp_elem = ceil(n_brc * rand(BOOT, n_brc));
for i = 1:BOOT
  cp_boot(:, :, i) = aux.conProb(conmat(cp_elem(i, :), :), ...
    idmat(cp_elem(i, :), :));
  cr_boot(:, :, i) = aux.conProb(conmat_reg(cp_elem(i, :), :), ...
    idmat_reg(cp_elem(i, :), :));
end
cp_boot_flat = reshape(cp_boot, [], BOOT);
cr_boot_flat = reshape(cr_boot, [], BOOT);
[~, cp_boot_flat_pval] = ttest( ...
  cp_boot_flat' - 1 / (n_sli - 2 * PC_SOMA - 1));
[~, cr_boot_flat_pval] = ttest( ...
  cr_boot_flat' - 1 / 2);
outstr.conProb_pval = reshape(cp_boot_flat_pval, n_sli, n_sli);
outstr.conProbReg_pval = reshape(cr_boot_flat_pval, 2, 2);

set.data.PCPC = outstr;
rng('shuffle');

end

