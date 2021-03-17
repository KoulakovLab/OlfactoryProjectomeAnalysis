function calcSrcPrjProb(o, lin)
%calcSrcPrjProb Do the bayesian matrices between source and proj.
if exist('lin', 'var')
    LINFIT = lin;
else
    LINFIT = 2;
end

% Probability; first clean data
o.cleanEmpty;

% Event is rolony selection;

% Uniform sampling per barcode;
% P(barcode): M(b, 1)
%   * rolony will be of barcode b
o.data.prob_brc = ones(o.nBrc, 1) / o.nBrc;

% Unitarity(?): sum(p_i) = 1
% P_barcode(slice): M(s, b)
%   * (given the rolony selected is of barcode b)
%   * rolony will be in slice s
o.data.prob_brcSrc = fillmissing( ...
    (o.srcImg ./ sum(o.srcImg, 2))', ...
    'constant', 0);
o.data.prob_brcPrj = fillmissing( ...
    (o.prjImg ./ sum(o.prjImg, 2))', ...
    'constant', 0);

% Chain rule; p(A) = p(B)*p(A|B)
% P(slice): M(s, 1) 
%   * rolony selected will be in slice s
o.data.prob_src = o.data.prob_brcSrc * o.data.prob_brc;
o.data.prob_prj = o.data.prob_brcPrj * o.data.prob_brc;

% Bayes rule; p(A)*p(B|A) = p(B)*p(A|B) [ := p(A,B) ] hence
%   p(B|A) = p(B) * p(A|B) / p(A)
% P_slice(barcode): M(b, s)
%   * (given the rolony selected is in slice s)
%   * rolony will belong to barcode b.
o.data.prob_srcBrc = fillmissing( ...
    (o.data.prob_brcSrc') .* o.data.prob_brc ./ (o.data.prob_src'), ...
    'constant', 0);
o.data.prob_prjBrc = fillmissing( ...
    (o.data.prob_brcPrj') .* o.data.prob_brc ./ (o.data.prob_prj'), ...
    'constant', 0);

% Chain rule: p(A|B) = p(C|B)*p(A|C)
% P(source|projection): M(s_prj,s_src)
%   * (given the rolony selected is of slice s_src in source)
%   rolony will be in slice s_prj in projection
o.data.prob_prjSrc = o.data.prob_brcSrc * o.data.prob_prjBrc;
o.data.prob_srcPrj = o.data.prob_brcPrj * o.data.prob_srcBrc;

%-------------------------------%
%-----DIMENSIONAL REDUCTION-----%
%-------------------------------%

% COMP: components are (region, lowdim)
% SCORE: scores are (barcode, lowdim)

%-----PCA-----%

% Do PCA on the slices depending on their probabilities
%   The probabilities do add up to 1, but not exactly (+-.2%); to account
%   for this; just do PCA on the full data; but recover N-1 components

[o.data.prob_pca_comp, o.data.prob_pca_sliScore] = pca( ...
    o.data.prob_srcPrj', ...
    'NumComponents', min(3, o.nSrcSli - 1));

%-----NNMF-----%
% Do non-negative matrix factorization into 3D; representative of the
%   underlying hypothesis.
[o.data.prob_nnmf_comp, o.data.prob_nnmf_sliScore] = nnmf( ...
    o.data.prob_srcPrj, min(3, o.nSrcSli - 1));
o.data.prob_nnmf_sliScore = o.data.prob_nnmf_sliScore';
% Scale so that the L1 norm of the components are 1
o.data.prob_nnmf_sliScore = o.data.prob_nnmf_sliScore .* ...
    sum(o.data.prob_nnmf_comp, 1);
o.data.prob_nnmf_comp = fillmissing(o.data.prob_nnmf_comp ./ ...
    sum(o.data.prob_nnmf_comp, 1), ...
    'constant', 0);

%--------------------%
%-----LINEAR FIT-----%
%--------------------%

if LINFIT
    % Get the different cut slices
    o.data.prob_fitSli = round(((1:LINFIT)') * o.nSrcSli / LINFIT);
    o.data.prob_fitSli = [ ...
        [1; 1 + (o.data.prob_fitSli(1:(end-1)))], o.data.prob_fitSli];
    
    % NNMF fit
    o.data.prob_nnmfFit = struct;
    o.data.prob_nnmfFit.vel = zeros(size(o.data.prob_nnmf_comp, 2), LINFIT);
    o.data.prob_nnmfFit.int = zeros(size(o.data.prob_nnmf_comp, 2), LINFIT);
    o.data.prob_nnmfFit.par = zeros(2, LINFIT);
    o.data.prob_nnmfFit.prob = cell(LINFIT, 1);
    for h = 1:LINFIT
        % Get coordinate slice
        cor = o.data.prob_nnmf_sliScore( ...
            o.data.prob_fitSli(h, 1):o.data.prob_fitSli(h, 2), :);
        [vel, int, par] = aux.fitLineMultiN(cor);
        % Fill in data
        o.data.prob_nnmfFit.vel(:, h) = vel;
        o.data.prob_nnmfFit.int(:, h) = int;
        o.data.prob_nnmfFit.par(:, h) = par;
        % Reconstruct the slice to probability values
        scores = vel .* (par') + int;
        probs = o.data.prob_nnmf_comp * scores;
        o.data.prob_nnmfFit.prob{h} = probs;
    end
    
    % PCA
    o.data.prob_pcaFit = struct;
    o.data.prob_pcaFit.vel = zeros(size(o.data.prob_pca_comp, 2), LINFIT);
    o.data.prob_pcaFit.int = zeros(size(o.data.prob_pca_comp, 2), LINFIT);
    o.data.prob_pcaFit.par = zeros(2, LINFIT);
    o.data.prob_pcaFit.prob = cell(LINFIT, 1);
    for h = 1:LINFIT
        % Get coordinate slice
        cor = o.data.prob_pca_sliScore( ...
            o.data.prob_fitSli(h, 1):o.data.prob_fitSli(h, 2), :);
        [vel, int, par] = aux.fitLineMultiN(cor);
        % Fill in data
        o.data.prob_pcaFit.vel(:, h) = vel;
        o.data.prob_pcaFit.int(:, h) = int;
        o.data.prob_pcaFit.par(:, h) = par;
        % Reconstruct the slice to probability values
        scores = vel .* (par') + int;
        probs = o.data.prob_pca_comp * scores + ...
             mean(o.data.prob_srcPrj, 2);
        o.data.prob_pcaFit.prob{h} = probs;
    end
end

end