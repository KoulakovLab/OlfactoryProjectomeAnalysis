function [crnDimTrp_cor, trp_cntId, cntDim_cor] = gen3ControlPts(N, P)
  %GENCONTROLPTS Generate triples of control points for classification
  % Used in a brute force search for voronoi triangulation of a 2-simplex
  %   embedded in 2D. (equilateral triangle, centered on origin and side
  %   length of sqrt(2))
  % Each control triple is a unique voronoi triangulation with the 
  %   triple-point (circumcenter) located inside the 2-simplex, and each
  %   node of the simplex is within a unique voronoi region.
  % INPUT:
  %   N: Circumcenter search number. There will be N*(N-1)/2 circumcenters
  %       tested.
  %   P: Angle search number. Control points will be searched, with a 
  %       minimum angular seperation of 2pi/P
  % OUTPUT:
  %   crnDimTrp_cor: Array of size (Corner=1:3, Dim=1:2, Triple)
  %       Gives the <Dim> coordinate of the point identifying <Corner>
  %   trp_cntId: The ID of the circumcenter of this triple
  %   cntDim_cor: The <Dim> coordinate of the <Id>'th circumcenter.
  
  % The distance of control points from the circumcenter
  RAD = sqrt(2) / 3;
  
  % Get the ordering; and convert to clockwise from 0 to 2pi
  [~, ~, crn] = aux.genSimplex(2);,
  crn_ang = atan2(crn(:, 2), crn(:, 1));
  crn_ang = crn_ang + (2*pi) * (crn_ang < 0);
  [~, ord] = sort(crn_ang);
  crn = crn(ord, :);
  
  % Get all circumCeNTer candidate COoRdinates
  cntDim_cor = aux.genEqlTriPts(N - 1) * sqrt(2) * (N - 1) / N;
  cN = size(cntDim_cor, 1);
  
  % Create list of ANGular DIVisors
  div_ang = ((0:(P-1))') * (2*pi/P);
  % Find the midway angle between two points (from->to)
  divDiv_angDelta = (div_ang') - (div_ang);
  divDiv_angDelta = divDiv_angDelta + (2*pi) * (divDiv_angDelta < 0);
  divDiv_angMean = div_ang + .5 * divDiv_angDelta;
  divDiv_angMean = divDiv_angMean - (2*pi) * (divDiv_angMean >= (2*pi));
  
  % Get the angle of the corners wrt the circumcenter
  cntCrnDim_dis = reshape(crn, 1, 3, 2) - reshape(cntDim_cor, cN, 1, 2);
  cntCrn_ang = atan2(cntCrnDim_dis(:, :, 2), cntCrnDim_dis(:, :, 1));
  cntCrn_ang = cntCrn_ang + (2*pi) * (cntCrn_ang < 0);
  
  % Find the (+)angle from center to the divisors for each center
  divCntCrn_ang = div_ang - reshape(cntCrn_ang, 1, cN, 3);
  divCntCrn_ang = divCntCrn_ang + (2*pi) * (divCntCrn_ang < 0);
  
  % Find the ID of the control point that the divisor is right after
  [~, divCnt_crnId] = min(divCntCrn_ang, [], 3);
  
  % Per each control point; break the divider index into ones following after
  cntCrn_divIds = arrayfun(@(x, y) find(divCnt_crnId(:, x) == y), ...
    ((1:cN)') .* ones(1, 3), ones(cN, 1) .* (1:3), 'UniformOutput', false);
  
  % For each control point; find 3 divisors such that they come rigt after
  % the nth corner
  cnt_trpCrn_divIds = cellfun(@(x, y, z) [ ...
    repmat(  repmat(x, length(y), 1), length(z), 1), ...
    repelem( repmat(y, length(z), 1), length(x), 1), ...
    repelem(repelem(z, length(y), 1), length(x), 1)], ...
    cntCrn_divIds(:, 1), cntCrn_divIds(:, 2), cntCrn_divIds(:, 3), ...
    'UniformOutput', false);
  
  % Get the mean angle between current and previous divider
  cnt_trpCrn_divAngMean = cellfun(@(x) ...
    divDiv_angMean(sub2ind([P, P], x(:, mod((1:3) - 2, 3) + 1), x)), ...
    cnt_trpCrn_divIds, 'UniformOutput', false);
  
  % Get the control coordinates
  cnt_crnDimTrp_cor = cellfun(@(x, y) [ ...
    cos(reshape(x', 3, 1, [])) + y(1), ...
    sin(reshape(x', 3, 1, [])) + y(2)], ...
    cnt_trpCrn_divAngMean, num2cell(cntDim_cor, 2), 'UniformOutput', false);
  
  % Create the array that houses all the triples
  crnDimTrp_cor = cat(3, cnt_crnDimTrp_cor{:});
  
  % Create the index array for circumcenter coordinates per triple
  cnt_numTrp = cellfun(@(x) size(x, 3), cnt_crnDimTrp_cor);
  trp_cntId = cell2mat(arrayfun(@(x, y) x * ones(y, 1), ...
    (1:cN)', cnt_numTrp, 'UniformOutput', false));
  
  % Reorder everything such that the initial corner order is satisfied
  crnDimTrp_cor(ord, :, :) = crnDimTrp_cor;
end

