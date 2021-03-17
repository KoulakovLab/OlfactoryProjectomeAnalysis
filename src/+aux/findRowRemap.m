function IND = findRowRemap(SRC, TGT)
%FINDROWREMAP Find the index mapping of the rows between two matrices
%   Granted that every row in TGT is found in SRC; then the result is s.t.
%   SRC(IND, :) is the same as TGT

distances = pdist2(TGT, SRC, 'euclidean');
[closestDist, IND] = min(distances, [], 2);

if any(closestDist ~= 0)
  warning('Rows not really equal');
end

end

