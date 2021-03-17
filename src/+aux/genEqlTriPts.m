function PTS = genEqlTriPts(N)
  % GENEQLTRIPTS Generate equally spaced points on the equilateral triangle
  % Each edge will contain N+1 points; and is of length 1
  % Total of (N+2)(N+1)/2 pts
  if (~isscalar(N)) || (round(N(1)) ~= N(1)) || (N(1) <= 0) 
    error('Need positive integer scalar N');
  end
  PTS_AGG = cell(N+1, 1);
  % Generate points on the bottom row; y=-1/2sqrt(3), x=[-1/2->1/2]
  PTS_AGG{1} = [(linspace(-1/2, 1/2, N+1)'), (-.5/sqrt(3)) * ones(N+1, 1)];
  VEC_DIS = [-1/2, sqrt(3)/2] ./ N;
  % Generate the top tier
  for i = 1:N
    PTS_AGG{i + 1} = VEC_DIS + PTS_AGG{i}(2:end, :);
  end
  % Concatanate
  PTS = vertcat(PTS_AGG{:});
end

