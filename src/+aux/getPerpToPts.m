function C = getPerpToPts(P, M)
  %GETPERPTOPTS Get intersection of bisector to connecting line
  if any(size(P) ~= [2, 2], 'all') || any(size(M) ~= [2, 2], 'all')
    error('Must give 2 sets of 2 points in 2D');
  end
  % Quantities
  DEL_P = P(2, :) - P(1, :);
  AVG_P = mean(P, 1)';
  DEL_M = M(2, :) - M(1, :);
  AVG_M = mean(M, 1)';
  % Solution
  C = [DEL_P * AVG_P, (DEL_M(1) * AVG_M(2) - DEL_M(2) * AVG_M(1))] ...
    / [DEL_P(1), -DEL_M(2); DEL_P(2), DEL_M(1)];
end

