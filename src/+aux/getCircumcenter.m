function C = getCircumcenter(X)
  %GETCIRCUMCENTER Calculate circumcenter of 3 points in 2D space
  if any(size(X) ~= [3, 2])
    error('Need 3 points in 2D');
  end
  del_1 = X(2, :) - X(1, :);
  del_2 = X(3, :) - X(2, :);
  avg_1 = mean(X(1:2, :), 1)';
  avg_2 = mean(X(2:3, :), 1)';
  C =  [del_1 * avg_1, del_2 * avg_2] / ([del_1; del_2]');
end

