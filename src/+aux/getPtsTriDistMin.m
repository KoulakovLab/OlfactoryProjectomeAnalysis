function DIS = getPtsTriDistMin(REF, TRI)
  %GETPTSTRIDISTMIN Get the shortest distance from a point to a triangle
  dis = zeros(3, 1);
  for cc = 1:3
    % Get the other two points
    c = mod(cc + (1:2) - 1, 3) + 1;
    % Get the delta vector, and the mean vector of these two points
    del = TRI(c(2), :) - TRI(c(1), :);
    avg = mean(TRI(c, :), 1);
    % Get the distance from end points to the referance point
    dif = REF - TRI(c, :);
    dds = sqrt(sum(dif .^ 2, 2));
    % Check sign of cosines
    scs = (dif ./ dds) * ((del ./ sqrt(sum(del .^ 2, 2)))');
    if prod(scs) >= 0
      % If the referance lies outside the plane (same sign for cosine between
      % delta of two triangle ends and the vectors from ends to referance) get
      % the minimum distance to an edge
      dis(cc) = min(dds);
    else
      % Solve for the point on the line connecting two edges and is 
      % closest to the referance point
      clo = [avg(1)*del(2) - avg(2)*del(1), del * (REF')] ...
      / [del(2), del(1); -del(1), del(2)];
    % Get the distance between closest point and referance
      dis(cc) = sqrt(sum((clo - REF) .^ 2, 2));
    end
  end
  DIS = min(dis);
end

