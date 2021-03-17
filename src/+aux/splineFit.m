function [m, y_int] = splineFit(x, y, x_int)
  %SPLINEFIT given x intervals; fit a linear spline that optimizes LSQ
  %  NOTE: Without the x_int argument; will do regular linear LSQ fit.
  %   INPUT
  %     x: x values of the data
  %     y: y values of the data
  %     x_int: the seperation boundaries
  %   OUTPUT
  %     m: slopes of the individual sections
  %     y_int: the y values of teh interval
  
  % Turn into column vectors
  x = x(:);
  y = y(:);
  
  if ~exist('x_int', 'var')
    x_int = [min(x) - 1e-15; max(x) + 1e-15];
  else
    x_int = x_int(:);
    
    if length(x_int) < 2
      error('Need at least 2 boundaries');
    end
  end
  
  % Get points within partitions created by the boundaries
  partid = (x > (x_int(1:(end - 1))')) & (x <= (x_int(2:end)'));
  partn = sum(partid, 1);
  partd = x_int(2:end) - x_int(1:(end - 1));
  partc = [0, cumsum(partn)];
  
  % matrix defined by x measurements
  A = spalloc(partc(end), length(x_int), partc(end) * 2);
  
  % vector for y measurements
  Y = zeros(partc(end), 1);
  
  if any(partn == 0)
    error('Empty partitions detected');
  end
  
  for i = 1:length(partn)
    % Get index of pts in this partition
    pi = partid(:, i);
    
    if sum(pi) == 0
      error('No points in interval');
    end
    
    % Get x and y values in this partition
    xi = x(pi);
    yi = y(pi);
    
    % Enter values into matrix
    A((partc(i) + 1):partc(i + 1), i) = ...
      1 + (x_int(i) - xi) / partd(i);
    A((partc(i) + 1):partc(i + 1), i + 1) = ...
      0 - (x_int(i) - xi) / partd(i);
    
    % Get y values
    Y((partc(i) + 1):partc(i + 1)) = yi;
    
  end
  
  % Get y intercepts
  y_int = A \ Y;
  
  % Calculate the slopes
  m = (y_int(2:end) - y_int(1:(end - 1))) ./ ...
    (x_int(2:end) - x_int(1:(end - 1)));
  
end

