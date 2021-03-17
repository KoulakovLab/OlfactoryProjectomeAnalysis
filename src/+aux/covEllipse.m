function h = covEllipse(mu, Sigma, p)
  % covEllipse draws the covariance ellipse
  if exist('p', 'var')
    s = -2 * log(1 - p);
  else
    s = -2 * log(.1);
  end
  
  [V, D] = eig(Sigma * s);
  
  t = linspace(0, 2 * pi);
  a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
  
  plot(a(1, :) + mu(1), a(2, :) + mu(2));
end