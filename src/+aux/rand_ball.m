function C = rand_ball(N, D)
% RANDO_BALL Generate uniformly random points within ball in D dim.
%   The ball of question is the B_D, embedded in R^D

% Input checking
assert(isscalar(D), 'Dimensions must be scalar');
assert(round(D) == D, 'Dimension must be integer');
assert(D > 1, 'Dimension must be >1');
assert(isscalar(N), 'Point number must be scalar');
assert(round(N) == N, 'Point number must be integer');
assert(N > 0, 'Number of points must be >0');

% Generate points on sphere, with none at the origin
C = aux.rand_sphere(N, D);

% Generate distances so that distribution is isovolumetric
C = C .* (rand(N, 1) .^ (1/D));

end