function C = rand_sphere(N, D)
% RAND_SPHERE Generate uniformly random points on the sphere inside D dim.
%   The sphere is S_{D-1}, and resides embedded in R^D

% Input checking
assert(isscalar(D), 'Dimensions must be scalar');
assert(round(D) == D, 'Dimension must be integer');
assert(D > 1, 'Dimension must be >1');
assert(isscalar(N), 'Point number must be scalar');
assert(round(N) == N, 'Point number must be integer');
assert(N > 0, 'Number of points must be >0');

% Generate points on sphere, with none at the origin
C = zeros(N, D);
R = zeros(N, 1);
while any(R == 0)
    I = (R == 0);
    C(I, :) = randn(sum(I), D);
    R(I) = sqrt(sum(C(I, :) .^ 2, 2));
end

% Project to S_(D-1)
C = C ./ R;

end