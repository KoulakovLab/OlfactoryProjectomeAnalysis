function [fromBaryo, toBaryo, subCor] = genSimplex(N)
%GENSIMPLEX Generate symplex transformations
% Function to tranlate an N-Simplex (N+1) pts from/to baryocentric form
%
% Baryocentric coordinates are when the N+1 points of the simplex are
% expressed in N+1 dimensions; with each vertex corresponding to a unit
% vector among a dimension.
%   t_0 = (1, 0, 0, ..., 0, 0)
%   t_1 = (0, 1, 0, ..., 0, 0)
%   ...
%   t_N = (0, 0, 0, ..., 0, 1)
%
% These functions convert from that sub-space in R^(N+1) to R^(N), while
% preserving distances. (each vertex is a distance sqrt(2) from other)
%   fromBaryo: Convert, from baryocentric 
if (~isscalar(N)) || (round(N(1)) ~= N(1)) || (N(1) <= 0)
  error('N must be scalar integer > 0');
end

% Procedurally generate the simplices
V = cell(N+1, 1);

% 0-simplex has one point in R^0
V{1} = zeros(1, 0);

% Recursively generate the R^(N) subspace of the N-simplex
for n = 1:N
  V{n + 1} = [ ...
    V{n}, (-1 / sqrt(n * (n + 1))) * ones(n, 1); ...
    zeros(1, n - 1), sqrt(n / (n + 1))];
end
subCor = V{N + 1};

% Create the square transformation matrix and the displacement row
D = V{N + 1}(N + 1, :);
V = V{N + 1}(1:N, :) - V{N + 1}(N + 1, :);

% The function to map from baryocentric to lower dimensional subspace
fromBaryo = @(x) (x(:, 1:N) * V) + D;

% The function to map from lower dimensional subspace to baryocentric
toBaryo = @(x) [((x - D) / V), 1 - sum((x - D) / V, 2)];

end

