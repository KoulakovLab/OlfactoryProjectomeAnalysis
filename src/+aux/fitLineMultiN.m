function [v, x, t] = fitLineMultiN(q)
%FITLINEMULTIN Fits a line to multidimensional data
% INPUTS
% * q:  Points in space; expect in format (instance, dimension)
% - If left empty; a random set is generated
%   Optional arguments are;
% OUTPUTS
% * x:  Starting point for parametrization
% * v:  Velocity of parametrization
% * t:  The distance, from x, of each point's projection on the line

% Do the fit
v = pca(q, 'NumComponents', 1)';
x = mean(q, 1);
t = (q - x) * (v');

% Make it so that t ranges from -1 to +1; generally in increasing fashion
%---ORDERING: if ~slope of t vs index is -; flip it around
% Least squares slope of t vs t_ind is comp.;
% m = [xy] - [x][y] / [x^2] - [x]^2; m ~ [xy] - (N+!)/2 [y]
% sum(n^2, 1, N) = N(N+1)(2N+1)/6
% sum(n^1, 1, N) = N(N+1)/2
% Slope > 0 iff [xy] > (N+1)/2 * [y]
if (t * (1:length(t))) > (sum(t) * (length(t) + 1)/2)
    v = -v;
    t = -t;
end
%---ORIGIN: Make the t-values centered
dt = .5 * (max(t) + min(t));
t = t - dt;
x = x + dt * v;
%---SCALING: Make t range from +-1
r = .5 * (max(t) - min(t));
t = t / r;
v = v * r;

end