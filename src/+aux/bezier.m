function OUT = bezier(PTS, SAM)
%BEZIER Given control points, and sampling size outputs the bezier curve
%   Don't provide sample point number to get a function handle that will
%   plot the bezier curve itself.

N = size(PTS, 1);
D = size(PTS, 2);
% Create function handle that will do the bezier curves
COEFF = zeros(1, N);
for n = 1:N
    COEFF(n) = nchoosek(N-1, n-1);
end
EXP1 = ((N-1):-1:0);
EXP2 = (0:1:(N-1));
OUT = @(x) (COEFF ...
    .* ((linspace(1, 0, x)') .^ EXP1) ...
    .* ((linspace(0, 1, x)') .^ EXP2)) * PTS;
% Evaluate to points on the bezier curve if requested.
if exist('SAM', 'var')
    OUT = OUT(SAM);
end

end

