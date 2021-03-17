function out = red2cyan(n)
%RED2CYAN Prints a colormap of n values; with ramps from red, white and
if ~exist('n', 'var')
    n = 256;
end
N1 = round(n/2);
N2 = n - N1;
% Ramp from red to white
out = [ ...
    hsv2rgb([0.0 * ones(N1, 1), (linspace(1, 0, N1)'), ones(N1, 1)]); ...
    hsv2rgb([0.5 * ones(N2, 1), (linspace(0, 1, N2)'), ones(N2, 1)])];
end

