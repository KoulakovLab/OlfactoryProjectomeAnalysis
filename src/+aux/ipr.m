function result = ipr(data)
%IPR Calculates Inverse-Participation-Ratio for a measurement
%   INPUT: data
%       data must be a non-negative matrix of size (N, M>1) where N is the
%       number of observations, and M is the number of data points per
%       observation.
%       * Negative values will be set to 0.
%       * NaN values will be considered 0.
%       * Rows with Inf will be calculated with the non-inf set to 0, and
%       the inf set to 1.
%       * Rows that sum to 0 will have 0 IPR.
%   OUTPUT: result
%       result is (N, 1) vector that contains the values of IPR for each
%       observation in data.
if iscolumn(data)
    error('Unsuitable input; must be size (N, M) with M > 1!');
end
if any(data < 0, 'all')
%     warning('Negative values in data will be set to 0!');
    data(data(:) < 0) = 0;
end
if any(isnan(data), 'all')
%     warning('NaN values in data will be set to 0');
    data(isnan(data(:))) = 0;
end
if any(isinf(data), 'all')
%     warning('Rows with Inf will be turned to one-shot measurements');
    data(any(isinf(data), 2), :) = isinf(data(any(isinf(data), 2), :));
end
IPRFUN = @(x) (sum(x, 2) .^ 2) ./ (sum(x .^2, 2) + (sum(x, 2) == 0));
result = IPRFUN(data);

end

