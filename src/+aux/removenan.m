function Y = removenan(X)
%REMOVENAN Takes in a vector; and removes nan values.
if ~isvector(X)
    warning('Non-vector input will be vectorized');
end
Y = X(~isnan(X(:)));
end

