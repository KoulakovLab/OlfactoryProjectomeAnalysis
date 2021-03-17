function y = limexpn(A, x)
% LIMEXPN Exponent, limited to cross the points (0, 0) and (1, 1).
%   This function evaluates y = (e^(Ax) - 1)/(e^A - 1), but takes care of
%   the removable discontinuity at A=0 by evaluating to y=x at A=0.
% INPUTS:
%   A:  Parameter that controls how sharp the exponent is. At =0, the 
%       exponent is a straight line.
%   x:  Values that the funciton is evaluated at.
% OUTPUTS:
%   y:  Result of the function

if A == 0
    y = x;
else
    y = (exp(A*x) - 1) ./ (exp(A) - 1);
end

end