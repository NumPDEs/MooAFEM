% asVector Return vector which is a column-wise representation of input.
%
%   y = asVector(x) takes an nxm matrix x and returns the n*m column-wise
%       representation. This conforms with the memory model of mooAFEM.

function y = asVector(x)

y = x(:);

end