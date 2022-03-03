% reactionPart Computes data for c * u * v.
%
%   See also: BilinearForm

function val = reactionPart(c, phi)

nc = size(c, Dim.Vector);
nD = size(phi, Dim.Vector);
assert(nc == 1, 'BilinearForm:wrongCoefficientDimension', 'c must be scalar.')

val = vectorProduct(phi, c, [nD,1], [1,1]);
val = vectorProduct(phi, val, [nD, 1], [nD,1]');

end
