% scalarPart Computes data for f * v.
%
%   See also: LinearForm

function val = scalarPart(f, phi)

nf = size(f, Dim.Vector);
nD = size(phi, Dim.Vector);
assert(nf == 1, 'LinearForm:wrongCoefficientDimension', 'f must be scalar.')

val = vectorProduct(phi, f, [nD,1], [1,1]);

end
