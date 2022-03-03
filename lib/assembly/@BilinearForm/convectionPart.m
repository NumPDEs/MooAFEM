% convectionPart Computes data for convection term b * \nabla u * v.
%
%   See also: BilinearForm

function val = convectionPart(b, Dphi, phi)

nb = size(b, Dim.Vector);
nD = size(phi, Dim.Vector);
assert(nb == 2, 'BilinearForm:wrongCoefficientDimension', 'b must be a 2-vector.')

val = vectorProduct(Dphi, b, [nD,2], [2,1]);
val = vectorProduct(phi, val, [nD,1], [nD,1]');

end
