% vectorPart Computes data for fvec * \nabla v.
%
%   See also: LinearForm

function val = vectorPart(fvec, Dphi)

nf = size(fvec, Dim.Vector);
nD = size(Dphi, Dim.Vector) / 2;
assert(nf == 2, 'LinearForm:wrongCoefficientDimension', 'fvec must be a 2-vector.')

val = vectorProduct(Dphi, fvec, [nD,2], [2,1]);

end
