% diffusionPart Computes data for diffusion term (a \nabla u) * \nabla v.
%
%   See also: BilinearForm

function val = diffusionPart(a, Dphi)

na = sqrt(size(a, Dim.Vector));
nD = size(Dphi, Dim.Vector) / 2;
n = 2 * nD / na;
assert(na == 1 || na == 2, 'BilinearForm:wrongCoefficientDimension', ...
    'a must be either scalar or [2,2]-matrix.')

val = vectorProduct(Dphi, a, [n,na], [na,na]');
val = vectorProduct(Dphi, val, [nD,2], [nD,2]');

end
