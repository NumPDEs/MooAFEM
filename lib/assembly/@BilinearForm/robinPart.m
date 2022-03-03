% robinPart Computes data for robin * u * v on Robin boundary.
%
%   See also: BilinearForm

function val = robinPart(robin, phi)

nr = size(robin, Dim.Vector);
nD = size(phi, Dim.Vector);
assert(nr == 1, 'BilinearForm:wrongCoefficientDimension', ...
    'Robin coefficient must be scalar.')

val = vectorProduct(phi, robin, [nD,1], [1,1]);
val = vectorProduct(phi, val, [nD, 1], [nD,1]');

end
