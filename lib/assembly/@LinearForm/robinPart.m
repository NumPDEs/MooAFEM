% robinPart Computes data for g * v on Robin boundary.
%
%   See also: LinearForm

function val = robinPart(robin, phi)

nf = size(robin, Dim.Vector);
nD = size(phi, Dim.Vector);
assert(nf == 1, 'LinearForm:wrongCoefficientDimension', ...
    'Robin function handle must be scalar.')

val = vectorProduct(phi, robin, [nD,1], [1,1]);

end
