% neumannPart Computes data for g * v on Neumann boundary.
%
%   See also: LinearForm

function val = neumannPart(neumann, phi)

nf = size(neumann, Dim.Vector);
nD = size(phi, Dim.Vector);
assert(nf == 1, 'LinearForm:wrongCoefficientDimension', ...
    'Neumann function handle must be scalar.')

val = vectorProduct(phi, neumann, [nD,1], [1,1]);

end
