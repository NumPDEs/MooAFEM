% integrateElement Integrate given function numerically element-wise.
%
%   val = integrateElement(fun, qr) integrates function fun over all elements
%       with quadrature rule qr.
%
%   val = integrateElement(fun, qr, idx) integrates only over elements given by
%       idx.

function val = integrateElement(fun, qr, idx)

arguments
    fun (1,1) Evaluable
    qr (1,1) QuadratureRule
    idx (1,:) {mustBeIndexVector} = ':'
end

assert(qr.dim == 2, 'Quadrature rule must be 2-dimensional.')

% accumulate weighted contributions of quadrature nodes
val = sum(eval(fun, qr.bary, idx) .* reshape(qr.weights, 1, 1, []), Dim.QuadratureNodes);

% apply correct scaling
trafo = getAffineTransformation(fun.mesh);
val = val .* reshape(trafo.area(idx), 1, []);

end
