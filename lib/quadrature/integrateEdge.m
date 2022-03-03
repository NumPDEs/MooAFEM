% integrateEdge Integrate given function numerically edge-wise.
%
%   val = integrateEdge(fun, qr) integrates function fun over all edges with
%       quadrature rule qr.
%
%   val = integrateEdge(fun, qr, idx) integrates only over edges given by idx.

function val = integrateEdge(fun, qr, idx)

arguments
    fun (1,1) Evaluable
    qr (1,1) QuadratureRule
    idx (1,:) {mustBeIndexVector} = ':'
end

assert(qr.dim == 1, 'Quadrature rule must be 1-dimensional.')

% accumulate weighted contributions of quadrature nodes
val = sum(evalEdge(fun, qr.bary, idx) .* reshape(qr.weights, 1, 1, []), Dim.QuadratureNodes);

% apply correct scaling
trafo = getAffineTransformation(fun.mesh);
val = val .* reshape(trafo.ds(idx), 1, []);

end
