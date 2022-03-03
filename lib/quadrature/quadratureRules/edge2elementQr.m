% edge2elementQr Convert 1D QuadratureRule on given edge to 2D QuadratureRule.
%
%   qr2d = edge2elementBarycentric(qr1d, edgeInd)

function qr2d = edge2elementQr(qr1d, edgeInd)

arguments
    qr1d (1,1) QuadratureRule
    edgeInd double {mustBeInRange(edgeInd, 1, 3)}
end

newBary = Barycentric2D.from1D(qr1d.bary, edgeInd);
qr2d = QuadratureRule(newBary, qr1d.weights);

end