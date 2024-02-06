% integrateTangentialJump Integrate tangential jump of given function over edges
%   numerically.
%
%   jump = integrateTangentialJump(fun, qr) returns tangential jump across edges using
%       1D quadrature rule qr. The sign of jump is "left - right". The
%       orientation of inner edges is arbitrary. For boundary edges "right" is
%       outside of the domain and this side is taken to be 0..
%   
%   jump = integrateTangentialJump(fun, qr,  [postprocHandle, postprocFuncs, idx], ...)
%       Additionally apply function handle postprocHandle to values at edges
%       idx. The arguments for the function handle are given in the cell array
%       postprocFuncs, where the first argument is reserved for the jump.
%       Multiple triples are possible.
%
%   See also: CompositeFunction, integrateJump

function edgeVal = integrateTangentialJump(fun, qr, varargin)

tangential = TangentialVector(fun.mesh);
edgeVal = integrateJump(fun, qr, ...
    @(j,t) vectorProduct(j, t, [size(j, Dim.Vector)/2,2], [2,1]), {tangential}, ':', ...
    varargin{:});

end