% integrateJump Integrate jump of given function over edges numerically.
%
%   jump = integrateJump(fun, qr) returns jump of fun across edges using 1D
%       quadrature rule qr. The sign of jump is "left - right". The orientation
%       of inner edges is arbitrary. For boundary edges "right" is outside of
%       the domain and this side is taken to be 0.
%   
%   jump = integrateJump(fun, qr, [postprocHandle, postprocFuncs, idx], ...)
%       Additionally apply function handle postprocHandle to values at edges
%       idx. The arguments for the function handle are given in the cell array
%       postprocFuncs, where the first argument is reserved for the jump.
%       Multiple triples are possible.
%
%   See also: CompositeFunction, integrateNormalJump

function jump = integrateJump(fun, qr, postprocHandle, postprocFuncs, idx)

arguments
    fun (1,1) Evaluable
    qr (1,1) QuadratureRule
end

arguments (Repeating)
    postprocHandle (1,1) function_handle
    postprocFuncs (1,:) cell
    idx (1,:) {mustBeIndexVector}
end

assert(qr.dim == 1, 'Quadrature rule must be 1-dimensional.')

% loop over edges with local numbering 1->2, 2->3, and 3->1 and do quadrature
mesh = fun.mesh;
for i = 1:3
    qrEl = edge2elementQr(qr, i);
    val = eval(fun, qrEl.bary);
    if i == 1
        sizeVal = size(val, Dim.Vector);
        jump = zeros(sizeVal, mesh.nEdges, qr.nNodes);
    end
    
    right = mesh.flipEdges(i,:);
    elementIsLeft = mesh.element2edges(i, ~right);
    elementIsRight = mesh.element2edges(i, right);
    
    jump(:, elementIsLeft, :) = ...
        jump(:, elementIsLeft, :) + broadcast(val, ~right);
    jump(:, elementIsRight, :) = ...
        jump(:, elementIsRight, :) - flip(broadcast(val, right), 3);
end

% apply correct scaling and normal vector / postprocessing if desired
trafo = getAffineTransformation(mesh);
for k = 1:numel(postprocHandle)
    edgeVals = cellfun(@(x) evalEdge(x, qr.bary, idx{k}), postprocFuncs{k}, 'UniformOutput', false);
    if isequal(idx{k}, ':')
        jump = postprocHandle{k}(jump, edgeVals{:});
    else
        jump(:,idx{k},:) = postprocHandle{k}(jump(:,idx{k},:), edgeVals{:});
    end
end
jump = sum(jump.*reshape(qr.weights, 1, 1, []), Dim.QuadratureNodes) .* trafo.ds;

end

%% local function
function res = broadcast(val, idx)
    if size(val, Dim.Elements) == 1
        res = val;
        return
    end
    res = val(:,idx,:);
end
