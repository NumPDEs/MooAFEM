% createInnerEdgeDofs Generate global DOF-numbers associated with inner edge
%   dofs.
%
%   dofs = createInnerEdgeDofs(obj) generates dofs for all edges.
%
%   dofs = createInnerEdgeDofs(obj, idx) generates dofs only for edges given by
%       index vector idx.

function dofs = createInnerEdgeDofs(obj, idx)

arguments
    obj
    idx (1,:) {mustBeIndexVector} = ':'
end

N = obj.mesh.nEdges;
nLocalDofs = getDofConnectivity(obj.finiteElement);

if nLocalDofs(2) == 0
    dofs = [];
    return
end

% dof numbers are simply given by vertex numbers
if isequal(idx, ':')
    dofs = 1:N;
elseif islogical(idx)
    assert(isequal(size(idx), [1, N]))
    dofs = find(idx);
else
    assert(all(0 < idx) && all(idx <= N))
    dofs = idx;
end

offset = obj.mesh.nCoordinates * nLocalDofs(1);
dofs = offset + (nLocalDofs(2)*(dofs-1) + (1:nLocalDofs(2))');

end
