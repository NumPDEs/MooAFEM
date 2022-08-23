% createVertexDofs Generate global DOF-numbers associated with vertex dofs.
%
%   dofs = createVertexDofs(obj) generates dofs for all vertices.
%
%   dofs = createVertexDofs(obj, idx) generates dofs only for vertices given by
%       index vector idx.

function dofs = createVertexDofs(obj, idx)

arguments
    obj
    idx (1,:) {mustBeIndexVector} = ':'
end

N = obj.mesh.nCoordinates;
nLocalDofs = getDofConnectivity(obj.finiteElement);

if nLocalDofs(1) == 0
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

end
