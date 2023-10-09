% createInnerElementDofs Generate global DOF-numbers associated with inner
%   element dofs.
%
%   dofs = createInnerElementDofs(obj) generates dofs for all elements.
%
%   dofs = createInnerElementDofs(obj, idx) generates dofs only for elements
%       given by index vector idx.

function dofs = createInnerElementDofs(obj, idx)

arguments
    obj
    idx (1,:) {mustBeIndexVector} = ':'
end

N = obj.mesh.nElements;
nLocalDofs = getDofConnectivity(obj.finiteElement);

if nLocalDofs(3) == 0
    dofs = [];
    return
end

% dof numbers are simply given by vertex numbers
if ischar(idx) && isequal(idx, ':')
    dofs = 1:N;
elseif islogical(idx)
    assert(isequal(size(idx), [1, N]))
    dofs = find(idx);
else
    assert(all(0 < idx) && all(idx <= N))
    dofs = idx;
end

offset = nLocalDofs(1:2) * [obj.mesh.nCoordinates; obj.mesh.nEdges];
dofs = offset + (nLocalDofs(3)*(dofs-1) + (1:nLocalDofs(3))');

end
