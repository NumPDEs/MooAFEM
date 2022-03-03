% assembleFreeDofs Create array of Dofs not on dirichlet boundary.
%
%   freeDofs = obj.assembleFreeDofs()

function freeDofs = assembleFreeDofs(obj)

dofs = getDofs(obj);
mesh = obj.mesh;

allDofs = true(dofs.nDofs, 1);
dirichletEdges = getCombinedBndEdges(mesh, obj.dirichlet);
allDofs(dofs.edge2Dofs(:,dirichletEdges)) = false;
freeDofs = find(allDofs);

end
