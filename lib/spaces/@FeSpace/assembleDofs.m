% assembleDofs Generate global DOF-connectivity arrays from local DOF locations.
%
%   dofs = obj.assembleDofs()

function dofs = assembleDofs(obj)
%% preallocate arrays
mesh = obj.mesh;
dofCoupling = getDofConnectivity(obj.finiteElement);
nDofs = dofCoupling * [3;3;1];
nEdgeDofs = dofCoupling * [2;1;0];
dofs.element2Dofs = zeros(nDofs, mesh.nElements);
dofs.edge2Dofs = zeros(nEdgeDofs, mesh.nEdges);

dofs.nDofs = dofCoupling * [mesh.nCoordinates; mesh.nEdges; mesh.nElements];

%% vertex Dofs
if dofCoupling(1) ~= 0
    dofNumbers = createVertexDofs(obj);
    dofs.element2Dofs(1:3,:) = dofNumbers(mesh.elements);
    dofs.edge2Dofs(1:2,:) = dofNumbers(mesh.edges);
end

%% edge dofs
if dofCoupling(2) ~= 0
    nInnerEdgeDofs = dofCoupling(2);
    dofNumbers = createInnerEdgeDofs(obj);
    idx = dofCoupling(1) * 3;
    for k = 1:3
        idx = idx(end) + (1:nInnerEdgeDofs);
        dofs.element2Dofs(idx,:) = dofNumbers(:,mesh.element2edges(k,:));
        dofs.element2Dofs(idx,mesh.flipEdges(k,:)) = flipud(dofs.element2Dofs(idx,mesh.flipEdges(k,:)));
    end
    idx = dofCoupling(1) * 2 + (1:nInnerEdgeDofs);
    dofs.edge2Dofs(idx,:) = dofNumbers;
end

%% inner dofs
if dofCoupling(3) ~= 0
    idx = dofCoupling(1:2) * [3;3] + (1:dofCoupling(3));
    dofs.element2Dofs(idx,:) = createInnerElementDofs(obj);
end

end
