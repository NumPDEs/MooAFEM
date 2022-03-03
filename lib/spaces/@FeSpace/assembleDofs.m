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
n = 0; % current number of local ...
k = 0; % ... local edge ...
N = 0; % ... and global Dofs

%% vertex Dofs
if dofCoupling(1) ~= 0
    dofs.element2Dofs(1:3,:) = mesh.elements;
    dofs.edge2Dofs(1:2,:) = mesh.edges;
    k = k + 2;
    n = n + 3;
    N = N + mesh.nCoordinates;
end

%% edge dofs
if dofCoupling(2) ~= 0
    nInnerEdgeDofs = dofCoupling(2);
    innerEdgeDofs = N + reshape(1:(nInnerEdgeDofs*mesh.nEdges), nInnerEdgeDofs, mesh.nEdges);
    for i = 1:3
        idx = (n+1):(n+nInnerEdgeDofs);
        dofs.element2Dofs(idx,:) = innerEdgeDofs(:,mesh.element2edges(i,:));
        dofs.element2Dofs(idx,mesh.flipEdges(i,:)) = flipud(dofs.element2Dofs(idx,mesh.flipEdges(i,:)));
        n = n + nInnerEdgeDofs;
    end
    dofs.edge2Dofs((k+1):(k+nInnerEdgeDofs),:) = innerEdgeDofs;
    N = N + nInnerEdgeDofs*mesh.nEdges;
end

%% inner dofs
if dofCoupling(3) ~= 0
    nInnerDofsPerElement = dofCoupling(3);
    dofs.element2Dofs((n+1):(n+nInnerDofsPerElement),:) ...
        = N + reshape(1:(nInnerDofsPerElement*mesh.nElements), nInnerDofsPerElement, mesh.nElements);
    N = N + nInnerDofsPerElement*mesh.nElements;
end

dofs.nDofs = N;

end
