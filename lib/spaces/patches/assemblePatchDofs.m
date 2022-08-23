function dofs = assemblePatchDofs(obj, ~)

%% preallocate arrays
mesh = obj.mesh;
dofs = getDofs(obj);
patchDofs = cell(mesh.nCoordinates, 1);

free = true(mesh.nEdges,1);
free(getCombinedBndEdges(mesh, obj.bnd.dirichlet)) = false;
free = find(free)';
[patchEdges, pEdges] = computePatchMembership(mesh.edges(:,free));
patchEdges = free(patchEdges);

[patchElements, pElements] = computePatchMembership(mesh.elements);

vertexLiesOnDirichletBnd = ...
    ismember(1:mesh.nCoordinates, mesh.edges(:,getCombinedBndEdges(mesh, obj.bnd.dirichlet)));

% TODO: get rid of this loop by means of mat2cell or similar
for k = 1:mesh.nCoordinates
    if vertexLiesOnDirichletBnd(k)
        vertexDof = [];
    else
        vertexDof = asVector(createVertexDofs(obj, k));
    end
    
    edges = patchEdges(pEdges(k):(pEdges(k+1)-1));
    elements = patchElements(pElements(k):(pElements(k+1)-1));
    patchDofs{k} = [vertexDof; ...
        asVector(createInnerEdgeDofs(obj, edges));...
        asVector(createInnerElementDofs(obj, elements))];
end

dofs.patch2Dofs = patchDofs;

end

function [patches, pointer] = computePatchMembership(objects)
    n = size(objects, 1);
    N = size(objects, 2);
    [sortedNodes, idx] = sort(objects(:));
    patches = repmat(1:N, [n 1]);
    patches = patches(idx);
    pointer = [0; find(diff(sortedNodes)); n*N] + 1;
end