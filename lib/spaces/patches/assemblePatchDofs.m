% assemblePatchDofs Assembles lists of unique inner dofs for each patch in
%   underlying mesh.
%
%   assemblePatchDofs(obj)

function patchDofs = assemblePatchDofs(obj)
%% preallocate arrays
mesh = obj.mesh;

%% compute patch membership of (free) vertices/edges/elements
dirichletEdges = getCombinedBndEdges(mesh, obj.bnd.dirichlet);
vertexNotOnDirichletBnd = ~ismember(1:mesh.nCoordinates, mesh.edges(:,dirichletEdges));

free = computeFreeEdges(mesh, dirichletEdges);
[patchEdges, pEdges] = computePatchMembership(mesh.edges(:,free));
patchEdges = free(patchEdges);

[patchElements, pElements] = computePatchMembership(mesh.elements);

%% combine vertex/edge/element dofs for each patch
vertexDofs = asVector(createVertexDofs(obj, vertexNotOnDirichletBnd));
edgeDofs = asVector(createInnerEdgeDofs(obj, patchEdges));
elementDofs = asVector(createInnerElementDofs(obj, patchElements));

nLocalDofs = getDofConnectivity(obj.finiteElement);
[patchDofs, nDofsPerPatch] = interleaveGroupedArrays(...
    vertexDofs, nLocalDofs(1)*double(vertexNotOnDirichletBnd), ...
    edgeDofs, nLocalDofs(2)*diff(pEdges)', ...
    elementDofs, nLocalDofs(3)*diff(pElements)');

patchDofs = mat2cell(patchDofs', nDofsPerPatch);

end

%% auxiliary functions

% compute free edges without setdiff -> linear complexity
function freeEdges = computeFreeEdges(mesh, dirichlet)
    freeEdges = true(mesh.nEdges,1);
    freeEdges(dirichlet) = false;
    freeEdges = find(freeEdges)';
end

% requires one global sort of connectivity arrays -> log-linear complexity
function [patches, pointer] = computePatchMembership(objects)
    n = size(objects, 1);
    N = size(objects, 2);
    [sortedNodes, idx] = sort(objects(:));
    patches = repmat(1:N, [n 1]);
    patches = patches(idx);
    pointer = [0; find(diff(sortedNodes)); n*N] + 1;
end

function [totalArray, totalGroupSize] = interleaveGroupedArrays(array, groupSize)
arguments (Repeating)
    array
    groupSize
end

totalGroupPointer = cumsum([1; reshape(vertcat(groupSize{:}), [], 1)])';
totalArray = zeros(1, totalGroupPointer(end)-1);
totalGroupSize = sum(vertcat(groupSize{:}), 1);

N = numel(groupSize);
for n = 1:N
    % add vector [0 1 2 0 1 0 1 2 3 ...] to group starting locations to
    % make space for all group members
    offset = totalGroupPointer(n:N:(end-1)) - cumsum([0, groupSize{n}(1:(end-1))]);
    idx = (1:numel(array{n})) + repelem(offset, groupSize{n}) - 1;
    totalArray(idx) = array{n};
end

end