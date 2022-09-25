% assemblePatchDofs Assembles lists of unique inner dofs for each patch in
%   underlying mesh.
%
%   assemblePatchDofs(obj)

function [patchDofs, patchElements] = assemblePatchDofs(obj)
%% preallocate arrays
fes = obj.fes;
mesh = fes.mesh;

%% compute patch membership of (free) vertices/edges/elements
dirichletEdges = getCombinedBndEdges(mesh, fes.bnd.dirichlet);
vertexNotOnDirichletBnd = ~ismember(1:mesh.nCoordinates, mesh.edges(:,dirichletEdges));

free = computeFreeEdges(mesh, dirichletEdges);
[patchEdges, nEdgesPerPatch] = computePatchMembership(mesh.edges(:,free));
patchEdges = free(patchEdges);

[patchElements, nElementsPerPatch] = computePatchMembership(mesh.elements);

%% combine vertex/edge/element dofs for each patch
vertexDofs = asVector(createVertexDofs(fes, vertexNotOnDirichletBnd));
edgeDofs = asVector(createInnerEdgeDofs(fes, patchEdges));
elementDofs = asVector(createInnerElementDofs(fes, patchElements));

nLocalDofs = getDofConnectivity(fes.finiteElement);
[patchDofs, nDofsPerPatch] = interleaveGroupedArrays(...
    vertexDofs, nLocalDofs(1)*double(vertexNotOnDirichletBnd), ...
    edgeDofs, nLocalDofs(2)*nEdgesPerPatch', ...
    elementDofs, nLocalDofs(3)*nElementsPerPatch');

patchDofs = mat2cell(patchDofs', nDofsPerPatch);
patchElements = mat2cell(patchElements, nElementsPerPatch);

end

%% auxiliary functions

% compute free edges without setdiff -> linear complexity
function freeEdges = computeFreeEdges(mesh, dirichlet)
    freeEdges = true(mesh.nEdges,1);
    freeEdges(dirichlet) = false;
    freeEdges = find(freeEdges)';
end

% requires one global sort of connectivity arrays -> log-linear complexity
function [patches, count] = computePatchMembership(objects)
    n = size(objects, 1);
    N = size(objects, 2);

    count = accumarray(objects(:), 1);
    ptr = cumsum([1; count]);
    patches = zeros(n*N, 1);

    for i = 1:N
        for k = 1:n
            idx = objects(k,i);
            patches(ptr(idx)) = i;
            ptr(idx) = ptr(idx) + 1;
        end
    end
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
