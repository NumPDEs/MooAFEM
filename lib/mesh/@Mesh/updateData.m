% updateData Set data of Mesh object
%
%   updateData(obj) updates the data stored in the calling mesh based on
%       refinement groups given by a refinement routine.
%
%   See also Mesh

function updateData(obj, bisecData)
%% create indices of ...
% new nodes
edgeNodes = obj.nCoordinates + (1:nnz(bisecData.bisectedEdges))';
innerNodes = cumsum([edgeNodes(end); bisecData.nRefinedElements.*bisecData.nInnerNodes]);

% new edges
innerEdges = cumsum([obj.nEdges + nnz(bisecData.bisectedEdges); bisecData.nRefinedElements.*bisecData.nInnerEdges]);
nChildEdges = 1 + bisecData.bisectedEdges;
oldEdge2parent = cumsum([1; nChildEdges]);

% new elements
nChildElems = ones(obj.nElements, 1);
for k = 1:bisecData.nBisecGroups
    nChildElems(bisecData.refinedElements(k)) = bisecData.nDescendants(k);
end
oldElement2parent = cumsum([1; nChildElems]);

%% create new coordinates
newCoordinates = [obj.coordinates, ...
    edgewiseCoordinates(obj, [1;1]/2, bisecData.bisectedEdges), ...
    zeros(2, innerNodes(end)-innerNodes(1))];

for k = find(bisecData.nInnerNodes ~= 0)'
    newCoordinates(:,innerNodes(k)+1:innerNodes(k+1)) = ...
        elementwiseCoordinates(obj, bisecData.bisecGroups{k}.innerNodes, bisecData.refinedElements(k));
end

%% refine existing edges
newEdges = [moveParents(obj.edges, oldEdge2parent), ...
    zeros(2, innerEdges(end) - innerEdges(1))];

children = getChildIndices(oldEdge2parent(bisecData.bisectedEdges), 2);
newEdges(:,children) = [obj.edges(1,bisecData.bisectedEdges), edgeNodes'; ...
    edgeNodes', obj.edges(2,bisecData.bisectedEdges)];

%% refine boundaries
newBoundaries = cell(obj.nBoundaries, 1);
for k = 1:obj.nBoundaries
    newBoundaries{k} = oldEdge2parent(obj.boundaries{k});
    bndChildren = newBoundaries{k}(bisecData.bisectedEdges(obj.boundaries{k})) + 1;
    newBoundaries{k} = [newBoundaries{k}; bndChildren];
end

%% allocate new element-wise arrays
oldEdgeNewIdx = oldEdge2parent(obj.element2edges);
newElements = moveParents(obj.elements, oldElement2parent);
newE2E = moveParents(oldEdgeNewIdx, oldElement2parent);
newFlip = moveParents(obj.flipEdges, oldElement2parent);

%% generate bisection rule dependent data
[left, right] = getLeftAndRightEdgePerElement(obj, bisecData.bisectedEdges, oldEdgeNewIdx);
edgeMidPts = zeros(obj.nEdges,1);
edgeMidPts(bisecData.bisectedEdges) = edgeNodes;
edgeMidPts = edgeMidPts(obj.element2edges);

for k = 1:bisecData.nBisecGroups
    idx = (innerEdges(k)+1):innerEdges(k+1);
    intEdges = reshape(idx, [], bisecData.nInnerEdges(k))';
    intNodes = reshape((innerNodes(k)+1:innerNodes(k+1))', bisecData.nInnerNodes(k), []);
    
    parentEdges = oldEdgeNewIdx(:, bisecData.refinedElements(k));
    parentElements = oldElement2parent(bisecData.refinedElements(k));
    children = getChildIndices(parentElements, bisecData.nDescendants(k));
    
    arg = restrictTo(bisecData.refinedElements(k), ...
        obj.elements, edgeMidPts, left, right, obj.flipEdges);
    
    newEdges(:,idx) = bisecData.bisecGroups{k}.createNewEdges(arg{1}, arg{2}, intNodes);
    newElements(:,children) = bisecData.bisecGroups{k}.refineElement(arg{1}, arg{2}, intNodes);
    newE2E(:,children) = bisecData.bisecGroups{k}.refineEdgeConnectivity(parentEdges, arg{3}, arg{4}, intEdges);
    newFlip(:,children) = bisecData.bisecGroups{k}.refineEdgeOrientation(arg{5});
end

%% assign new arrays
obj.coordinates = newCoordinates;
obj.elements = newElements;
obj.edges = newEdges;
obj.element2edges = newE2E;
obj.flipEdges = newFlip;
obj.boundaries = newBoundaries;

end

%% local functions
function children = getChildIndices(parents, nDescendants)
    children = reshape((0:(nDescendants-1)) + parents, [], 1);
end

function newArray = moveParents(oldArray, parentPositions)
    newArray = zeros(size(oldArray, Dim.Vector), parentPositions(end)-1, class(oldArray));
    newArray(:,parentPositions(1:(end-1))) = oldArray;
end

function restrictedArrays = restrictTo(idx, varargin)
    restrictedArrays = cellfun(@(x) x(:,idx), varargin, 'UniformOutput', false);
end

function [left, right] = getLeftAndRightEdgePerElement(mesh, markedEdges, oldEdgeNewIdx)
    % here, we use the fact that left and right edge index only differ by 1
    bisectedEdges = markedEdges(mesh.element2edges);
    left = oldEdgeNewIdx.*bisectedEdges;
    right = left + bisectedEdges;
    left = left + mesh.flipEdges;
    right = right - mesh.flipEdges;
end
