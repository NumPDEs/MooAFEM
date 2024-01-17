% updateData Set data of Mesh object
%
%   updateData(obj) updates the data stored in the calling mesh based on
%       refinement groups given by a refinement routine.
%
%   See also Mesh

function updateData(obj, bisecData)
%% create old-to-new connectivity
oldEdge2parent = cumsum([1; 1+bisecData.bisectedEdges]);
oldElement2parent = cumsum([1; getNChildrenPerElement(bisecData)]);

%% create new coordinates
% new nodes = old nodes | nodes on edges | interior nodes
idxN = cumsum([1; obj.nCoordinates; nnz(bisecData.bisectedEdges); ...
    bisecData.nRefinedElements.*bisecData.nInnerNodes]);

newCoordinates = zeros(2, idxN(end)-1);
newCoordinates(:,idxN(1):idxN(2)-1) = obj.coordinates;
newCoordinates(:,idxN(2):idxN(3)-1) = edgewiseCoordinates(obj, [1;1]/2, bisecData.bisectedEdges);

for k = find(bisecData.nInnerNodes ~= 0)'
    newCoordinates(:,idxN(k+2):idxN(k+3)-1) = ...
        elementwiseCoordinates(obj, bisecData.bisection{k}.innerNodes, getRefinedElementIdx(bisecData, k));
end

%% refine existing edges
% new edges = old edges + children | interior edges
idxS = cumsum([1; obj.nEdges + nnz(bisecData.bisectedEdges); ...
    bisecData.nRefinedElements.*bisecData.nInnerEdges]);

newEdges = zeros(2, idxS(end)-1);
newEdges(:,oldEdge2parent(1:(end-1))) = obj.edges;

edgeNodes = (idxN(2):idxN(3)-1)';
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

%% refine element data
newElements = zeros(3, oldElement2parent(end)-1);
newEl2Ed = zeros(3, oldElement2parent(end)-1);
newFlip = false(3, oldElement2parent(end)-1);

oldEdgeNewIdx = oldEdge2parent(obj.element2edges);
elems = bisecData.getRefinedElementIdx(0);
newElements(:,oldElement2parent(elems)) = obj.elements(:,elems);
newEl2Ed(:,oldElement2parent(elems)) = oldEdgeNewIdx(:,elems);
newFlip(:,oldElement2parent(elems)) = obj.flipEdges(:,elems);

%% get left and right edge per element
% here, we use the fact that left and right edge index only differ by 1
bisectedEdges = bisecData.bisectedEdges(obj.element2edges);
left = oldEdgeNewIdx.*bisectedEdges;
right = left + bisectedEdges;
left = left + obj.flipEdges;
right = right - obj.flipEdges;

%% refine elements + add inner edges
edgeMidPts = zeros(obj.nEdges,1);
edgeMidPts(bisecData.bisectedEdges) = edgeNodes;
edgeMidPts = edgeMidPts(obj.element2edges);

for k = find(bisecData.nRefinedElements)'
    idx = idxS(k+1):idxS(k+2)-1;
    elems = getRefinedElementIdx(bisecData, k);
    intEdges = reshape(idx, [], bisecData.nInnerEdges(k))';
    intNodes = reshape((idxN(k+2):idxN(k+3)-1)', bisecData.nInnerNodes(k), []);
    
    parentEdges = oldEdgeNewIdx(:, elems);
    parentElements = oldElement2parent(elems);
    children = getChildIndices(parentElements, bisecData.nDescendants(k));
    
    arg = restrictTo(elems, obj.elements, edgeMidPts, left, right, obj.flipEdges);
    
    newEdges(:,idx) = bisecData.bisection{k}.createNewEdges(arg{1}, arg{2}, intNodes);
    newElements(:,children) = bisecData.bisection{k}.refineElement(arg{1}, arg{2}, intNodes);
    newEl2Ed(:,children) = bisecData.bisection{k}.refineEdgeConnectivity(parentEdges, arg{3}, arg{4}, intEdges);
    newFlip(:,children) = bisecData.bisection{k}.refineEdgeOrientation(arg{5});
end

%% recompute adjacency information on the edges
newEd2El = Mesh.computeEdgeAdjacency(newEl2Ed, newFlip);

%% assign new arrays
obj.coordinates = newCoordinates;
obj.elements = newElements;
obj.edges = newEdges;
obj.element2edges = newEl2Ed;
obj.edge2elements = newEd2El;
obj.flipEdges = newFlip;
obj.boundaries = newBoundaries;

end

%% local functions
function children = getChildIndices(parents, nDescendants)
    children = reshape((0:(nDescendants-1)) + parents, [], 1);
end

function restrictedArrays = restrictTo(idx, varargin)
    restrictedArrays = cellfun(@(x) x(:,idx), varargin, 'UniformOutput', false);
end
