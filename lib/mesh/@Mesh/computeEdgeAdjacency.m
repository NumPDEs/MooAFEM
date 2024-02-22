% computeEdgeAdjacency Compute information on elements neighbouring an edge
%   respecting the edge orientation; second entry is zero in case of a
%   boundary edge
%
%   edge2elements = computeEdgeAdjacency(element2edges, flipEdges)
%       computes neighboring information from connectivity and orientation
%       arrays.
%
%   See also: Mesh

function edge2elements = computeEdgeAdjacency(element2edges, flipEdges)

arguments
    element2edges (3,:) double
    flipEdges (3,:)
end

% abbreviate numbering of elements
nElems = size(element2edges, Dim.Elements);
elemNumbers = repmat((1:nElems), 3, 1);

% determine the two elements adjacent to an edge
pickFirst = @(x) x(1);
neighbor1 = accumarray(element2edges(:), elemNumbers(:), [], pickFirst);
neighbor2 = accumarray(element2edges(:), elemNumbers(:)) - neighbor1;

% check whether an orientation needs to be flipped
isExchanged = accumarray(element2edges(:), flipEdges(:), [], pickFirst);
% isInterior = accumarray(element2edges(:), flipEdges(:), [], @max);

% assemble edge2elements array and correct depending on the orientation of
% the edge
edge2elements = transpose([neighbor1, neighbor2]);
edge2elements(isExchanged) = flipud(edge2elements(isExchanged));

end