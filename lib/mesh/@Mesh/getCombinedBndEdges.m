% getCombinedBndEdges Combined edges of given boundary parts.
%
%   edges = getCombinedBndEdges(mesh, idx) returns the combined edges of all
%       boundary parts given by idx.
%
%   See also: Mesh

function edges = getCombinedBndEdges(obj, idx)

arguments
    obj
    idx (:,1) {mustBeIndexVector} = ':'
end

assert(isequal(idx, ':') | all(0 < idx & idx <= obj.nBoundaries), ...
    'Index vector must contain valid indices of boundary parts.')

allEdges = false(obj.nEdges, 1);
allEdges(vertcat(obj.boundaries{idx})) = true;
edges = find(allEdges);

end
