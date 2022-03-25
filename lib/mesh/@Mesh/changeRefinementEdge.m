% changeRefinementEdge Rotate internal data structures such that given
%   edge is new refinement edge per element.
%
% mesh.changeRefinementEdge(newRefinementEdge) changes the elements,
%   element2edges, and flipEdges data structures of the mesh such that the
%   new edge is the edge given in newRefinementEdge per element. New edges
%   are given as index in {1,2,3}, according to the old element2edges data
%   structure.

function changeRefinementEdge(obj, newRefinementEdge)

arguments
    obj
    newRefinementEdge (1,:) double {mustBeInRange(newRefinementEdge, 1, 3)}
end

for k = 2:3
    idx = find(newRefinementEdge == k);
    obj.element2edges(:,idx) = circshift(obj.element2edges(:,idx), 4-k, 1);
    obj.elements(:,idx) = circshift(obj.elements(:,idx), 4-k, 1);
    obj.flipEdges(:,idx) = circshift(obj.flipEdges(:,idx), 4-k, 1);
end

end