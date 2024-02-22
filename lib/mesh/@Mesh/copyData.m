% copyData Copy data from other Mesh object
%
%   copyData(obj, otherMesh) updates the data stored in the calling mesh based
%       on other mesh.
%
%   See also Mesh

function copyData(obj, otherMesh)

assert(isa(otherMesh, 'Mesh'))

% if it is a mesh, data must be valid -> just copy
obj.coordinates = otherMesh.coordinates;
obj.elements = otherMesh.elements;
obj.edges = otherMesh.edges;
obj.element2edges = otherMesh.element2edges;
obj.edge2elements = otherMesh.edge2elements;
obj.flipEdges = otherMesh.flipEdges;
obj.boundaries = otherMesh.boundaries;

obj.notify('CompletelyChanged')

end
