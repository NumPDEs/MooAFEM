% clone Create new mesh with same data
%
%   newMesh = clone(obj)
%
%   See also Mesh

function newMesh = clone(obj)

newMesh = Mesh([], [], []);
copyData(newMesh, obj);

end
