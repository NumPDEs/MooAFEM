% computePatchAreas Compute areas of patches in mesh.
%
%   patchAreas = computePatchAreas(mesh)
%
%   See also: Mesh

function patchAreas = computePatchAreas(obj)

trafo = getAffineTransformation(obj);
patchAreas = accumarray(obj.elements(:), repelem(trafo.area, 3), [obj.nCoordinates, 1])';

end
