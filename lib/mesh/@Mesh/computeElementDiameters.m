% computeElementDiameters Compute diameters of all elements in mesh.
% 
%   diam = computeElementDiameters(mesh)
%
%   See also: Mesh

function diam = computeElementDiameters(obj)

trafo = getAffineTransformation(obj);
diam = max(trafo.ds(obj.element2edges), [], 1);

end
