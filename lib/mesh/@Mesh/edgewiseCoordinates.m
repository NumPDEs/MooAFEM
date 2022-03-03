% edgewiseCoordinates Compute edge-wise coordinates from barycentric
%   coordinates.
%
%   Given barycentric coordinates (b1,b2) and edges with nodes (z1,z2), compute
%   the cartesian coordinate z1*b1+z2*b2 for certain edges. If more than one
%   barycentric coordinate is given, their respective coordinates are stored
%   along the 3rd dimension of the resulting array.
%
%   points = edgewiseCoordinates(mesh, bary) takes all edges.
%
%   points = edgewiseCoordinates(mesh, bary, idx) only takes edges given by idx.
%
%   See also: Mesh

function points = edgewiseCoordinates(obj, bary, idx)

arguments
    obj
    bary (1,1) Barycentric1D
    idx (:,1) {mustBeIndexVector} = ':'
end

%% compute cartesian coordinates
points = vectorProduct(...
    reshape(obj.coordinates(:,obj.edges(:,idx)), 4, []), ...
    bary.arangeAlongDim(3), ...
    [2,2], [2,1]);

end
