% elementwiseCoordinates Compute element-wise coordinates from barycentric
%   coordinates.
%
%   Given barycentric coordinates (b1,b2,b3) and elements with nodes (z1,z2,z3),
%   compute the cartesian coordinate z1*b1+z2*b2+z3*b3 on certain elements. If
%   more than one barycentric coordinate is given, their respective coordinates
%   are stored along the 3rd dimension of the resulting array.
%
%   points = elementwiseCoordinates(mesh, bary) takes all elements.
%
%   points = elementwiseCoordinates(mesh, bary, idx) Take only elements
%       given by idx.
%
%   See also: Mesh

function points = elementwiseCoordinates(obj, bary, idx)

arguments
    obj
    bary (1,1) Barycentric2D
    idx (:,1) {mustBeIndexVector} = ':'
end

%% compute cartesian coordinates
points = vectorProduct(...
    reshape(obj.coordinates(:,obj.elements(:,idx)), 6, []), ...
    bary.arangeAlongDim(3), ...
    [2,3], [3,1]);

end
