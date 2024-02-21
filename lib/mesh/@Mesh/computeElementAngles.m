% computeElementAngles Compute the three interior angles of all elements in mesh.
% 
%   diam = computeElementAngles(mesh)
%
%   See also: Mesh

function angles = computeElementAngles(obj)

trafo = getAffineTransformation(obj);

% Compute angles opposite of each edge by law of cosines
lengths = trafo.ds(obj.element2edges);
angles = acos((lengths([2 3 1],:).^2 + lengths([3 1 2],:).^2 - lengths([1 2 3],:).^2) ...
         ./ (2*lengths([2 3 1],:).*lengths([3 1 2],:)));

% Rearrange angles according to vertices of elements
angles = angles([2 3 1],:);

end
