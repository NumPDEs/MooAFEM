% unitTriangle Construct Mesh object on unit triangle.
%
%   mesh = Mesh.unitTriangle() constructs a mesh with one element and
%       homogeneous dirichlet boundary conditions on the unit triangle. In
%       particular, coordinates, elements, and edges are the same as in the
%       mesh constructed by mesh = Mesh.loadGeometry('unittriangle'). This
%       method, however, is faster.
%
%   See also: Mesh

function obj = unitTriangle()

coordinates = [0, 1, 0; 0, 0, 1];
elements = [2;3;1];
boundaries = {[1, 2, 3; 2, 3, 1]};

obj = Mesh(coordinates, elements, boundaries);

end