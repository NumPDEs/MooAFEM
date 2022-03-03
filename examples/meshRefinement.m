% ******************************************************************************
% Example that showcases how meshes can be refined.
% ******************************************************************************

%% local refinement
% A mesh can be refined locally, i.e., only certain marked elements are refined.
% The marked elements are given by their index in the mesh.elements array. After
% refining the marked elements, the mesh closure is built to ensure a conforming
% mesh without hanging nodes.
geo = 'Lshape';
mesh = Mesh.loadFromGeometry(geo);
figure()
plot(mesh, 'labelElements', true, 'labelEdges', true)
title('Reference Mesh')
marked = 1:min(3, mesh.nElements);
mesh.refineLocally(marked);
figure()
plot(mesh, 'labelElements', true)
title('Element based refinement')

% The default method is Newest Vertex Bisection, where all 3 edges in a triangle
% are bisected. Other element based methods are available as additional option:
% - NVB{1,3,5}: Element based Newest Vertex Bisecton with 1, 3, or 5 bisected
%       edges per step.
% - RGB: Element based Red-Green-Blue refinement.
% e.g., mesh.refineLocally(marked, 'NVB1');


% Instead of elements, also edges can be marked for refinement with the method
% 'NVBEdge'
mesh = Mesh.loadFromGeometry(geo);
marked = 1:min(3, mesh.nEdges);
mesh.refineLocally(marked, 'NVBEdge');
figure()
plot(mesh, 'labelEdges', true)
title('Edge based refinement')

%% uniform refinement
% A mesh can also be refined uniformly, i.e., all elements/edges are refined n 
% times given a refinement method.
mesh = Mesh.loadFromGeometry(geo);
mesh.refineUniform(3, 'RGB');
figure()
plot(mesh)
title('3 times uniform RGB-refinement')
