% ******************************************************************************
% Example that showcases how meshes can be generated, modified and plotted.
% ******************************************************************************

%% generation
% In general, it is best to store meshes as coordinate, element and boundary
% files in a subfolder of lib/@Mesh/geometries. Then, the Mesh can be
% constructed easily, given the name of the subfolder.
mesh = Mesh.loadFromGeometry('Lshape');

% A mesh can also be cloned.
referenceMesh = clone(mesh);

%% modification
% The data of a mesh can be transformed (translate, scale, rotate, skew) easily.
% However, due to performance reasons, the fields of a mesh cannot be altered
% directly. Instead, a copy must be made.
newCoordinates = scaleCoordinates(mesh.coordinates, 2, [0;0]);
newCoordinates = rotateCoordinates(newCoordinates, 3*pi/4);
newCoordinates(:,5) = [0;0.5];
boundaries = cellfun(@(bnd) mesh.edges(:,bnd), mesh.boundaries, 'UniformOutput', false);
mesh = Mesh(newCoordinates, mesh.elements, boundaries);

%% plotting
% A mesh can be plotted via the plot method. It accepts all name-value options
% from the built-in patch plotting command.
figure()
plot(mesh, 'lineWidth', 2, 'FaceColor', 'b');

% Also, enumerations of coordinates and elements can be displayed, as well as a
% coloring of the different boundary parts. All combinations of patch- and
% custom-options are valid.
figure()
plot(referenceMesh, 'labelCoordinates', true, ...
                    'labelElements', true, ...
                    'labelEdges', true, ...
                    'colorBoundaries', true);

% Meshes can also be exported in a format suitable as \input into a TikZpicture
% environment within LaTeX.
mesh.saveTikzConforming()
delete('./mesh.tex')