% loadFromGeometry Construct Mesh object from geometry.
%
%   mesh = Mesh.loadFromGeometry(geometryIdentifier) loads geometry given by the
%       string geometryIdentifier. The geometry must be given in a subfolder of
%       @Mesh/geometries containing coordinates.dat, elements.dat and
%       boundary<N>.dat files, where N are successive integers beginning at 1.
%
%   See also: Mesh

function obj = loadFromGeometry(geometryIdentifier)

arguments
    geometryIdentifier (1,1) string
end

%% access folder
[meshDirectory, ~, ~] = fileparts(mfilename('fullpath'));
geometryDirectory = fullfile(meshDirectory, 'geometries', geometryIdentifier);
assert(isfolder(geometryDirectory), 'Can''t find directory %s.', geometryDirectory)

%% load coordinates / elements
coordinatesFile = fullfile(geometryDirectory, 'coordinates.dat');
assert(isfile(coordinatesFile), 'Can''t find file ''coordinates.dat'' in %s.', geometryDirectory)
elementsFile = fullfile(geometryDirectory, 'elements.dat');
assert(isfile(elementsFile), 'Can''t find file ''elements.dat'' in %s.', geometryDirectory)
coordinates = load(coordinatesFile)';
elements = load(elementsFile)';

%% load boundary parts
bndFiles = dir(fullfile(geometryDirectory, 'boundary*.dat'));
assert(~isempty(bndFiles), 'No boundary parts stored in %s.', geometryDirectory)
boundaries = cell(numel(bndFiles),1);
for i = 1:numel(bndFiles)
    % bndFiles is sorted alphabetically -> boundaryN.dat is assigned to boundaries{N}
    file = fullfile(geometryDirectory, bndFiles(i).name);
    boundaries{i} = load(file)';
end

obj = Mesh(coordinates, elements, boundaries);

end