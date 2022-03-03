% saveTikzConforming Save mesh data for later use in TikZ (with LaTeX).
%
%   Save mesh data in a file suitable for \input{} within a tikzpicture
%   environment in LaTeX.
%   The drawing style is defined in the first line of the output file for
%   convenient manual editing.
%
%   obj.saveTikzConforming() saves meshdata in file 'mesh.tex'.
%
%   obj.saveTikzConforming(fileName) saves meshdata in file with custom name.
%
%   See also Mesh

function saveTikzConforming(obj, fileName)

arguments
    obj (1,1) Mesh
    fileName (1,1) string = 'mesh.tex'
end

%% open file
fileID = fopen(fileName, 'a');
assert(fileID ~= -1, 'Cannot open file %s.', fileName)
fprintf(fileID, '\\tikzset{meshStyle/.style={very thin, black, fill=none}}\n\n');

%% define coordinates in tikz
fprintf(fileID, '\\coordinate (P%d) at (%.4g,%.4g);\n', [1:obj.nCoordinates; obj.coordinates]);
fprintf(fileID, '\n');
    
%% draw mesh lines
fprintf(fileID, '\\draw[meshStyle] (P%d) -- (P%d) -- (P%d) -- cycle;\n', obj.elements);
fclose(fileID);

end
