function plotLevelTriangulation(obj, ax, jLevel)
%%PLOTLEVELTRIANGULATION auxiliary private function for creation of plots
%of triangulations, it extracts the triangulation data in various possible
%formats (compatible with octAFEM and MooAFEM)
%   PLOTLEVELTRIANGULATION(obj, ax, jLevel)

% Copyright 2023 Philipp Bringmann
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%


    % This functions plots a single level only
    assert(numel(jLevel) == 1);

    % Determine names of variables containing mesh information
    arrayDataFound = false;
    transposeArrays = false;
    if isfield(obj.level(jLevel), 'c4n')
        % Arrays in octAFEM style
        coordinatesName = 'c4n';
        elementsName = 'n4e';
        arrayDataFound = true;
    elseif isfield(obj.level(jLevel), 'coordinates')
        % Arrays in MooAFEM style
        coordinatesName = 'coordinates';
        elementsName = 'elements';
        arrayDataFound = true;
        transposeArrays = true;
    end

    % Extract mesh information
    if arrayDataFound
        % Extract data present as coordinates / elements arrays
        c4n = obj.level(jLevel).(coordinatesName);
        n4e = obj.level(jLevel).(elementsName);
        if transposeArrays
            c4n = permute(c4n, [2 1]);
            n4e = permute(n4e, [2 1]);
        end
        % Reorganise mesh information
        c4e = permute(reshape(c4n(n4e,:), length(n4e), 3, 2), ...
                      [3 2 1]);
    else
        % Extract data present as mesh object
        for k = 1:length(obj.label)
            data = obj.level(jLevel).(obj.label(k));
            if isa(data, 'MeshData')
                % octAFEM MeshData object
                c4e = data.c4e;
                break
            elseif isa(data, 'Mesh')
                % MooAFEM Mesh object
                c4n = permute(obj.level(jLevel).(coordinatesName), [2 1]);
                n4e = permute(obj.level(jLevel).(elementsName), [2 1]);
                % Reorganise mesh information
                c4e = permute(reshape(c4n(n4e,:), length(n4e), 3, 2), ...
                              [3 2 1]);
            end
        end
    end

    % Plot mesh
    edgecolor = [0 0.4470 0.7410];
    patch(ax, permute(c4e(1,:,:), [2 3 1]), permute(c4e(2,:,:), [2 3 1]), ...
          'white', 'EdgeColor', edgecolor);

    % Format plot
    title(ax, {'Mesh plot'; [num2str(size(c4e, 3)), ' triangles']});
    axis(ax, 'equal');
    axis(ax, 'tight');
end