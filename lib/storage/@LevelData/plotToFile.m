function plotToFile(obj, xVariable, yVariable)
%%PLOTTOFILE creates plots of various scalar variables (specified in 
%yVariable) with respect %to the variable xVariable and stores it to a file,
%filename is generated automatically from information in LevelData object
%   PLOTTOFILE(obj, xVariable, yVariable)
%
%   See also LevelData/plot

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

    arguments
        obj
        xVariable {mustBeTextScalar}
    end

    arguments (Repeating)
        yVariable {mustBeTextScalar}
    end

    if isempty(yVariable)
        yVariable = setdiff(obj.label, xVariable);
    end

    ensureFolderExists(obj.foldername);
    h = createStandardFigure();
    obj.plot(xVariable, yVariable{:});

    % Export plot
    print(h, '-dpng', '-r600', obj.foldername + '/' + obj.filename + '.png');
    close(h);
end
