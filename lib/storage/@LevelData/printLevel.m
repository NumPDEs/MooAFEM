function printLevel(obj, jLevel)
%%PRINTLEVEL prints information on the given levels to command line
%   PRINTLEVEL(obj) defaults to final level
%   PRINTLEVEL(obj, jLevel)
%
%   See also LevelData/printHeader

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
        jLevel (1,:) double = obj.nLevel
    end

    % Print header in case of plotting the first level
    if obj.isInitialRun()
        obj.printHeader();
    end

    % Iterate over given list of levels
    specifier = obj.getFormatSpecifier('  ');
    for k = 1:length(jLevel)
        % Extract data of all variables
        data = cell(obj.nVariable, 1);
        ind = 0;
        for jVariable = 1:obj.nVariable
            if obj.isScalar(jVariable)
                ind = ind + 1;
                data{ind} = obj.level(jLevel(k)).(obj.label{jVariable});
            end
        end
        % Print information on current level to command line
        fprintf(specifier, data{1:ind});
    end
end
