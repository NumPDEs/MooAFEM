function data = get(obj, jLevel, variableName)
%%GET extracts data from this LevelData object on a given list jLevel of 
%level numbers. One or more variables to be returned can be specified by their
%names.
%   data = GET(obj, jLevel, variableName, ...)
%   data = GET(obj, ':', variableName, ...)

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
        jLevel {mustBeIndexVector} = ':'
    end

    arguments (Repeating)
        variableName {mustBeTextScalar}
    end

    if strcmp(jLevel, ':')
        jLevel = 1:obj.nLevel;
    end

    % Initialise return variable
    data = nan(length(jLevel), length(variableName));
    containsCharOnly = true;

    % filter non-existing variables
    existingVariables = ismember(variableName, obj.label);
    if any(~existingVariables)
        warning(['Variable(s) ', strjoin(variableName(~existingVariables), ', '), ' not found']);
    end
    variableName = variableName(existingVariables);

    % Iterate over all variables
    for jVariable = 1:length(variableName)
        name = variableName{jVariable};

        % Determine index of current variable
        idx = obj.getIndex(name);
        if obj.isScalar(idx)
            % Extract scalar variables for each level
            value = zeros(length(jLevel), 1);
            for k = 1:length(jLevel)
                value(k) = obj.level(jLevel(k)).(name);
            end
        else
            % Extract non-scalar variables only in case of a single level
            if length(jLevel) > 1
                warning('Export of level-oriented scalar variables only');
                value = [];
            else
                data(1,jVariable) = obj.level(jLevel).(name);
                return
            end
        end

        % Save extracted data to return variable
        data(:,jVariable) = value;
        % Post-process character arrays
        if obj.type.(name).rawType == RawType.TEXT
            data = mat2cell(char(data), ones(length(data), 1));
        end
    end
end
