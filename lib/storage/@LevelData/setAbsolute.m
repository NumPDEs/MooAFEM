function setAbsolute(obj, jLevel, variableName, value)
%%SETABSOLUTE stores specified data of absolute-niveau type (constants / 
%eigenvalues / efficiency indices / ...) to this LevelData object for the
%specified list jLevel of levels, repeated pairs of variable names and values,
%first dimension of data must correspond to number of levels
%   SETABSOLUTE(obj, jLevel, variableName, value, ...)
%   SETABSOLUTE(obj, ':', variableName, value, ...)
%
%   See also LevelData/set

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
        value
    end

    if strcmp(jLevel, ':')
        jLevel = 1:obj.nLevel;
    end

    % Determine number of arguments
    nArgument = length(variableName);

    % Store data to LevelData object
    for j = 1:nArgument
        name = variableName{j};
        val = value{j};
        % Determine type of data
        if isa(val, 'ValueDetails')
            % Copy type argument
            currentType = val;
            valueList = repmat({nan}, length(jLevel), 1);
        else
            % Determine type by example data
            [valueList, currentType] = splitIntoLevelWiseData(length(jLevel), val);
        end

        if ~ismember(name, obj.absoluteVariable)
            % Store type
            obj.type(name) = currentType.adaptPrintWidthTo(name);
            % Store variable as absolute variable
            obj.absoluteVariable{end+1} = name;
        end
        % Store level-oriented data
        for k = 1:length(jLevel)
            obj.level(jLevel(k)).(name) = valueList{k};
        end
    end
end