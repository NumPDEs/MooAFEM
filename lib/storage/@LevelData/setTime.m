function setTime(obj, jLevel, variableName, value)
%%SETTIME stores specified data of increasing type (time / ndof / ...) to 
%this LevelData object for the specified list jLevel of levels, repeated pairs
%of variable names and values, first dimension of data must correspond to number
%of levels
%   SETTIME(obj, jLevel, variableName, value, ...)
%   SETTIME(obj, ':', variableName, value, ...)
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
        if isa(val, 'Type')
            % Copy type argument
            currentType = val;
            valueList = repmat({nan}, length(jLevel), 1);
        else
            % Determine type by example data
            [currentType, valueList] = ...
                Type.determineTypeValue(length(jLevel), name, val);
        end
        % Store type
        if ~ismember(name, obj.label)
            obj.type.(name) = currentType;
        end
        % Store variable as absolute variable
        if ~ismember(name, obj.timeVariable)
            obj.timeVariable{end+1} = name;
        end
        % Store level-oriented data
        for k = 1:length(jLevel)
            obj.level(jLevel(k)).(name) = valueList{k};
        end
    end
end