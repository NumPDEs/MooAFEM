function set(obj, jLevel, variableName, value)
%%SET stores specified data to this LevelData object for the specified list
%jLevel of levels. The variableName/value pairs can be repeating, their first
%dimension must correspond to number of levels
%   SET(obj, jLevel, variableName, value, ...)
%   SET(obj, ':', variableName, value, ...)

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
        % Store type
        if ~ismember(name, obj.label)
            obj.type(name) = currentType.adaptPrintWidthTo(name);
        end
        % Store level-oriented data
        for k = 1:length(jLevel)
            obj.level(jLevel(k)).(name) = valueList{k};
        end
    end
end