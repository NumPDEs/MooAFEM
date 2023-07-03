function setTime(obj, jLevel, varargin)
%%SETTIME stores specified data of increasing type (time / ndof / ...) to 
%this LevelData object for the specified list jLevel of levels, the cell
%array varargin must contain pairs of variable names and data, dimension 
%one of data must correspond to number of levels
%   SETTIME(obj, jLevel, varargin)
%   SETTIME(obj, ':', varargin)
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


    % Proceed input
    if strcmp(jLevel, ':')
        jLevel = 1:obj.nLevel;
    end

    % Check number of variable input
    assert(mod(length(varargin), 2) == 0, ...
           'Invalid number of arguments');

    % Determine number of arguments
    nArgument = length(varargin) / 2;

    % Store data to LevelData object
    for j = 1:nArgument
        variableName = varargin{2*j-1};
        % Determine type of data
        if isa(varargin{2*j}, 'Type')
            % Copy type argument
            currentType = varargin{2*j};
            valueList = repmat({nan}, length(jLevel), 1);
        else
            % Determine type by example data
            [currentType, valueList] = ...
                Type.determineTypeValue(length(jLevel), ...
                                        variableName, ...
                                        varargin{2*j});
        end
        % Store type
        if ~ismember(variableName, obj.label)
            obj.type.(variableName) = currentType;
        end
        % Store variable as absolute variable
        if ~ismember(variableName, obj.timeVariable)
            obj.timeVariable{end+1} = variableName;
        end
        % Store level-oriented data
        for k = 1:length(jLevel)
            obj.level(jLevel(k)).(variableName) = valueList{k};
        end
    end
end