function data = get(obj, jLevel, varargin)
%%GET extracts data from this LevelData object on a given list jLevel of 
%level numbers, the cell array varargin must contain the names of the
%variables to return
%   data = GET(obj, jLevel, varargin)
%   data = GET(obj, ':', varargin)

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

    % Initialise return variable
    data = nan(length(jLevel), length(varargin));
    containsCharOnly = true;

    % Iterate over all variables
    for jVariable = 1:length(varargin)
        variableName = varargin{jVariable};
        % Check for existence of specified variable
        if ~ismember(variableName, obj.label)
            warning(['Variable ', variableName, ' not found']);
            value = [];
        else
            % Check for strings/character arrays
            if ~strcmp(obj.type.(variableName).type, 'c') && ...
                    ~strcmp(obj.type.(variableName).type, 's')
                containsCharOnly = false;
            end
            % Determine index of current variable
            idx = obj.getIndex(variableName);
            if obj.isScalar(idx)
                % Extract scalar variables for each level
                value = zeros(length(jLevel), 1);
                for k = 1:length(jLevel)
                    value(k) = obj.level(jLevel(k)).(variableName);
                end
            else
                % Extract non-scalar variables only in case of a single level
                if length(jLevel) > 1
                    warning('Export of level-oriented scalar variables only');
                    value = [];
                else
                    data(1,jVariable) = obj.level(jLevel).(variableName);
                    return
                end
            end
        end
        % Save extracted data to return variable
        data(:,jVariable) = value;
        % Post-process character arrays
        if containsCharOnly
            data = mat2cell(char(data), ones(length(data), 1));
        end
    end
end