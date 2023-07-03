function output = get(obj, jItem, varargin)
%%GET extracts data from this LevelDataCollection object for a given list
%jItem of item numbers, the cell array varargin must contain the names of 
%the variables to return
%   output = GET(obj, jItem, varargin)
%   output = GET(obj, ':', varargin)

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
    if strcmp(jItem, ':')
        jItem = 1:obj.nItem;
    end

    % Initialise output variable
    output = cell(1, length(varargin));

    % Determine number of levels (uses first item!)
    nLevel = obj.item{1}.nLevel;

    % Iterate over variables
    for jVariable = 1:length(varargin)
        variableName = varargin{jVariable};
        % Check for presence of variable
        if ~ismember(variableName, obj.timeVariable)
            warning(['Variable ', variableName, ' not found']);
            data = [];
        else
            data = zeros(nLevel, length(jItem));
            % Iterate over items
            for k = 1:length(jItem)
                data(:,k) = obj.item{jItem(k)}.get(':', variableName);
            end
        end
        % Store to return variable
        output{jVariable} = data;
    end
end