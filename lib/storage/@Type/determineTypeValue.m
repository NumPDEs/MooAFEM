function [type, valueList] = determineTypeValue(nLevel, variableName, value)
%%DETERMINETYPEVALUE auxiliary function for determining type from an array
%of values
%    [type, valueList] = DETERMINETYPEVALUE(nLevel, variableName, value)

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


    if nLevel == 1
        % Array value represents data on a single level
        type = Type(variableName, value);
        valueList = {value};
    else
        if isvector(value)
            % Array value represents a vector of scalar values for each
            % level
            assert(length(value) == nLevel, ...
                   'Length of argument invalid');
            type = Type(variableName, value(1));
            valueList = mat2cell(value(:), ones(nLevel, 1));
        else
            % Array value represents an higher dimensional objects for each
            % level
            dims = size(value);
            if any(dims == nLevel)
                levelDim = find(dims == nLevel);
                remainingDim = 1:ndims(value);
                remainingDim(levelDim) = [];
                if sum(dims == nLevel) > 1
                    % Several dimensions could represent the level
                    % information
                    warning(['Unclear dimensions of argument. ', ...
                             'Assume level-oriented data in dimension ', ...
                             num2str(levelDim)]);
                end
                % Extract and rearrange values
                value = permute(value, [levelDim, remainingDim]);
                valueList = cell(nLevel, 1);
                ind = repmat({':'}, length(remainingDim), 1);
                type = Type(variableName, value(1, ind{:}));
                for k = 1:nLevel
                    valueList{k} = value(k, ind{:});
                end
            else
                error('Unable to extract data from argument');
            end
        end
    end
end