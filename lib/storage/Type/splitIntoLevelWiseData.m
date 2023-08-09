
function [valueList, type] = splitIntoLevelWiseData(nLevel, value)
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
        type = ValueDetails.fromExample(value);
        valueList = {value};
    else
        if isvector(value)
            % Array value represents a vector of scalar values for each % level
            assert(length(value) == nLevel, 'Length of argument invalid');
            type = ValueDetails.fromExample(value(1));
            valueList = mat2cell(value(:), ones(nLevel, 1));
        else
            % Array value represents an higher dimensional objects for each % level
            dims = size(value);
            [levelDim, idx] = extractSingleLevelDimension(dims, nLevel);
            % Extract values
            valueList = cell(nLevel, 1);
            for k = 1:nLevel
                idx(levelDim) = k;
                valueList{k} = subsref(value, idx);
            end
            type = ValueDetails.fromExample(valueList{end});
        end
    end
end

function [levelDim, idx] = extractSingleLevelDimension(dims, nLevel)
    levelDim = find(dims == nLevel);
    if isempty(levelDim)
        error('Unable to extract data from argument');
    end

    if length(levelDim) > 1
        % Several dimensions could represent the level information
        levelDim = levelDim(1);
        warning('Unclear dimensions of argument. Assume level-oriented data in dimension %d.', levelDim);
    end

    % index that is ':' for all dimensions, and where the index in the level
    % dimension can be changed
    idx.subs = repmat({':'}, 1, ndims(value));
    idx.type = '()';
end