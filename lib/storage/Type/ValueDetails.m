classdef ValueDetails
%%ValueDetails class representing the type and printing parameters of a level
%oriented data

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


    properties (GetAccess=public, SetAccess=protected)
        % Raw type
        rawType
        % Shape of data (scalar, vector, or array)
        shape
        % Number of characters required for printing
        printWidth
    end

    properties (Dependent)
        % Boolean determining if type is scalar
        isScalar (1,1) boolean
    end

    methods
        %% CONSTRUCTOR
        function obj = ValueDetails(name, value)
            %%TYPE creates an object of the Type class given an example
            %value of the corresponding data
            %   obj = TYPE(name, value)

            % Determine data type and width
            if isnan(value)
                % NaN is treated as integer to avoid problems with plotting
                obj.rawType = RawType.INT;
                obj.printWidth = 3;
            elseif isinteger(value)
                obj.rawType = RawType.INT;
                S = cellstr(num2str(value));
                obj.printWidth = max(cellfun(@(x) length(x), S(:)));
            elseif isfloat(value)
                obj.rawType = RawType.FLOAT;
                floatPrecision = str2double(RawType.FLOAT.formatSpec(2));
                obj.printWidth = 6 + floatPrecision;
            elseif ischar(value)
                obj.rawType = RawType.TEXT;
                obj.printWidth = length(value);
            elseif isstring(value)
                obj.rawType = RawType.TEXT;
                obj.printWidth = max(cellfun(@(x) length(x), value(:)));
            else
                error('Invalid data type');
            end

            % Determine shape
            if isscalar(value)
                obj.shape = 's';
            elseif isvector(value)
                obj.shape = 'v';
            else
                obj.shape = 'a';
            end

            % Update width considering length of the name
            obj.printWidth = max(obj.printWidth, length(name));
        end

        function bool = get.isScalar(obj)
            bool = strcmp(obj.shape, 's');
        end
    end

    % AUXILIARY
    methods (Static)
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
                type = ValueDetails(variableName, value);
                valueList = {value};
            else
                if isvector(value)
                    % Array value represents a vector of scalar values for each
                    % level
                    assert(length(value) == nLevel, ...
                        'Length of argument invalid');
                    type = ValueDetails(variableName, value(1));
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
                        type = ValueDetails(variableName, value(1, ind{:}));
                        for k = 1:nLevel
                            valueList{k} = value(k, ind{:});
                        end
                    else
                        error('Unable to extract data from argument');
                    end
                end
            end
        end
    end
end
