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
        % Format specification for printing
        formatSpec
        % Boolean determining if type is scalar
        isScalar (1,1) boolean
    end

    methods
        %% CONSTRUCTOR
        function obj = ValueDetails(rawType, shape, printWidth)
            %%VALUEDETAILS creates an instance of ValueDetails given an example
            %value of the corresponding data
            %   obj = VALUEDETAILS(rawType, shape, printWidth)

            obj.rawType = rawType;
            obj.shape = shape;
            obj.printWidth = printWidth;
        end

        function obj = adaptPrintWidthTo(obj, name)
            %ADAPTPRINTWIDTHTO updates the print width considering the length of
            %variable name 
            obj.printWidth = max(obj.printWidth, length(name));
        end

        function formatSpec = get.formatSpec(obj)
            formatSpec = obj.rawType.formatSpec;
        end

        function bool = get.isScalar(obj)
            bool = strcmp(obj.shape, 's');
        end
    end

    % AUXILIARY
    methods (Static)
        function details = fromExample(value)
            %%FROMEXAMPLE creates an instance of ValueDetails given an example
            %value of the corresponding data
            %   obj = FROMEXAMPLE(value)

            % Determine data type and width
            if isnan(value)
                % NaN is treated as integer to avoid problems with plotting
                rawType = RawType.INT;
                printWidth = 3;
            elseif isinteger(value)
                rawType = RawType.INT;
                S = cellstr(num2str(value));
                printWidth = max(cellfun(@(x) length(x), S(:)));
            elseif isfloat(value)
                rawType = RawType.FLOAT;
                floatPrecision = str2double(RawType.FLOAT.formatSpec(2));
                printWidth = 6 + floatPrecision;
            elseif ischar(value)
                rawType = RawType.TEXT;
                printWidth = length(value);
            elseif isstring(value)
                rawType = RawType.TEXT;
                printWidth = max(cellfun(@(x) length(x), value(:)));
            else
                error('Invalid data type');
            end

            % Determine shape
            if isscalar(value)
                shape = 's';
            elseif isvector(value)
                shape = 'v';
            else
                shape = 'a';
            end

            details = ValueDetails(rawType, shape, printWidth);
        end
    end
end
