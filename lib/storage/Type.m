classdef Type < handle
%%TYPE class reprenting the type of a level oriented data

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


    properties
        % Data type as in format specifier
        type
        % Shape of data (scalar, vector, or array)
        shape
        % Number of characters required for printing
        width
    end

    properties (Dependent)
        % Boolean determining if type is float
        isFloat (1,1) boolean
        % Boolean determining if type is scalar
        isScalar (1,1) boolean
    end

    properties (Access=private)
        % String representing type and formatting of floats
        STR_FLOAT = '.5e';
    end

    methods
        %% CONSTRUCTOR
        function obj = Type(name, value)
            %%TYPE creates an object of the Type class given an example
            %value of the corresponding data
            %   obj = TYPE(name, value)

            % Determine data type and width
            if isnan(value)
                obj.type = 'd';
                obj.width = 3;
            elseif isinteger(value)
                obj.type = 'd';
                S = cellstr(num2str(value));
                obj.width = max(cellfun(@(x) length(x), S(:)));
            elseif isfloat(value)
                obj.type = obj.STR_FLOAT;
                obj.width = 6 + str2double(obj.STR_FLOAT(2));
            elseif ischar(value)
                obj.type = 'c';
                obj.width = length(value);
            elseif isstring(value)
                obj.type = 's';
                obj.width = max(cellfun(@(x) length(x), value(:)));
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
            obj.width = max(obj.width, length(name));
        end

        function bool = get.isFloat(obj)
            bool = strcmp(obj.type, obj.STR_FLOAT);
        end

        function bool = get.isScalar(obj)
            bool = strcmp(obj.shape, 's');
        end
    end
end
