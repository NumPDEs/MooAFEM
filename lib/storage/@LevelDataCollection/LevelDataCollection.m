classdef LevelDataCollection < handle
%%LEVELDATACOLLECTION class representing a collection of LevelData objects

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
        % Identifier for computation
        identifier = ''
        % Root folder for file storage
        root = 'results'
        % Cell array of LevelData objects
        item = cell(0)
    end

    properties (Dependent)
        % Name of problem
        problem (1,:) char
        % Name of domain
        domain (1,:) char
        % Name of method
        method (1,:) char
        % Path to folder for file storage
        foldername (1,:) char
        % Name of file for storage
        filename (1,:) char
        % Number of items
        nItem (1,1) integer
        % Cell array of names of timing variables
        timeVariable (1,:) cell
        % Number of time variables
        nTimeVariable (1,1) int32
        % Boolean indicating whether data has been stored for a single level
        isInitialRun (1,1) boolean
    end

    properties (Dependent, Access=private)
        % String specifying output format in printItem
        formatSpecifier (1,:) char
        % String specifying header format in prinItem
        headerSpecifier (1,:) char
    end

    properties (Access=protected)
        % Separator for printing to command line
        separator
        % Minimal width of variables for printing in command line
        minimalWidth = 8;
    end

    methods
        %% CONSTRUCTOR
        function obj = LevelDataCollection(rootpath, identifier)
            %%LEVELDATACOLLECTION

            % Set root path for file storage
            if nargin >= 1
                obj.root = rootpath;
            end
            makeDirectory(obj.root);

            % Set identifier
            if nargin >= 2
                obj.identifier = identifier;
            end
        end

        %% GET FUNCTIONS FOR DEPENDENT VARIABLES
        function nItem = get.nItem(obj)
            nItem = length(obj.item);
        end

        function timeVariable = get.timeVariable(obj)
            if isempty(obj.item)
                timeVariable = {};
            else
                timeVariable = obj.item{1}.timeVariable;
            end
        end

        function nVariable = get.nTimeVariable(obj)
            nVariable = length(obj.timeVariable);
        end

        function bool = get.isInitialRun(obj)
            bool = (obj.nItem <= 1);
        end

        function name = get.problem(obj)
            if isempty(obj.item)
                name = '';
            else
                name = [obj.item{1}.problem];
            end
        end

        function name = get.domain(obj)
            if isempty(obj.item)
                name = '';
            else
                name = [obj.item{1}.domain];
            end
        end

        function name = get.method(obj)
            if isempty(obj.item)
                name = '';
            else
                name = [obj.item{1}.method];
            end
        end

        function path = get.foldername(obj)
            rootpath = obj.root;
            if numel(rootpath) > 0 && ~strcmp(rootpath(end), '/')
                rootpath = [rootpath, '/'];
            end
            path = [rootpath, obj.problem, '_', obj.domain, '_', obj.method];
        end

        function file = get.filename(obj)
            file = [obj.identifier, '_collection'];
        end

        function spec = get.headerSpecifier(obj)
            % Creates formatting string for the header of the output 
            % to command line
            spec = [' run', obj.separator];
            for j = 1:obj.nTimeVariable
                t = obj.item{1}.type.(obj.timeVariable{j});
                if j < obj.nTimeVariable
                    spec = [spec, '%', obj.getWidth(t), 's', obj.separator]; %#ok<*AGROW>
                else
                    spec = [spec, '%', obj.getWidth(t), 's\n'];
                end
            end
            if obj.nTimeVariable == 0
                spec = [spec, '\n'];
            end
        end

        function spec = get.formatSpecifier(obj)
            % Creates formatting string for printing to command line
            spec = ['%4d', obj.separator];
            for j = 1:obj.nTimeVariable
                t = obj.item{1}.type.(obj.timeVariable{j});
                if j < obj.nTimeVariable
                    spec = [spec, '%', obj.getWidth(t), t.type, obj.separator];
                else
                    spec = [spec, '%', obj.getWidth(t), t.type, '\n'];
                end
            end
            if obj.nTimeVariable == 0
                spec = [spec, '\n'];
            end
        end

        %% READ ITEM DATA
        output = get(obj, jItem, varargin)

        %% MODIFY ITEMS
        set(obj, jItem, varargin)

        function append(obj, varargin)
            %%APPEND simplified addition specified data to this list, the
            %cell array varargin must contain LevelData objects
            %   APPEND(obj, varargin)
            assert(all(isa(varargin{:}, 'LevelData')), ...
                   'Arguments must be of class LevelData');
            indices = obj.nItem + (1:length(varargin));
            obj.set(indices, varargin{:});
        end

        function remove(obj, indices)
            %%REMOVE removes the specified list of objects
            %   REMOVE(obj, indices)
            indices = intersect(indices, 1:obj.nItem);
            obj.item{indices} = [];
        end

        %% OUTPUT ITEM DATA
        printHeader(obj)
        printItem(obj, jItem)
        printStatistics(obj, variableName)
        plotStatistics(obj, xVariable, variableName)

        %% EXPORT COMPLETE DATA
        saveToFile(obj, folder, file)
        saveToTable(obj, separator)
    end

    methods (Access = private)
        function width = getWidth(obj, type)
            width = num2str(max(type.width, obj.minimalWidth));
        end

        printTable(obj, fid, variableName, data)
    end
end
