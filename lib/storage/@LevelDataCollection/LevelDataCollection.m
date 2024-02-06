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
        metaData
        % Root folder for file storage
        root = 'results'
        % Cell array of LevelData objects
        item = cell(0)
    end

    properties (Dependent)
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

            obj.metaData = dictionary(...
                "problem", "problem", ...
                "domain", "domain", ...
                "method", "method", ...
                "identifier", "main", ...
                "hostname", string(getHostname()), ...
                "timestamp", string(datetime('now', 'Format', 'yyyy-MM-dd_HH:mm:ss')));

            % Set identifier
            if nargin >= 2
                obj.metaData("identifier") = identifier;
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

        function path = get.foldername(obj)
            rootpath = obj.root;
            if numel(rootpath) > 0 && ~strcmp(rootpath(end), '/')
                rootpath = [rootpath, '/'];
            end
            path = rootpath + strjoin([obj.metaData("problem"), ...
                                       obj.metaData("domain"), ...
                                       obj.metaData("method")], '_') + '/';
        end

        function file = get.filename(obj)
            file = obj.metaData("identifier") + '_collection';
        end

        function spec = getHeaderSpecifier(obj, separator)
            % Creates formatting string for the header of the output to command line
            spec = cell(1, obj.nTimeVariable+1);
            spec{1} = ' run';
            for j = 1:obj.nTimeVariable
                t = obj.item{1}.type(obj.timeVariable{j});
                spec{j+1} = assembleSpecifier(obj.getWidth(t), 's');
            end
            spec = strjoin(spec, separator) + "\n";
        end

        function spec = getFormatSpecifier(obj, separator)
            % Creates formatting string for printing to command line
            spec = cell(1, obj.nTimeVariable+1);
            spec{1} = '%4d';
            for j = 1:obj.nTimeVariable
                t = obj.item{1}.type(obj.timeVariable{j});
                spec{j+1} = assembleSpecifier(obj.getWidth(t), t.formatSpec);
            end
            spec = strjoin(spec, separator) + "\n";
        end

        %% READ ITEM DATA
        output = get(obj, jItem, varargin)

        %% MODIFY ITEMS
        set(obj, jItem, varargin)

        function append(obj, data)
            %%APPEND simplified addition specified of data to this list, one or
            %more LevelData objects
            %   APPEND(obj, data)

            arguments
                obj
            end

            arguments (Repeating)
                data LevelData
            end

            indices = obj.nItem + (1:length(data));
            obj.set(indices, data{:});
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
            width = num2str(max(type.printWidth, obj.minimalWidth));
        end

        printTable(obj, fid, variableName, data, separator)
    end
end

function spec = assembleSpecifier(width, format)
    spec = ['%', num2str(width), format];
end