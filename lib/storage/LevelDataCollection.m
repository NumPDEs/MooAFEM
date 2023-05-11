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
 
        %% FILE STORAGE OF COMPLETE OBJECT
        function saveToFile(obj, folder, file)
            %%SAVETOFILE saves this object to a file with automatically
            %generated folder and file names if no optional arguments are
            %specified
            %   obj.SAVETOFILE()
            %   obj.SAVETOFILE(folder, file)

            % Proceed optional input
            if nargin < 2
                folder = obj.foldername;
            end
            if nargin < 3
                file = obj.filename;
            end
            % Create problem- and method-specific folder
            folderExists = makeDirectory(folder);
            if ~folderExists
                error('Could not create folder and save to file');
            end
            % Save this object to file
            if isOctave()
                save([folder, '/', file, '.mat'], 'obj');
            else
                save([folder, '/', file, '.mat'], 'obj', '-v7.3');
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
            spec = [' run', obj.separator];
            for j = 1:obj.nTimeVariable
                t = obj.item{1}.type.(obj.timeVariable{j});
                if j < obj.nTimeVariable
                    spec = [spec, '%', obj.getWidth(t), 's', obj.separator]; %#ok<*AGROW>
                else
                    spec = [spec, '%', obj.getWidth(t), 's\n'];
                end
            end
        end

        function spec = get.formatSpecifier(obj)
            spec = ['%4d', obj.separator];
            for j = 1:obj.nTimeVariable
                t = obj.item{1}.type.(obj.timeVariable{j});
                if j < obj.nTimeVariable
                    spec = [spec, '%', obj.getWidth(t), t.type, obj.separator];
                else
                    spec = [spec, '%', obj.getWidth(t), t.type, '\n'];
                end
            end
        end

        %% READ ITEM DATA
        function output = get(obj, jItem, varargin)
            if strcmp(jItem, ':')
                jItem = 1:obj.nItem;
            end
            output = cell(1, length(varargin));
            nLevel = obj.item{1}.nLevel;
            for jVariable = 1:length(varargin)
                variableName = varargin{jVariable};
                if ~ismember(variableName, obj.timeVariable)
                    warning(['Variable ', variableName, ' not found']);
                    data = [];
                else
                    data = zeros(nLevel, length(jItem));
                    for k = 1:length(jItem)
                        data(:,k) = obj.item{jItem(k)}.get(':', variableName);
                    end
                end
                output{jVariable} = data;
            end
        end

        %% MODIFY ITEMS
        function set(obj, jItem, varargin)
            assert(all(isa(varargin{:}, 'LevelData')), ...
                   'Arguments must be of class LevelData');
            assert(length(jItem) == length(varargin), ...
                   'Number of indices must equal number of given arguments');
            for k = 1:length(jItem)
                obj.item{jItem(k)} = varargin{k};
            end
        end

        function append(obj, varargin)
            assert(all(isa(varargin{:}, 'LevelData')), ...
                   'Arguments must be of class LevelData');
            indices = obj.nItem + (1:length(varargin));
            obj.set(indices, varargin{:});
        end

        function remove(obj, indices)
            indices = intersect(indices, 1:obj.nItem);
            obj.item{indices} = [];
        end

        %% OUTPUT ITEM DATA
        function printHeader(obj)
            obj.separator = '  ';
            header = sprintf(obj.headerSpecifier, obj.timeVariable{:});
            fprintf(header);
            fprintf([repmat('-', 1, length(header)-1), '\n']);
        end

        function printItem(obj, jItem)
            if nargin < 2
                jItem = obj.nItem;
            end
            obj.separator = '  ';
            for k = 1:length(jItem)
                data = cell(obj.nTimeVariable+1, 1);
                data{1} = jItem(k);
                ind = 1;
                for l = 1:obj.nTimeVariable
                    ind = ind + 1;
                    data{ind} = obj.item{jItem(k)}.level.(obj.timeVariable{l});
                end
                fprintf(obj.formatSpecifier, data{1:ind});
            end
        end

        function printStatistics(obj, variableName)
            obj.separator = '  ';
            if nargin < 2
                variableName = obj.timeVariable;
            else
                variableName = intersect(variableName, obj.timeVariable);
            end
            data = obj.get(':', variableName{:});
            fprintf('\n');
            for j = 1:length(variableName)
                obj.printTable(1, variableName{j}, data{j});
                fprintf('\n\n');
            end
        end

        function plotStatistics(obj, xVariable, variableName)
            obj.separator = '  ';
            if nargin < 3
                variableName = obj.timeVariable;
            else
                variableName = intersect(variableName, obj.timeVariable);
            end
            xData = obj.item{1}.get(':', xVariable);
            data = obj.get(':', variableName{:});
            ax = obj.plotData(xData, variableName, data);
            % Add title
            title(ax, 'Runtime plot');
            % Add axes labels
            xlabel(ax, xVariable);
            ylabel(ax, 'runtime (in sec.)');
            % Update legend
            configureLegend(ax, 'northwest');
        end

        %% EXPORT ITEM DATA
        function saveToTable(obj, separator)
            % Proceed optional input
            if nargin < 2
                obj.separator = ',';
            else
                obj.separator = separator;
            end
            % Create problem- and method-specific folder
            folderExists = makeDirectory(obj.foldername);
            if ~folderExists
                error('Could not create folder and save to file');
            end
            % Save table to file
            data = obj.get(':', obj.timeVariable{:});
            for j = 1:obj.nTimeVariable
                fid = fopen([obj.foldername, '/', obj.filename, '_', ...
                             obj.timeVariable{j},'.csv'], 'w');
                obj.printTable(fid, obj.timeVariable{j}, data{j});
                fclose(fid);
            end
        end
    end


    methods (Access = private)
        function width = getWidth(obj, type)
            width = num2str(max(type.width, obj.minimalWidth));
        end

        function printTable(obj, fid, variableName, data)
            TITLE = ['# TIME STATISTICS - ', variableName];
            HEADLINE = sprintf(['%5s', repmat([obj.separator, '%11s'], 1, 5)],...
                                'level', 'mean', 'median', 'std', 'min', 'max');
            if fid == 1
                fprintf(fid, [TITLE, '\n\n', HEADLINE, '\n', ...
                              repmat('-', 1, length(HEADLINE)), '\n']);
            else
                fprintf(fid, [HEADLINE, '\n']);
            end
            statisticalData = [mean(data, 2),...
                               median(data, 2),...
                               std(data, 0, 2),...
                               min(data, [], 2),...
                               max(data, [], 2)];
            fprintf(fid, ['%5d', repmat([obj.separator, '%8.5e'], 1, 5), '\n'],...
                    [1:size(data, 1); permute(statisticalData, [2 1])]);
        end
    end

    methods (Static, Access=private)
        function ax = plotData(xData, variableName, data)
            ax = gca;
            COLOURORDER = getPlotStyle();
            for j = 1:length(variableName)
                yMean = mean(data{j}, 2);
                yMin = min(data{j}, [], 2);
                yMax = max(data{j}, [], 2);
                % Extract label for legend from dictionary
                variableLabel = variableName{j};
                % Create plot
                if isOctave()
                    % fmt = [OCTAVE_FORMATS{nLines}, ';', variableLabel, ';'];
                    % errorbar(ax, xData, yMean, ...
                    %          yMean - yMin, yMax - yMean, fmt);
                else
                    errorbar(ax, xData, yMean, ...
                             yMean - yMin, yMax - yMean, ...
                             'LineStyle', '-', 'Marker', '.', ...
                             'Color', COLOURORDER(j,:), ...
                             'DisplayName', variableLabel);
                end
                % Use double-logarithmic scales
                set(ax, 'XScale', 'log', 'YScale', 'log');
                % Add new line into the current figure when calling plotConvergence again
                hold(ax, 'on');
            end
        end
    end
end