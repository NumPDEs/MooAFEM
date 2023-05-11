classdef LevelData < handle
%%LEVELDATA class representing results from level-oriented computations

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
        % Name of problem
        problem = 'problem'
        % Name of domain
        domain = 'domain'
        % Name of method
        method = 'method'
        % Identifier for computation
        identifier = 'main'
        % Root folder for file storage
        root = 'results'
        % Structure array storing string representations of variables
        dictionary
        % Time stamp of creation
        timestamp
        % Hostname of computing machine
        hostname
        % Cell array of level data
        level (1,:) struct
        % Structure array of data types
        type (1,1) struct
        % Cell array of names of absolute value variables
        absoluteVariable (1,:) cell
        % Cell array of names of timing variables
        timeVariable (1,:) cell
    end

    properties (Access=private)
        % Separator for printing to command line
        separator
        % Minimal width of variables for printing in command line
        minimalWidth = 8;
    end

    properties (Dependent)
        % Path to folder for file storage
        foldername (1,:) char
        % Name of file for storage
        filename (1,:) char
        % Cell array of variable names
        label (1,:) cell
        % Number of levels
        nLevel (1,1) int32
        % Number of variables
        nVariable (1,1) int32
        % Number of scalar variables
        nScalarVariable (1,1) int32
        % Cell array of names of scalar variables
        scalarVariable (1,:) cell
        % Number of absolute value variables
        nAbsoluteVariable (1,1) int32
        % Number of time variables
        nTimeVariable (1,1) int32
        % Boolean indicating whether data has been stored for a single level
        isInitialRun (1,1) boolean
    end

    properties (Dependent, Access=private)
        % Boolean vector determining whether variable is a scalar value
        isScalar (1,:) boolean
        % String specifying output format in printLevel
        formatSpecifier (1,:) char
        % String specifying header format in printLevel
        headerSpecifier (1,:) char
    end

    methods
        %% CONSTRUCTOR
        function obj = LevelData(rootpath)
            %%LEVELDATA

            % Set root path for file storage
            if nargin >= 1
                obj.root = rootpath;
            end
            makeDirectory(obj.root);
            % Set hostname
            obj.hostname = getHostname();
            % Save time of creation
            if isOctave()
                obj.timestamp = datestr(now, 'yyyy-MM-dd_HH:mm:ss'); %#ok<TNOW1,DATST>
            else
                obj.timestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HH:mm:ss'));
            end
            % Initialise dictionary with some default values
            obj.dictionary = getDefaultDictionary();
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
                save([folder, '/', file, '.mat'], 'obj', '-v7');
            else
                save([folder, '/', file, '.mat'], 'obj', '-v7.3');
            end
        end

        %% SET GLOBAL VARIABLES
        function set.problem(obj, name)
            assert(ischar(name), 'Insert character array as problem name.');
            obj.problem = name;
        end

        function set.domain(obj, name)
            assert(ischar(name), 'Insert character array as domain name.');
            obj.domain = name;
        end

        function set.method(obj, name)
            assert(ischar(name), 'Insert character array as method name.');
            obj.method = name;
        end

        function set.root(obj, path)
            assert(ischar(path), 'Insert character array as path name.');
            obj.root = path;
        end

        %% GET GLOBAL VARIABLES
        function label = get.label(obj)
            if obj.nLevel > 0
                label = fieldnames(obj.level(1));
            else
                label = {};
            end
        end

        function nLevel = get.nLevel(obj)
            nLevel = length(obj.level);
        end

        function nVariable = get.nVariable(obj)
            nVariable = length(obj.label);
        end

        function path = get.foldername(obj)
            rootpath = obj.root;
            if numel(rootpath) > 0 && ~strcmp(rootpath(end), '/')
                rootpath = [rootpath, '/'];
            end
            path = [rootpath, obj.problem, '_', obj.domain, '_', obj.method];
        end

        function file = get.filename(obj)
            file = [obj.identifier];
        end

        function bool = get.isScalar(obj)
            bool = false(1, obj.nVariable);
            for j = 1:obj.nVariable
                t = obj.type.(obj.label{j});
                if strcmp(t.shape, 's')
                    bool(j) = true;
                end
            end
        end

        function variables = get.scalarVariable(obj)
            variables = cell(obj.nScalarVariable, 1);
            k = 1;
            for j = 1:obj.nVariable
                if obj.isScalar(j)
                    variables{k} = obj.label{j};
                    k = k + 1;
                end
            end
        end

        function variables = get.timeVariable(obj)
            variables = obj.timeVariable;
        end

        function nVariable = get.nScalarVariable(obj)
            nVariable = nnz(obj.isScalar);
        end

        function nVariable = get.nAbsoluteVariable(obj)
            nVariable = length(obj.absoluteVariable);
        end

        function nVariable = get.nTimeVariable(obj)
            nVariable = length(obj.timeVariable);
        end

        function bool = get.isInitialRun(obj)
            bool = (obj.nLevel <= 1);
        end

        function spec = get.headerSpecifier(obj)
            spec = '';
            for j = 1:obj.nScalarVariable
                t = obj.type.(obj.scalarVariable{j});
                if j < obj.nScalarVariable
                    spec = [spec, '%', obj.getWidth(t), 's', obj.separator]; %#ok<*AGROW>
                else
                    spec = [spec, '%', obj.getWidth(t), 's\n'];
                end
            end
        end

        function spec = get.formatSpecifier(obj)
            spec = '';
            for j = 1:obj.nScalarVariable
                t = obj.type.(obj.scalarVariable{j});
                if j < obj.nScalarVariable
                    spec = [spec, '%', obj.getWidth(t), t.type, obj.separator];
                else
                    spec = [spec, '%', obj.getWidth(t), t.type, '\n'];
                end
            end
        end

        %% READ LEVEL DATA
        function data = get(obj, jLevel, varargin)
            if strcmp(jLevel, ':')
                jLevel = 1:obj.nLevel;
            end
            data = nan(length(jLevel), length(varargin));
            containsCharOnly = true;
            for jVariable = 1:length(varargin)
                variableName = varargin{jVariable};
                if ~ismember(variableName, obj.label)
                    warning(['Variable ', variableName, ' not found']);
                    value = [];
                else
                    if ~strcmp(obj.type.(variableName).type, 'c') && ...
                            ~strcmp(obj.type.(variableName).type, 's')
                        containsCharOnly = false;
                    end
                    idx = obj.getIndex(variableName);
                    if obj.isScalar(idx)
                        value = zeros(length(jLevel), 1);
                        for k = 1:length(jLevel)
                            value(k) = obj.level(jLevel(k)).(variableName);
                        end
                    else
                        if length(jLevel) > 1
                            warning('Export of level-oriented scalar variables only');
                            value = [];
                        else
                            data(1,jVariable) = obj.level(jLevel).(variableName);
                            return
                        end
                    end
                end
                data(:,jVariable) = value;
                if containsCharOnly
                    data = mat2cell(char(data), ones(length(data), 1));
                end
            end
        end

        %% MODIFY LEVEL DATA
        function set(obj, jLevel, varargin)
            assert(mod(length(varargin), 2) == 0, ...
                   'Invalid number of arguments');
            if strcmp(jLevel, ':')
                jLevel = 1:obj.nLevel;
            end
            nArgument = length(varargin) / 2;
            for j = 1:nArgument
                variableName = varargin{2*j-1};
                value = varargin{2*j};
                if ~ismember(variableName, obj.label)
                    obj.type.(variableName) = Type(variableName, value(1,:,:));
                end
                for k = 1:length(jLevel)
                    obj.level(jLevel(k)).(variableName) = value(k,:,:);
                end
            end
        end

        function setTime(obj, jLevel, varargin)
            assert(mod(length(varargin), 2) == 0, 'Invalid number of arguments');
            if strcmp(jLevel, ':')
                jLevel = 1:obj.nLevel;
            end
            nArgument = length(varargin) / 2;
            for j = 1:nArgument
                variableName = varargin{2*j-1};
                value = varargin{2*j};
                assert(length(jLevel) == size(value, 1), ...
                       ['Invalid size of argument: ', variableName])
                if ~ismember(variableName, obj.label)
                    obj.type.(variableName) = Type(variableName, value(jLevel(1),:,:));
                end
                if ~ismember(variableName, obj.timeVariable)
                    obj.timeVariable{end+1} = variableName;
                end
                for k = 1:length(jLevel)
                    obj.level(jLevel(k)).(variableName) = value(k);
                end
            end
        end

        function setAbsolute(obj, jLevel, varargin)
            assert(mod(length(varargin), 2) == 0, 'Invalid number of arguments');
            nArgument = length(varargin) / 2;
            for j = 1:nArgument
                variableName = varargin{2*j-1};
                value = varargin{2*j};
                if ~ismember(variableName, obj.label)
                    obj.type.(variableName) = Type(variableName, value);
                end
                if ~ismember(variableName, obj.absoluteVariable)
                    obj.absoluteVariable{end+1} = variableName;
                end
                for k = 1:length(jLevel)
                    obj.level(jLevel(k)).(variableName) = value(k);
                end
            end
        end

        function append(obj, varargin)
            assert(mod(length(varargin), 2) == 0, ...
                   'Invalid number of arguments');
            nlvl = obj.nLevel + 1;
            obj.set(nlvl, varargin{:});
        end

        function remove(obj, varargin)
            obj.level = rmfield(obj.level, varargin{:});
            obj.type = rmfield(obj.type, varargin{:});
        end

        function removeLevel(obj, jLevel)
            for k = 1:length(jLevel)
                obj.level(jLevel(k)) = [];
            end
        end

        function removeNonscalar(obj)
            variableNonScalar = setdiff(obj.label, obj.scalarVariable);
            obj.remove(variableNonScalar);
        end

        %% OUTPUT LEVEL DATA
        function printHeader(obj)
            obj.separator = '  ';
            header = sprintf(obj.headerSpecifier, obj.scalarVariable{:});
            fprintf(header);
            fprintf([repmat('-', 1, length(header)-1), '\n']);
        end

        function printLevel(obj, jLevel)
            if nargin < 2
                jLevel = obj.nLevel;
            end
            if obj.isInitialRun()
                obj.printHeader();
            end
            obj.separator = '  ';
            for k = 1:length(jLevel)
                data = cell(obj.nVariable, 1);
                ind = 0;
                for jVariable = 1:obj.nVariable
                    if obj.isScalar(jVariable)
                        ind = ind + 1;
                        data{ind} = obj.level(jLevel(k)).(obj.label{jVariable});
                    end
                end
                fprintf(obj.formatSpecifier, data{1:ind});
            end
        end

        function ax = plot(obj, xVariable, varargin)
            if nargin >= 3
                variableName = setdiff(varargin, xVariable);
            else
                variableName = setdiff(obj.label, xVariable);
            end
            variableName = intersect(variableName, obj.scalarVariable);
            variableName = setdiff(variableName, obj.absoluteVariable);
            variableName = setdiff(variableName, obj.timeVariable);
            ax = obj.plotLevel(xVariable, variableName);
            % Add title
            title(ax, 'Convergence history plot');
            % Add axes labels
            xlabel(ax, xVariable);
            ylabel(ax, 'error');
            % Update legend
            configureLegend(ax, 'northeast');
        end

        function ax = plotTime(obj, xVariable, varargin)
            if nargin >= 3
                variableName = setdiff(varargin, xVariable);
            else
                variableName = setdiff(obj.label, xVariable);
            end
            variableName = intersect(variableName, obj.scalarVariable);
            variableName = intersect(variableName, obj.timeVariable);
            ax = obj.plotLevel(xVariable, variableName);
            % Add title
            title(ax, 'Time plot');
            % Add axes labels
            xlabel(ax, xVariable);
            ylabel(ax, 'runtime');
            % Update legend
            configureLegend(ax, 'northwest');
        end

        function ax = plotAbsolute(obj, xVariable, varargin)
            if nargin >= 3
                variableName = setdiff(varargin, xVariable);
            else
                variableName = setdiff(obj.label, xVariable);
            end
            variableName = intersect(variableName, obj.scalarVariable);
            variableName = intersect(variableName, obj.absoluteVariable);
            ax = obj.plotLevelAbsolute(xVariable, variableName);
            % Add title
            title(ax, 'Value plot');
            % Add axes labels
            xlabel(ax, xVariable);
            ylabel(ax, 'value');
            % Update legend
            configureLegend(ax, 'northeast');
        end

        function ax = plotTriangulation(obj, jLevel)
            if nargin < 2
                jLevel = 1:obj.nLevel;
            elseif strcmp(jLevel, 'end')
                jLevel = obj.nLevel;
            end
            % Compute time per frame to obtain an overall duration
            % of 20 seconds
            DURATION = 20;
            timePerLevel = DURATION / length(jLevel);
            % Create axes object
            ax = gca;
            % Plot every frame
            for k = 1:length(jLevel)
                ax = obj.plotLevelTriangulation(ax, jLevel(k));
                pause(timePerLevel);
            end
        end

        %% Export level data
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
            fid = fopen([obj.foldername, '/', obj.filename, '.csv'], 'w');
            fprintf(fid, obj.headerSpecifier, obj.scalarVariable{:});
            for k = 1:obj.nLevel
                data = cell(obj.nScalarVariable, 1);
                for j = 1:obj.nScalarVariable
                    data{j} = obj.level(k).(obj.scalarVariable{j});
                end
                fprintf(fid, obj.formatSpecifier, data{:});
            end
            fclose(fid);
        end

        function plotToFile(obj, xVariable, varargin)
            % Create problem- and method-specific folder
            folderExists = makeDirectory(obj.foldername);
            if ~folderExists
                error('Could not create folder and save to file');
            end
            % Create plot
            h = figure('Units', 'centimeters', 'Position',3 + [0 0 32 24]);
            if nargin < 3
                obj.plot(xVariable);
            else
                obj.plot(xVariable, varargin);
            end
            % Export plot
            print(h, '-dpng', '-r600', [obj.foldername, '/', obj.filename, '.png']);
            % Close figure
            close(h);
        end

        function plotTriangulationToFile(obj, jLevel)
            if nargin < 2
                jLevel = 1:obj.nLevel;
            elseif strcmp(jLevel, 'end')
                jLevel = obj.nLevel;
            end
            % Create problem- and method-specific folder
            folderExists = makeDirectory(obj.foldername);
            if ~folderExists
                error('Could not create folder and save to file');
            end
            if length(jLevel) == 1
                % Create single image
                % Create figure
                h = figure('Units', 'centimeters', 'Position',3 + [0 0 32 24]);
                obj.plotLevelTriangulation(gca(), jLevel);
                % Export plot
                print(h, '-dpng', '-r600', [obj.foldername, '/', obj.filename, '.png']);
                % Close figure
                close(h);
            else
                % Set parameters
                DURATION = 20;
                if isOctave()
                    obj.createTriangulationGIF(jLevel, DURATION);
                else
                    obj.createTriangulationAVI(jLevel, DURATION);
                end
            end
        end

    end

    %% AUXILIARY1
    methods (Access = private)
        function label = getDictionaryEntry(obj, variableName)
            if isfield(obj.dictionary, variableName)
                label = obj.dictionary.(variableName);
            else
                label = variableName;
            end
        end

        function idx = getIndex(obj, variableName)
            [~, idx] = ismember(variableName, obj.label);
        end

        function width = getWidth(obj, type)
            width = num2str(max(type.width, obj.minimalWidth));
        end

        function ax = plotLevel(obj, xVariable, variableName)
            ax = gca;
            xValue = obj.get(1:obj.nLevel, xVariable);
            [COLOURORDER, MARKER] = getPlotStyle();
            for j = 1:length(variableName)
                if obj.type.(variableName{j}).isFloat
                    yValue = obj.get(1:obj.nLevel, variableName{j});
                    % Extract label for legend from dictionary
                    variableLabel = obj.getDictionaryEntry(variableName{j});
                    % Create plot
                    loglog(ax, xValue, yValue, '-', ...
                           'Marker', MARKER{mod(j, length(MARKER))}, ...
                           'Color', COLOURORDER(j, :), ...
                           'LineWidth', 1.5, ...
                           'DisplayName', variableLabel);
                    % Add new line into the current figure when calling plotConvergence again
                    hold(ax, 'on');
                end
            end
        end

        function ax = plotLevelAbsolute(obj, xVariable, variableName)
            ax = gca;
            xValue = obj.get(1:obj.nLevel, xVariable);
            [COLOURORDER, MARKER] = getPlotStyle();
            for j = 1:length(variableName)
                if obj.type.(variableName{j}).isFloat
                    yValue = obj.get(1:obj.nLevel, variableName{j});
                    % Extract label for legend from dictionary
                    variableLabel = obj.getDictionaryEntry(variableName{j});
                    % Create plot
                    semilogx(ax, xValue, yValue, '-', ...
                             'Marker', MARKER{mod(j, length(MARKER))}, ...
                             'Color', COLOURORDER(j, :), ...
                             'LineWidth', 1.5, ...
                             'DisplayName', variableLabel);
                    % Add new line into the current figure when calling plotConvergence again
                    hold(ax, 'on');
                end
            end
        end

        function ax = plotLevelTriangulation(obj, ax, jLevel)
            assert(numel(jLevel) == 1);
            % Extract mesh information
            if isfield(obj.level(jLevel), 'c4n')
                c4n = obj.level(jLevel).('c4n');
                n4e = obj.level(jLevel).('n4e');
                % Reorganise mesh information
                c4e = permute(reshape(c4n(n4e,:), length(n4e), 3, 2), ...
                              [3 2 1]);
            else
                for k = 1:length(obj.label)
                    data = obj.level(jLevel).(obj.label(k));
                    if isa(data, 'MeshData')
                        c4e = data.c4e;
                        break
                    end
                end
            end
            % Plot mesh
            edgecolor = [0 0.4470 0.7410];
            patch(ax, permute(c4e(1,:,:), [2 3 1]), permute(c4e(2,:,:), [2 3 1]), ...
                  'white', 'EdgeColor', edgecolor);
            % Format plot
            title(ax, {'Mesh plot'; [num2str(size(c4e, 3)), ' triangles']});
            axis(ax, 'equal');
            axis(ax, 'tight');
        end

        function createTriangulationAVI(obj, jLevel, duration)
            % Set video parameter
            FPS = 30;
            % Create video writer
            videowriter = VideoWriter([obj.foldername, '/', obj.filename, '.avi']);
            set(videowriter, 'FrameRate', FPS);
            open(videowriter);
            % Compute time per frame to obtain a fixed overall duration
            framesPerLevel = fix(FPS * DURATION / length(jLevel));
            ax = gca();
            set(ax, 'NextPlot', 'replaceChildren');
            % Plot every frame
            for k = 1:length(jLevel)
                obj.plotLevelTriangulation(ax, jLevel(k));
                for l = 1:framesPerLevel
                    writeVideo(videowriter, getframe(h));
                end
            end
            % Close video writer
            close(videowriter);
            % Close figure
            close(h);
        end

        function createTriangulationGIF(obj, jLevel, duration)
            error('not working yet')
            % Create figure
            h = figure('Units', 'centimeters', 'Position',3 + [0 0 32 24]);
            % Assign plot to a frame
            frame = getframe(h);
            % Convert frame to RGB image (3 dimensional)
            im = frame2im(frame);
            % Transform RGB samples to 1 dimension with a color map "cm".
            [imind, cm] = rgb2ind(im);
            % Compute time per frame to obtain a fixed overall duration
            framesPerLevel = fix(FPS * DURATION / length(jLevel));
            imwrite(imind,cm,filename,'gif','DelayTime', DelayTime , 'Compression' , 'lzw');
            ax = gca();
            set(ax, 'NextPlot', 'replaceChildren');
            % Plot every frame
            for k = 1:length(jLevel)
                obj.plotLevelTriangulation(ax, jLevel(k));
                for l = 1:framesPerLevel
                    writeVideo(videowriter, getframe(h));
                end
                % Create GIF file
                % Add each new plot to GIF
                imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', DelayTime , 'Compression' , 'lzw');
            end
        end
    end
end
