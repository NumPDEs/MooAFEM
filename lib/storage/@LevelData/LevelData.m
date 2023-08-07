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
        isInitialRun (1,1) logical
    end

    properties (Dependent, Access=private)
        % Boolean vector determining whether variable is a scalar value
        isScalar (1,:) logical
        % String specifying output format in printLevel
        formatSpecifier (1,:) char
        % String specifying header format in printLevel
        headerSpecifier (1,:) char
    end

    methods
        %% CONSTRUCTOR
        function obj = LevelData(rootpath)
            %%LEVELDATA creates an object for storing level-oriented data,
            %file storage is initialised in the specified optional path
            %   leveldata = LEVELDATA(rootpath)

            % Set root path for file storage
            if nargin >= 1
                obj.root = rootpath;
            end
            ensureFolderExists(obj.root);
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
            % List of names of all variables
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
            % Name of folder for file storage
            rootpath = obj.root;
            if numel(rootpath) > 0 && ~strcmp(rootpath(end), '/')
                rootpath = [rootpath, '/'];
            end
            path = [rootpath, obj.problem, '_', obj.domain, '_', obj.method];
        end

        function file = get.filename(obj)
            % Name of file for file storage
            file = [obj.identifier];
        end

        function bool = get.isScalar(obj)
            % Creates boolean vector determining whether a variable
            % contains a scalar value per level only
            bool = false(1, obj.nVariable);
            for j = 1:obj.nVariable
                t = obj.type.(obj.label{j});
                if strcmp(t.shape, 's')
                    bool(j) = true;
                end
            end
        end

        function variables = get.scalarVariable(obj)
            % Returns list of names of all scalar variables
            variables = cell(obj.nScalarVariable, 1);
            k = 1;
            for j = 1:obj.nVariable
                if obj.isScalar(j)
                    variables{k} = obj.label{j};
                    k = k + 1;
                end
            end
        end

        function nVariable = get.nScalarVariable(obj)
            nVariable = nnz(obj.isScalar);
        end

        function variables = get.timeVariable(obj)
            % Returns list of names of all time variables
            variables = obj.timeVariable;
        end

        function nVariable = get.nTimeVariable(obj)
            nVariable = length(obj.timeVariable);
        end

        function variables = get.absoluteVariable(obj)
            % Returns list of names of all absolute variables
            variables = obj.absoluteVariable;
        end

        function nVariable = get.nAbsoluteVariable(obj)
            nVariable = length(obj.absoluteVariable);
        end

        function bool = get.isInitialRun(obj)
            % True if only one level is stored
            bool = (obj.nLevel <= 1);
        end

        function spec = get.headerSpecifier(obj)
            % Creates formatting string for the header of the output
            % to command line
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
            % Creates formatting string for printing to command line
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
        data = get(obj, jLevel, varargin)

        %% MODIFY LEVEL DATA
        set(obj, jLevel, varargin)
        setTime(obj, jLevel, varargin)
        setAbsolute(obj, jLevel, varargin)

        function append(obj, varargin)
            %%APPEND simplified addition of specified data for the
            %subsequent level, the cell array varargin must contain
            %pairs of variable names and data (numeric or string /
            %characters)
            %   APPEND(obj, varargin)
            assert(mod(length(varargin), 2) == 0, ...
                   'Invalid number of arguments');
            nlvl = obj.nLevel + 1;
            obj.set(nlvl, varargin{:});
        end

        function remove(obj, varargin)
            %%REMOVE removes the specified list of variables
            %   REMOVE(obj, varargin)
            obj.level = rmfield(obj.level, varargin{:});
            obj.type = rmfield(obj.type, varargin{:});
        end

        function removeLevel(obj, jLevel)
            %%REMOVELEVEL removes all data for the specified list jLevel
            %of level numbers
            %   REMOVELEVEL(obj, jLevel)
            for k = 1:length(jLevel)
                obj.level(jLevel(k)) = [];
            end
        end

        function removeNonscalar(obj)
            %%REMOVENONSCALAR removes all non-scalar data completely;
            %e.g., to reduce storage size
            %   REMOVENONSCALAR(obj)
            variableNonScalar = setdiff(obj.label, obj.scalarVariable);
            obj.remove(variableNonScalar);
        end

        %% OUTPUT LEVEL DATA
        printHeader(obj)
        printLevel(obj, jLevel)
        ax = plot(obj, xVariable, varargin)
        ax = plotTime(obj, xVariable, varargin)
        ax = plotAbsolute(obj, xVariable, varargin)
        ax = plotTriangulation(obj, jLevel)

        %% EXPORT LEVEL DATA
        saveToFile(obj, folder, file)
        saveToTable(obj, separator)
        plotToFile(obj, xVariable, varargin)
        plotTriangulationToFile(obj, jLevel)
    end

    %% AUXILIARY FUNCTIONS
    methods (Access = private)
        function idx = getIndex(obj, variableName)
            [~, idx] = ismember(variableName, obj.label);
        end

        function width = getWidth(obj, type)
            width = num2str(max(type.width, obj.minimalWidth));
        end

        ax = plotLevel(obj, plotFunction, xVariable, variableName)
        plotLevelTriangulation(obj, ax, jLevel)
    end
end
