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
        % map for general metadata
        metaData (1,1) dictionary
        % Root folder for file storage
        root
        % Structure array storing string representations of variables
        dictionary
        % Array of level data
        level (1,:) struct
        % Dictionary of data types
        type (1,1) dictionary
        % Cell array of names of absolute value variables
        absoluteVariable (1,:) cell
        % Cell array of names of timing variables
        timeVariable (1,:) cell
    end

    properties (Access=private)
        % Dictionary of categories for cataloguing and plotting
        category (1,1) dictionary
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
    end

    methods
        %% CONSTRUCTOR
        function obj = LevelData(rootpath)
            %%LEVELDATA creates an object for storing level-oriented data,
            %file storage is initialised in the specified optional path
            %   leveldata = LEVELDATA(rootpath)

            arguments
                rootpath {mustBeTextScalar} = 'results'
            end
            obj.root = rootpath;

            % TODO: set output of hostname to string and remove cast
            obj.metaData = dictionary(...
                "problem", "problem", ...
                "domain", "domain", ...
                "method", "method", ...
                "identifier", "main", ...
                "hostname", string(getHostname()));

            % Save time of creation
            if isOctave()
                obj.metaData("timestamp") = datestr(now, 'yyyy-MM-dd_HH:mm:ss'); %#ok<TNOW1,DATST>
            else
                obj.metaData("timestamp") = char(datetime('now', 'Format', 'yyyy-MM-dd_HH:mm:ss'));
            end

            obj.category = dictionary();

            ensureFolderExists(obj.root);
            % Initialise dictionary with some default values
            obj.dictionary = getDefaultDictionary();
        end

        %% SET GLOBAL VARIABLES
        function set.root(obj, path)
            arguments
                obj
                path {mustBeTextScalar}
            end
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
            path = rootpath + strjoin([obj.metaData("problem"), obj.metaData("domain"), obj.metaData("method")], '_');
        end

        function file = get.filename(obj)
            % Name of file for file storage
            file = [obj.metaData("identifier")];
        end

        function bool = get.isScalar(obj)
            % Creates boolean vector determining whether a variable
            % contains a scalar value per level only
            bool = false(1, obj.nVariable);
            for j = 1:obj.nVariable
                t = obj.type(obj.label{j});
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


        %% READ LEVEL DATA
        data = get(obj, jLevel, variableName)

        %% MODIFY LEVEL DATA
        set(obj, jLevel, variableName, value)
        setTime(obj, jLevel, variableName, value)
        setAbsolute(obj, jLevel, variableName, value)

        function append(obj, varargin)
            %%APPEND simplified addition of specified data for the
            %subsequent level, the arguments must be specified as
            %pairs of variable names and data (numeric or string /
            %characters)
            %   APPEND(obj, variableName, value, ...)

            arguments
                obj
            end

            arguments (Repeating)
                varargin
            end

            assert(mod(length(varargin), 2) == 0, ... 
                   'Invalid number of arguments');
            nlvl = obj.nLevel + 1;
            obj.set(nlvl, varargin{:});
        end

        function remove(obj, variableName)
            %%REMOVE removes the specified list of variables
            %   REMOVE(obj, variableName, ...)

            arguments
                obj
            end

            arguments (Repeating)
                variableName {mustBeTextScalar}
            end

            obj.level = rmfield(obj.level, variableName);
            obj.type(variableName) = [];
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
            obj.remove(variableNonScalar{:});
        end

        %% OUTPUT LEVEL DATA
        printHeader(obj)
        printLevel(obj, jLevel)
        ax = plot(obj, xVariable, yVariable)
        ax = plotTime(obj, xVariable, yVariable)
        ax = plotAbsolute(obj, xVariable, yVariable)
        ax = plotTriangulation(obj, jLevel)

        %% EXPORT LEVEL DATA
        saveToFile(obj, folder, file)
        saveToTable(obj, separator)
        plotToFile(obj, xVariable, yVariable)
        plotTriangulationToFile(obj, jLevel)
    end

    %% AUXILIARY FUNCTIONS
    methods (Access = private)
        function idx = getIndex(obj, variableName)
            [~, idx] = ismember(variableName, obj.label);
        end

        function width = getWidth(~, type)
            minimalWidth = 8;
            width = num2str(max(type.printWidth, minimalWidth));
        end

        function spec = getHeaderSpecifier(obj, separator)
            % Creates formatting string for the header of the output to command line
            spec = cell(1, obj.nScalarVariable);
            for j = 1:obj.nScalarVariable
                t = obj.type(obj.scalarVariable{j});
                spec{j} = assembleSpecifier(obj.getWidth(t), 's');
            end
            spec = strjoin(spec, separator) + "\n";
        end

        function spec = getFormatSpecifier(obj, separator)
            % Creates formatting string for printing to command line
            spec = cell(1, obj.nScalarVariable);
            for j = 1:obj.nScalarVariable
                t = obj.type(obj.scalarVariable{j});
                spec{j} = assembleSpecifier(obj.getWidth(t), t.formatSpec);
            end
            spec = strjoin(spec, separator) + "\n";
        end

        ax = plotLevel(obj, plotFunction, xVariable, variableName)
        plotLevelTriangulation(obj, ax, jLevel)
    end
end

function spec = assembleSpecifier(width, format)
    spec = ['%', num2str(width), format];
end