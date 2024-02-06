function leveldatacollection = TimeIt(identifier, nRun, storeResults, functionName, arguments)
%%TIMEIT function wrapper for multiple runs of functions storing timing
%data in LevelData objects
%   leveldatacollection = TIMEIT(identifier, nRun, functionName, arguments, ...)

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

    arguments
        identifier {mustBeTextScalar}
        nRun (1,1) double
        storeResults (1,1) logical
        functionName {mustBeTextScalar}
    end

    arguments (Repeating)
        arguments
    end

    % Welcome statement
    fprintf('\n#\n# MooAFEM - TIME MEASUREMENT\n');
    fprintf('# Current Time: %s\n#\n\n', string(datetime));
    fprintf('This is TimeIt wrapper for function "%s"\n\n', functionName);

    % Initialisation of collection for LevelData objects
    leveldatacollection = LevelDataCollection();
    leveldatacollection.metaData("identifier") = identifier;

    % Run experiments
    for j = 1:nRun
        % Run current experiment
        temporaryIdentifier = sprintf('timingRun%04d', j);
        outputList = fevalc(functionName, arguments);
        leveldata = outputList{1};
        leveldata.metaData("identifier") = temporaryIdentifier;
        % Remove unused information
        leveldata.removeNonscalar();
        % Store information
        leveldatacollection.append(leveldata);
        if storeResults
            % Save intermediate collection to file
            leveldatacollection.saveToFile();
        end
        % Print information on current run
        leveldatacollection.printItem();
    end

    % Copy metadata from leveldata
    leveldatacollection.metaData("problem") = leveldata.metaData("problem");
    leveldatacollection.metaData("domain") = leveldata.metaData("domain");
    leveldatacollection.metaData("method") = leveldata.metaData("method");

    % Closing
    fprintf('\n');
end


function output = fevalc(functionName, arguments) %#ok<INUSD>
%%FEVALC suppresses output to commandline

    % Create function call
    functioncall = 'feval(functionName, arguments{:})';
    % Call function
	[~, output] = evalc(functioncall);
    % Store output in cell variable
    if ~iscell(output)
        output = {output};
    end
end
