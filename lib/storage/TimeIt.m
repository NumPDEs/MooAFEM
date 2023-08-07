function leveldatacollection = TimeIt(identifier, nRun, functionName, varargin)
%%TIMEIT function wrapper for multiple runs of functions storing timing
%data in LevelData objects
%   leveldatacollection = TIMEIT(identifier, nRun, functionName, varargin)

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

    % Proceed input
    if nargin < 4
        varargin = cell(0);
    end

    % Welcome statement
    fprintf('\n#\n# MooAFEM - TIME MEASUREMENT\n');
    fprintf('# Current Time: %s\n#\n\n', datestr(now));
    fprintf('This is TimeIt wrapper for function "%s"\n\n', functionName);

    % Initialisation of collection for LevelData objects
    leveldatacollection = LevelDataCollection();
    leveldatacollection.metaData("identifier") = identifier;

    % Run experiments
    for j = 1:nRun
        % Run current experiment
        temporaryIdentifier = sprintf('timingRun%04d', j);
        outputList = fevalc(functionName, varargin{:});
        leveldata = outputList{1};
        leveldata.metaData("identifier") = temporaryIdentifier;
        % Remove unused information
        leveldata.removeNonscalar();
        % Store information
        leveldatacollection.append(leveldata);
        % Save intermediate collection to file
        leveldatacollection.saveToFile();
        % Print information on current run
        leveldatacollection.printItem();
    end

    % Closing
    fprintf('\n');
end


function output = fevalc(functionName, varargin) %#ok<INUSD>
%%FEVALC suppresses output to commandline

    % Create function call
    functioncall = 'feval(functionName, varargin{:})';
    % Call function
    if isOctave()
        evalc(['output = ' functioncall]);
    else
        [~, output] = evalc(functioncall);
    end
    % Store output in cell variable
    if ~iscell(output)
        output = {output};
    end
end
