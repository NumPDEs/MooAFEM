function plotToFile(obj, xVariable, varargin)
%%PLOTTOFILE creates plots of various scalar variables (specified in 
%varargin) with respect %to the variable xVariable and stores it to a file,
%filename is generated automatically from information in LevelData object
%   PLOTTOFILE(obj, xVariable, varargin)
%
%   See also LevelData/plot

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


    % Create problem- and method-specific folder
    folderExists = makeDirectory(obj.foldername);
    if ~folderExists
        error('Could not create folder and save to file');
    end

    % Create figure object
    h = createStandardFigure();

    % Create plot
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