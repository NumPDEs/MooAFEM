function saveToFile(obj, folder, file)
%%SAVETOFILE saves this object to a file with automatically
%generated folder and file names if no optional arguments are
%specified
%   obj.SAVETOFILE()
%   obj.SAVETOFILE(folder, file)

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


    % Proceed optional input
    if nargin < 2
        folder = obj.foldername;
    end
    if nargin < 3
        file = obj.filename;
    end

    % Create problem- and method-specific folder
    ensureFolderExists(folder);

    % Save this object to file
	save(folder + '/' + file + '.mat', 'obj', '-v7.3');
end
