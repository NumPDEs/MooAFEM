function folderExists = makeDirectory(path)
%%MAKEDIRECTORY creates a directory if it does not exist yet
%   folderExists = MAKEDIRECTORY(path)

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


    % Create folder
    if isfolder(path)
        folderExists = true;
    else
        folderExists = mkdir(path);
    end

    % Throw warning
    if ~folderExists
        warning('Error creating folder')
    end
end