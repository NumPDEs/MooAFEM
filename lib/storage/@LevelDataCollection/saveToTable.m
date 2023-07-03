function saveToTable(obj, separator)
%%SAVETOTABLE stores statistical information for all time variables to a
%text file, defaults to csv-format (separator = ','), filename is generated
%automatically from information in LevelData object
%   SAVETOTABLE(obj)
%   SAVETOTABLE(obj, separator)

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
        obj.separator = ',';
    else
        obj.separator = separator;
    end

    % Create problem- and method-specific folder
    folderExists = makeDirectory(obj.foldername);
    if ~folderExists
        error('Could not create folder and save to file');
    end

    % Save data for each variable to a separate file
    data = obj.get(':', obj.timeVariable{:});
    for j = 1:obj.nTimeVariable
        fid = fopen([obj.foldername, '/', obj.filename, '_', ...
                     obj.timeVariable{j},'.csv'], 'w');
        obj.printTable(fid, obj.timeVariable{j}, data{j});
        fclose(fid);
    end
end