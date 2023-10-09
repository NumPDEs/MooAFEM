function saveToTable(obj, separator)
%%SAVETOTABLE stores all level information in a text file, defaults to
%csv-format (separator = ','), filename is generated automatically from 
%information in LevelData object
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

    arguments
        obj
        separator {mustBeTextScalar} = ','
    end

    % Open file
    ensureFolderExists(obj.foldername);
    fid = fopen(obj.foldername + '/' + obj.filename + '.csv', 'w');

    % Save header to file
    specifier = obj.getHeaderSpecifier(separator);
    fprintf(fid, specifier, obj.scalarVariable{:});

    % Save information on each level to file
    specifier = obj.getFormatSpecifier(separator);
    for k = 1:obj.nLevel
        data = cell(obj.nScalarVariable, 1);
        for j = 1:obj.nScalarVariable
            data{j} = obj.level(k).(obj.scalarVariable{j});
        end
        fprintf(fid, specifier, data{:});
    end

    fclose(fid);
end
