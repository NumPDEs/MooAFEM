function printStatistics(obj, variable)
%%PRINTSTATISTICS prints statististical information of scalar time
%variables (specified in variable)
%   PRINTSTATISTICS(obj, variable, ...)

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
    end

    arguments (Repeating)
        variable {mustBeTextScalar}
    end

    % Choose time variables for plotting
    if isempty(variable)
        variable = obj.timeVariable;
    else
        variable = intersect(variable, obj.timeVariable);
    end

    % Extract data from each item and on each level
    data = obj.get(':', variable{:});

    % Print separator between tables
    fprintf('\n');

    % Print table for each given variable
    for j = 1:length(variable)
        obj.printTable(1, variable{j}, data{j}, '  ');
        fprintf('\n\n');
    end
end


