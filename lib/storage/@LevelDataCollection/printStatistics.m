function printStatistics(obj, varargin)
%%PRINTSTATISTICS prints statististical information of scalar time
%variables (specified in varargin)
%   PRINTSTATISTICS(obj, varargin)

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


    % Choose time variables for plotting
    if nargin < 2
        varargin = obj.timeVariable;
    else
        varargin = intersect(varargin, obj.timeVariable);
    end

    % Specify separator for table
    obj.separator = '  ';

    % Extract data from each item and on each level
    data = obj.get(':', varargin{:});

    % Print separator between tables
    fprintf('\n');

    % Print table for each given variable
    for j = 1:length(varargin)
        obj.printTable(1, varargin{j}, data{j});
        fprintf('\n\n');
    end
end


