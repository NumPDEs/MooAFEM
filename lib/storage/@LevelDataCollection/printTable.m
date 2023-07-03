function printTable(obj, fid, variableName, data)
%%PRINTTABLE auxiliary private function for printing statististical
%information for data of scalar time variables variableName to command
%line (fid=1) or file (specified by file identifier fid)
%   PRINTTABLE(obj, fid, variableName, data)

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


    % Define formatting strings for title and headline
    TITLE = ['# TIME STATISTICS - ', variableName];
    HEADLINE = sprintf(['%5s', repmat([obj.separator, '%11s'], 1, 5)],...
                        'level', 'mean', 'median', 'std', 'min', 'max');

    if fid == 1
        % Print to command line (including title)
        fprintf(fid, [TITLE, '\n\n', HEADLINE, '\n', ...
                      repmat('-', 1, length(HEADLINE)), '\n']);
    else
        % Print to file
        fprintf(fid, [HEADLINE, '\n']);
    end

    % Compute statistical information
    statisticalData = [mean(data, 2),...
                       median(data, 2),...
                       std(data, 0, 2),...
                       min(data, [], 2),...
                       max(data, [], 2)];

    % Print information
    fprintf(fid, ['%5d', repmat([obj.separator, '%8.5e'], 1, 5), '\n'],...
            [1:size(data, 1); permute(statisticalData, [2 1])]);
end