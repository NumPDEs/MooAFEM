function printItem(obj, jItem)
%%PRINTITEM prints information on the given list of item indices to command
%line
%   PRINTITEM(obj) defaults to final level
%   PRINTITEM(obj, jItem)
%
%   See also LevelDataCollection/printHeader

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


    % Print final item by default
    if nargin < 2
        jItem = obj.nItem;
    end

    % Print header in case of plotting the first level
    if obj.isInitialRun()
        obj.printHeader();
    end

    % Set separator variable
    obj.separator = '  ';

    % Iterate over given list of item indices
    for k = 1:length(jItem)
        % Extract data of all time variables
        data = cell(obj.nTimeVariable+1, 1);
        data{1} = jItem(k);
        ind = 1;
        for l = 1:obj.nTimeVariable
            ind = ind + 1;
            data{ind} = obj.item{jItem(k)}.level.(obj.timeVariable{l});
        end
        % Print information on current item to command line
        fprintf(obj.formatSpecifier, data{1:ind});
    end
end