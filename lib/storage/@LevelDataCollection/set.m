function set(obj, jItem, data)
%%SET stores specified list of LevelData objects to the indices determined
%in the jItem array
%   SET(obj, jItem, data)

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
        jItem
    end

    arguments (Repeating)
        data LevelData
    end

    % Number of items and item indices must coincide
    assert(length(jItem) == length(data), ...
           'Number of indices must equal number of given arguments');

    % Store items
    for k = 1:length(jItem)
        obj.item{jItem(k)} = data{k};
    end
end