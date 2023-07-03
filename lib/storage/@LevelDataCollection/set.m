function set(obj, jItem, varargin)
%%SET stores specified list of LevelData objects to the indices determined
%in the jItem array
%   SET(obj, jItem, varargin)

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


    % Check type of input
    assert(all(isa(varargin{:}, 'LevelData')), ...
           'Arguments must be of class LevelData');

    % Number of items and item indices must coincide
    assert(length(jItem) == length(varargin), ...
           'Number of indices must equal number of given arguments');

    % Store items
    for k = 1:length(jItem)
        obj.item{jItem(k)} = varargin{k};
    end
end