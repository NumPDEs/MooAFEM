function h = createStandardFigure()
%%CREATESTANDARDFIGURE creates a new figure object with default choice of
%size, the modification of this file allows a unified formatting of all
%figures for the output in the LevelData class
%   h = CREATESTANDARDFIGURE()

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

    % Create figure
    h = figure('Units', 'centimeters', 'Position', 3 + [0 0 32 24]);

end