function configureLegend(ax, location)
%%CONFIGURELEGEND formats the legend of the current axis and displays it in
%the specified location
%   CONFIGURELEGEND(ax, location)

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


    % Create handle to legend in current axis
    leg = legend(ax, 'show');

    % Specify location
    set(leg, 'location', location);

    % Set interpreter to latex in Matlab
    if isOctave()
        set(leg, 'interpreter', 'tex');
    else
        set(leg, 'interpreter', 'latex');
    end
end