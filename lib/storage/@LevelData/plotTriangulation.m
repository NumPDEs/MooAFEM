function ax = plotTriangulation(obj, jLevel)
%%PLOTTRIANGULATION creates (a series of) plots for triangulations (must be
%stored as pair c4n (double: nNodes x 2) and n4e (int: nElem x 3) or as
%pair coordinates (double: 2 x nNodes) and elements (int: 3 x nElem)
%   ax = PLOTTRIANGULATION(obj) plots all levels
%   ax = PLOTTRIANGULATION(obj, ':')
%   ax = PLOTTRIANGULATION(obj, jLevel)
%   ax = PLOTTRIANGULATION(obj, 'end')

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


    % Proceed input
    if nargin < 2
        jLevel = 1:obj.nLevel;
    elseif strcmp(jLevel, ':')
        jLevel = 1:obj.nLevel;
    elseif strcmp(jLevel, 'end')
        jLevel = obj.nLevel;
    end

    % Compute time per frame to obtain an overall duration of 20 seconds
    DURATION = 20;
    timePerLevel = DURATION / length(jLevel);

    % Create axes object
    ax = gca;

    % Plot every frame
    for k = 1:length(jLevel)
        obj.plotLevelTriangulation(ax, jLevel(k));
        pause(timePerLevel);
    end
end