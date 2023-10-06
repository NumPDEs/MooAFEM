function ax = plotTriangulation(obj, jLevel, opt)
%%PLOTTRIANGULATION creates (a series of) plots for triangulations (must be
%stored as pair c4n (double: nNodes x 2) and n4e (int: nElem x 3) or as
%pair coordinates (double: 2 x nNodes) and elements (int: 3 x nElem)
%   ax = PLOTTRIANGULATION(obj) plots all levels
%   ax = PLOTTRIANGULATION(obj, jLevel) plots levels specified in jLevel
%   ax = PLOTTRIANGULATION(obj, jLevel, opt) additionally specifies options:
%      - timePerLevel (default 1s)

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
        jLevel (1,:) {mustBeIndexVector} = ':'
        opt.timePerLevel {mustBeNumeric} = 1
    end

    if strcmp(jLevel, ':')
        jLevel = 1:obj.nLevel;
    end

    % Create axes object
    ax = gca;

    % Plot every frame
    % TODO: store whole mesh such that also boundary can be plotted?
    for k = jLevel
        coordinates = obj.level(k).coordinates;
        elements = obj.level(k).elements;
        mesh = Mesh(coordinates, elements, {});
        plot(mesh, 'FaceAlpha', 0.0);
        pause(opt.timePerLevel);
    end
end