function ax = plotAbsolute(obj, xVariable, varargin)
%%PLOTABSOLUTE plots various scalar absolute variables (specified in
%varargin) with respect to the variable xVariable, returns handle to axis
%object
%   ax = PLOTABSOLUTE(obj, xVariable, varargin)
%
%   See also LevelData/plot, LevelData/plotTime

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
    if nargin >= 3
        variableName = setdiff(varargin, xVariable);
    else
        variableName = setdiff(obj.label, xVariable);
    end

    % Plot only scalar variables that do belong to absolute variables
    variableName = intersect(variableName, obj.scalarVariable);
    variableName = intersect(variableName, obj.absoluteVariable);

    % Creates semilogx plot 
    ax = obj.plotLevel('semilogx', xVariable, variableName);

    % Add title
    title(ax, 'Value plot');

    % Add axes labels
    xlabel(ax, xVariable);
    ylabel(ax, 'value');

    % Update legend
    configureLegend(ax, 'northeast');
end