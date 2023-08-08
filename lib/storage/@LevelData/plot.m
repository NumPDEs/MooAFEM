function ax = plot(obj, xVariable, yVariable)
%%PLOT plots various scalar variables (specified in yVariable) with respect
%to the variable xVariable, returns handle to axis object
%   ax = PLOT(obj, xVariable, yVariable, ...)
%
%   See also LevelData/plotTime, LevelData/plotAbsolute

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
        xVariable {mustBeTextScalar}
    end

    arguments (Repeating)
        yVariable {mustBeTextScalar}
    end

    if isempty(yVariable)
        yVariable = setdiff(obj.label, xVariable);
    end

    % TODO: introduce plot specification to allow for different plot types (time, absolute, level)

    % Plot only scalar variables that do not belong to time or absolute
    % variables
    yVariable = intersect(yVariable, obj.scalarVariable);
    yVariable = setdiff(yVariable, obj.absoluteVariable);
    yVariable = setdiff(yVariable, obj.timeVariable);

    % Creates double logarithmic splot 
    ax = obj.plotLevel(@loglog, xVariable, yVariable);

    % Add title
    title(ax, 'Convergence history plot');

    % Add axes labels
    xlabel(ax, xVariable);
    ylabel(ax, 'error');

    % Update legend
    configureLegend(ax, 'northeast');
end