function ax = plotLevel(obj, category, xVariable, yVariable)
%%PLOTLEVEL auxiliary private function for creation of plots for usage in
%LevelData/plot, LevelData/plotTime, and LevelData/plotAbsolute, uses the
%specified DataCategory for the plot of variables in yVariable with
%respect to xVariable, returns handle to axis object
%   ax = PLOTLEVEL(obj, category, xVariable, yVariable, ...)

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
        category DataCategory
        xVariable {mustBeTextScalar}
        yVariable (1,:) cell
    end

    % Extract variables with correct category, type, and shape for y-axis
    idx = (obj.category(yVariable) == category) ...
        & ([obj.type(yVariable).rawType] == RawType.FLOAT) ...
        & [obj.type(yVariable).isScalar];
    yVariable = yVariable(idx); 

    % Create handle to currently active axis object
    ax = gca;

    % Extract variable for x-axis
    xValue = obj.get(1:obj.nLevel, xVariable);

    % Create unified list for formatting plot lines
    [COLOURORDER, MARKER] = getPlotStyle();

    % Iterate over given variables
    for j = 1:length(yVariable)
        % Extract value for y-axis
        yValue = obj.get(1:obj.nLevel, yVariable{j});

        % Extract label for legend from dictionary
        if obj.dictionary.isKey(yVariable{j})
            variableLabel = obj.dictionary(yVariable{j});
        else
            variableLabel = yVariable{j};
        end

        % Create plot
        category.plotFunction( ...
                ax, xValue, yValue, '-', ...
                'Marker', MARKER{mod(j, length(MARKER))}, ...
                'Color', COLOURORDER(j, :), ...
                'LineWidth', 1.5, ...
                'DisplayName', variableLabel);
        % Add new line into the current figure when calling plotConvergence again
        hold(ax, 'on');
    end

    % Set plot appearance
    title(ax, category.title);
    xlabel(ax, xVariable);
    ylabel(ax, category.yLabel);
    configureLegend(ax, category.legendLocation);
end
