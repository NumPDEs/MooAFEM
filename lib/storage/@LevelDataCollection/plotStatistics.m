function ax = plotStatistics(obj, xVariable, varargin)
%%PLOTSTATISTICS plots with statistical information (mean value with min
%and max value) of scalar time variables (specified in varargin) with
%respect to the variable xVariable, returns handle to axis object
%   PLOTSTATISTICS(obj, xVariable, varargin)

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


    % Choose time variables for plotting
    if nargin < 3
        varargin = obj.timeVariable;
    else
        varargin = intersect(varargin, obj.timeVariable);
    end

    % Specify separator for table
    obj.separator = '  ';

    % Extract data for x-axis
    xData = obj.item{1}.get(':', xVariable);

    % Extract data from each item and on each level
    data = obj.get(':', varargin{:});

    % Create plot
    ax = plotData(xData, varargin, data);

    % Add title
    title(ax, 'Runtime plot');

    % Add axes labels
    xlabel(ax, xVariable);
    ylabel(ax, 'runtime (in sec.)');

    % Update legend
    configureLegend(ax, 'northwest');
end


function ax = plotData(xData, variableName, data)
%%PLOTDATA auxiliary private function for creation of plots, returns handle
%to axis object
%   ax = PLOTDATA(xData, variableName, data)

    % Create handle to currently active axis object
    ax = gca;

    % Create unified list for formatting plot lines
    [COLOURORDER, MARKER] = getPlotStyle();

    % Iterate over given variables
    for j = 1:length(variableName)
        % Compute statistics
        yMean = mean(data{j}, 2);
        yMin = min(data{j}, [], 2);
        yMax = max(data{j}, [], 2);
        % Extract label for legend from dictionary
        variableLabel = variableName{j};
        % Create plot
		errorbar(ax, xData, yMean, ...
				 yMean - yMin, yMax - yMean, ...
				 'LineStyle', '-', 'Marker', '.', ...
				 'Color', COLOURORDER(j,:), ...
				 'DisplayName', variableLabel);
        % Use double-logarithmic scales
        set(ax, 'XScale', 'log', 'YScale', 'log');
        % Add new line into the current figure when calling plotConvergence again
        hold(ax, 'on');
    end
end
