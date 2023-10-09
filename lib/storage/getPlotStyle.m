function [colorOrder, marker] = getPlotStyle()
%%GETPLOTSTYLE creates lists of default colours and markers for plotting
%   [colorOrder, marker] = GETPLOTSTYLE()

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


    % Generate a sequence of colours for the graphs
    % Matlab default
    colorOrder = [0      0.4470 0.7410;
                  0.8500 0.3250 0.0980;
                  0.9290 0.6940 0.1250;
                  0.4940 0.1840 0.5560;
                  0.4660 0.6740 0.1880;
                  0.3010 0.7450 0.9330;
                  0.6350 0.0780 0.1840];

    % ALTERNATIVE: TU corporate design colours
    % colorOrder = [0      0.4    0.6;
    %               0.3922 0.3882 0.3882;
    %               0      0.4941 0.4431;
    %               0.7294 0.2745 0.5098;
    %               0.8824 0.5373 0.1333];

    % ALTERNATIVE: variable length colour sequence
    % nColor = max(ceil(nVariable^(1/3)), 2);
    nColor = 3;
    [red, green, blue] = meshgrid(linspace(0, 1, nColor));
    % colorOrder = [red(:), green(:), blue(:)];
    % Append
    colorOrder = [colorOrder; red(:), green(:), blue(:)];
    % Remove white as line colour
    colorOrder = colorOrder(1:end-1,:);

    % Define sequence of markers
    marker = {'o', '+', 'x', 'square', 'diamond', '*', '^', 'v', ...
              'pentagram', 'hexagram', '<', '>', '.', '_', '|'};
end
