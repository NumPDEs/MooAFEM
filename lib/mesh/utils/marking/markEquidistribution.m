% markEquidistribution Employ the equidistribution marking strategy.
%
%   [marked, total] = markEquidistribution(indicators, theta) marks elements
%       where the element-wise indicators are >= (1-theta) * the mean of the
%       indicators.

function [marked, total] = markEquidistribution(indicators, theta)

arguments
    indicators (:,1) double
    theta (1,1) double {mustBeInRange(theta, 0, 1)}
end

total = sum(indicators);
n = numel(indicators);
marked = find(indicators >= (1-theta)*total/n);

end