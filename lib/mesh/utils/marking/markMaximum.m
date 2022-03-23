% markMaximum Employ the maximum marking strategy.
%
%   [marked, total] = markMaximum(indicators, theta) marks all elements where
%       the element-wise indicators are >= (1-theta) * the maximal indicator.

function [marked, total] = markMaximum(indicators, theta)

arguments
    indicators (:,1) double
    theta (1,1) double {mustBeInRange(theta, 0, 1)}
end

total = sum(indicators);
maxIndicator = max(indicators);
marked = find(indicators >= (1-theta)*maxIndicator);

end