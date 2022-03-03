% markDoerflerSorting Employ Doerfler marking by sorting indicators.
%
%   [marked, total] = markDoerflerSorting(indicators, theta) marks elements
%       based on element-wise indicators and marking parameter theta.

function [marked, total] = markDoerflerSorting(indicators, theta)

arguments
    indicators (:,1) double
    theta (1,1) double {mustBeInRange(theta, 0, 1)}
end

[indicators, idx] = sort(indicators, 'descend');
sumIndicators = cumsum(indicators);
total = sumIndicators(end);
N = find(theta*total <= sumIndicators, 1, 'first');
marked = idx(1:N);

end