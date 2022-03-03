% markDoerflerBinning Employ Doerfler marking by binning.
%
%   [marked, total] = markDoerflerBinning(indicators, theta) marks elements
%       based on element-wise indicators and marking parameter theta. Returns a
%       set of marked elements with up to twice minimal cardinality.
%
%   [marked, total] = markDoerflerBinning(indicators, theta, C) returns a set of
%       marked elements with up to C > 1 times minimal cardinality.

function [marked, total] = markDoerflerBinning(indicators, theta, C)

arguments
    indicators (:,1) double
    theta (1,1) double {mustBeInRange(theta, 0, 1)}
    C (1,1) double {mustBeGreaterThan(C, 1)} = 2
end

% preparation
N = length(indicators);
M = max(indicators);
total = sum(indicators);
K = floor(log(N*M/((1-theta)*total)) / log(C));

% find bins and rearrange values
bins = findgroups(min(floor(-log(indicators/M) / log(C)) + 1, K+2));
permutation = cell2mat(splitapply(@(x) {x}, (1:N)', bins));

% find subset such that Doerfler marking is fullfilled
sumIndicators = cumsum(indicators(permutation));
L = find(theta*total <= sumIndicators, 1, 'first');
marked = permutation(1:L);

end