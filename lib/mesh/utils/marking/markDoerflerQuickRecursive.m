% markDoerflerQuick Employ Doerfler marking with minimal cardinality by 
% quick-select-type algorithm of average complexit O(N).
%
%   [marked, total] = markDoerflerSorting(indicators, theta) marks elements
%       based on element-wise indicators and marking parameter theta.

function marked = markDoerflerQuickRecursive(indicators, theta)

    arguments
        indicators (:,1) double
        theta (1,1) double {mustBeInRange(theta, 0, 1)}
    end
    
    % Compute threshold for break condition
    total = sum(indicators);
    threshold = theta * total;

    % Initialize permutation
    lower = 1;
    upper = length(indicators);
    permutation = 1:length(indicators);

    % Recursively update the permutation to gather the largest contributions
    % in the beginning
    [permutation, nSelected] = recursiveCall(indicators, permutation, lower, upper, threshold);

    % Return the marked elements
    marked = reshape(permutation(1:nSelected), [], 1);
end


function [permutation, nSelected] = recursiveCall(indicators, permutation, lower, upper, threshold)

    % Compute median (MATLAB built-in with linear complexity)
    % m = median(x, varargin{:});
    m = matlab.internal.math.columnmedian(indicators(permutation(lower:upper)), false);

    % Choose larger one of two values closest to the median
    % and determine its position
    position = find(abs(abs(indicators(permutation(lower:upper)) - m) ...
                            - min(abs(indicators(permutation(lower:upper)) - m))) < 1e-12);
    [~, idx] = min(x(position));
    position = position(idx);

    % Choose pivot element of the overall permuation
    pivot = lower - 1 + position;
    etaPivot = indicators(permutation(pivot));

    % Update permutation
    [permutation, greater, smaller] = partition(indicators, permutation, lower, upper, etaPivot);

    % Compute new sum of selected elements
    selectedSum = sum(indicators(permutation(lower:greater)));
    selectedSumWithPivots = selectedSum + (smaller - greater - 1) * etaPivot;

    % Check for break condition
    if selectedSum >= threshold
        % Threshold is reached, but set is not minimal
        % Store accepted part of permutation
        [permutation, nSelected] = recursiveCall(indicators, permutation, lower, greater, threshold);
    elseif selectedSumWithPivots >= threshold
        % Threshold is reached and set is minimal
        nSelected = greater + ceil((threshold - selectedSum) / etaPivot);
    else
        % Threshold is not reached
        newThreshold = threshold - selectedSumWithPivots;
        [permutation, nSelected] = recursiveCall(indicators, permutation, smaller, upper, newThreshold);
    end
end


function [permutation, greater, smaller] = partition(indicators, permutation, lower, upper, etaPivot)

    % Initialize all possibly modified indices
    indices = lower:upper;

    % Determine indices for smaller and larger values
    permutationGreater = permutation(indices(indicators(permutation(indices)) > etaPivot));
    permutationSmaller = permutation(indices(indicators(permutation(indices)) < etaPivot));
    permutationEqual = permutation(indices(indicators(permutation(indices)) == etaPivot));

    smaller = upper + 1 - length(permutationSmaller);
    greater = lower - 1 + length(permutationGreater);

    % Return the updated permutation
    permutation = [permutation(1:lower-1), ...
                   permutationGreater, ...
                   permutationEqual, ...
                   permutationSmaller, ...
                   permutation(upper+1:end)];
end
