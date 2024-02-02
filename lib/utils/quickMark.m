function [permutation, nSelected] = quickMark(eta, permutation, lower, upper, threshold)

    % Compute median as pivot element
    [~, idx] = quickMedian(eta(permutation(lower:upper)));
    pivot = lower - 1 + idx;
    med = eta(permutation(pivot));

    % Update permutation
    [permutation, greater, smaller] = partition(eta, permutation, lower, upper, pivot);

    selectedSum = sum(eta(permutation(lower:greater)));
    selectedSumWithPivots = selectedSum + (smaller - greater - 1) * med;

    if selectedSum >= threshold
        [permutation, nSelected] = quickMark(eta, permutation, lower, greater, threshold);
    elseif selectedSumWithPivots >= threshold
        nSelected = greater + ceil((threshold - selectedSum) / med);
    else
        newThreshold = threshold - selectedSumWithPivots;
        [permutation, nSelected] = quickMark(eta, permutation, smaller, upper, newThreshold);
    end
end


function [m, position] = quickMedian(x, varargin)
    if nargin < 2
        varargin = {};
    end
    % Compute median (MATLAB built-in with linear complexity)
    m = median(x, varargin{:});
    % Choose two values closest to the median
    position = find(abs(abs(x - m) - min(abs(x - m))) < 1e-12);
    xmedian = x(position);
    % Choose larger one of them
    [~, idx] = min(xmedian);
    position = position(idx);
end


function [permutation, greater, smaller] = partition(eta, permutation, lower, upper, pivot)

    % Initialize all possibly modified indices
    indices = lower:upper;

    etaPivot = eta(permutation(pivot));

    % Select indices for smaller values
    permutationGreater = permutation(indices(eta(permutation(indices)) > etaPivot));
    permutationSmaller = permutation(indices(eta(permutation(indices)) < etaPivot));
    permutationEqual = permutation(indices(eta(permutation(indices)) == etaPivot));

    smaller = upper + 1 - length(permutationSmaller);
    greater = lower - 1 + length(permutationGreater);

    % Return the updated permutation
    permutation = [permutation(1:lower-1), ...
                   permutationGreater, ...
                   permutationEqual, ...
                   permutationSmaller, ...
                   permutation(upper+1:end)];

end
