% markDoerflerQuick Employ Doerfler marking with minimal cardinality by 
% quick-select-type algorithm of average complexit O(N).
%
%   [marked, total] = markDoerflerSorting(indicators, theta) marks elements
%       based on element-wise indicators and marking parameter theta.

function [marked, total] = markDoerflerQuick(indicators, theta)

    arguments
        indicators (:,1) double
        theta (1,1) double {mustBeInRange(theta, 0, 1)}
    end

    % Compute threshold for break condition
    total = sum(indicators);
    threshold = theta * total;

    % Initialize permutation of values in indicators
    lower = 1;  % first invalid entry of permutation
    upper = length(indicators);  % final invalid entry of permutation
    permutation = nan(size(indicators));
    activePermutation = lower:upper;

    % Iterative realization of the quickMark algorithm
    minimalSetMarked = false;
    while ~minimalSetMarked
        % Abbreviate active indices
        indices = 1:length(activePermutation);

        % Compute median (MATLAB built-in with linear complexity)
        m = matlab.internal.math.columnmedian(indicators(activePermutation), false);

        % Choose exact median or larger one of the two values closest to the median
        pivot = find(abs(abs(indicators(activePermutation) - m) - min(abs(indicators(activePermutation) - m))) < 1e-12, 1);
        etaPivot = indicators(activePermutation(pivot));

        % Compare all active indices to the pivot value (median)
        IdxGreater = indicators(activePermutation) > etaPivot;
        IdxEqual = indicators(activePermutation) == etaPivot;
        IdxSmaller = ~IdxGreater & ~IdxEqual;

        greater = nnz(IdxGreater);
        smaller = greater + nnz(IdxEqual);

        % Update permutation by overwriting active indices
        activePermutation = activePermutation([indices(IdxGreater), ...
                                               indices(IdxEqual), ...
                                               indices(IdxSmaller)]);

        % Compute the sum of the selected values
        selectedSum = sum(indicators(activePermutation(1:greater)));
        selectedSumWithPivots = selectedSum + nEqual * etaPivot;

        % Check if the threshold is met
        if selectedSum >= threshold
            % Threshold is reached, but set is not minimal
            % Store accepted part of permutation
            permutation(upper-nSmaller-nEqual+1:upper) = activePermutation(greater+1:end);
            upper = upper - nSmaller - nEqual;
            % Restrict active part of permutation to accepted values
            activePermutation = activePermutation(1:greater);
        elseif selectedSumWithPivots >= threshold
            % Threshold is reached and set is minimal
            % Store remaining part of active permutation
            permutation(lower:upper) = activePermutation;
            % Determine number of finally selected contributions
            nSelected = lower - 1 + greater + ceil((threshold - selectedSum) / etaPivot);
            % Break iteration
            minimalSetMarked = true;
        else
            % Threshold is not reached
            % Store all large contributions
            permutation(lower:lower+nGreater+nEqual-1) = activePermutation(1:smaller-1);
            lower = lower+nGreater+nEqual;
            % Restrict active part of permutation to accepted values
            activePermutation = activePermutation(smaller:end);
            % Update new threshold
            threshold = threshold - selectedSumWithPivots;
        end
    end

    % Select the marked indices
    marked = reshape(permutation(1:nSelected), [], 1);
end
