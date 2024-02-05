function marked = quickMarkOptimized(eta, theta)

    % Compute threshold for break condition
    threshold = theta * sum(eta);

    % Initialize permutation of values in eta
    lower = 1;  % first invalid entry of permutation
    upper = length(eta);  % final invalid entry of permutation
    permutation = nan(size(eta));
    activePermutation = lower:upper;

    % Iterative realization of the quickMark algorithm
    minimalSetMarked = false;
    while ~minimalSetMarked
        % Abbreviate active indices
        indices = 1:length(activePermutation);

        % Compute median (MATLAB built-in with linear complexity)
        m = matlab.internal.math.columnmedian(eta(activePermutation), false);

        % Choose exact median or larger one of the two values closest to the median
        pivot = find(abs(abs(eta(activePermutation) - m) - min(abs(eta(activePermutation) - m))) < 1e-12, 1);
        etaPivot = eta(activePermutation(pivot));

        % Compare all active indices to the pivot value (median)
        IdxGreater = eta(activePermutation) > etaPivot;
        IdxEqual = eta(activePermutation) == etaPivot;
        IdxSmaller = ~IdxGreater & ~IdxEqual;

        greater = nnz(IdxGreater);
        smaller = greater + nnz(IdxEqual);

        % Update permutation by overwriting active indices
        activePermutation = activePermutation([indices(IdxGreater), ...
                                               indices(IdxEqual), ...
                                               indices(IdxSmaller)]);

        % Compute the sum of the selected values
        selectedSum = sum(eta(activePermutation(1:greater)));
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




