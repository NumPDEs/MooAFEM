function marked = quickMark(eta, theta)

    % Compute threshold for break condition
    threshold = theta * sum(eta);

    % Initialize permutation of values in eta
    lower = 1;
    upper = length(eta);
    permutation = lower:upper;

    % Iterative realization of the quickMark algorithm
    minimalSetMarked = false;
    while ~minimalSetMarked
        % Abbreviate active indices
        indices = lower:upper;

        % Compute median (MATLAB built-in with linear complexity)
        m = matlab.internal.math.columnmedian(eta(permutation(lower:upper)), false);

        % Choose exact median or larger one of the two values closest to the median
        position = find(abs(abs(eta(permutation(lower:upper)) - m) - min(abs(eta(permutation(lower:upper)) - m))) < 1e-12, 1);
        pivot = lower - 1 + position;
        etaPivot = eta(permutation(pivot));

        % Compare all active indices to the pivot value (median)
        IdxGreater = eta(permutation(lower:upper)) > etaPivot;
        IdxEqual = eta(permutation(lower:upper)) == etaPivot;
        IdxSmaller = ~IdxGreater & ~IdxEqual;

        % Compute boundary indices for the partition
        smaller = upper + 1 - nnz(IdxSmaller);
        greater = lower - 1 + nnz(IdxGreater);

        % Update permutation by overwriting active indices
        permutation(lower:upper) = permutation([indices(IdxGreater), ...
                                                indices(IdxEqual), ...
                                                indices(IdxSmaller)]);

        % Compute the sum of the selected values
        selectedSum = sum(eta(permutation(lower:greater)));
        selectedSumWithPivots = selectedSum + (smaller - greater - 1) * etaPivot;

        % Check if the threshold is reached
        if selectedSum >= threshold
            % Threshold is reached, but set is not minimal
            upper = greater;
        elseif selectedSumWithPivots >= threshold
            % Threshold is reached and set is minimal
            nSelected = greater + ceil((threshold - selectedSum) / etaPivot);
            minimalSetMarked = true;
        else
            % Threshold is not reached
            threshold = threshold - selectedSumWithPivots;
            lower = smaller;
        end
    end

    % Select the marked indices
    marked = reshape(permutation(1:nSelected), [], 1);
end

