function marked = markDoerflerQuick(eta, theta)

    threshold = theta * sum(eta);

    lower = 1;
    upper = length(eta);

    permutation = 1:length(eta);

    [permutation, nSelected] = quickMark(eta, permutation, lower, upper, threshold);
    marked = reshape(permutation(1:nSelected), [], 1);

end