function leveldata = TimeQuickMark(nLevel)

    N = cumprod(repmat(4, nLevel, 1));

    eta = rand(N(end), 1);
    theta = 0.5;

    leveldata = LevelData();

    for j = 1:nLevel
        tic;
        quickMark(eta(1:N(j)), theta);
        runtimeQuick = toc;
        tic;
        quickMarkOptimized(eta(1:N(j)), theta);
        runtimeInPlace = toc;
        tic;
        markDoerflerSorting(eta(1:N(j)), theta);
        runtimeSort = toc;
        leveldata.append("jLevel", j, "N", N(j));
        leveldata.setTime(leveldata.nLevel, ...
                           "runtimeQuick", runtimeQuick, ...
                           "runtimeInPlace", runtimeInPlace, ...
                           "runtimeSort", runtimeSort);
        leveldata.printLevel();
    end


end