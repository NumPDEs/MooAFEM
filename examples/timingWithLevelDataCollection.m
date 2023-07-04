% ******************************************************************************
% Example of timing with LevelDataCollection
% ******************************************************************************


%% parameters
nRuns = 10;

%% run time measurements
leveldatacollection = ...
    TimeIt('debugTiming', nRuns, 'storingLevelOrientedData', false);

%% print statistical analysis of timing results
leveldatacollection.printStatistics();

%% plot results
figure();
leveldatacollection.plotStatistics('ndof');

%% save results to comma-separated file
leveldatacollection.saveToTable();
