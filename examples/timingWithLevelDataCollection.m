% ******************************************************************************
% Example of timing with LevelDataCollection
% ******************************************************************************


%% parameters
nRuns = 10;

%% run time measurements
storeResults = false;
leveldatacollection = ...
    TimeIt('debugTiming', nRuns, storeResults, 'storingLevelOrientedData', false, false);

%% print statistical analysis of timing results
leveldatacollection.printStatistics();

%% plot results
figure();
leveldatacollection.plotStatistics('ndof');

%% save results to comma-separated file
% leveldatacollection.saveToTable();
