% *************************************************************************
% Setup file. Adds all subfolders to path and compiles mex-functions.
% *************************************************************************

%% create + move to directory for compiled mex-files
[rootDirectory, ~, ~] = fileparts(mfilename('fullpath'));
mexDirectory = fullfile(rootDirectory, 'bin', 'mex');
if ~isfolder(mexDirectory)
   mkdir(mexDirectory)
end
cd(mexDirectory)

%% find C-files and compile if not already compiled
cFiles = dir(fullfile(rootDirectory, '**', '*.c'));
try
    for i = 1:numel(cFiles)
        [~, filename, ~] = fileparts(cFiles(i).name);
        if ~isfile(strcat(filename, '.', mexext))
            mex(fullfile(cFiles(i).folder, cFiles(i).name));
        end
    end
catch me
    warning('MEX compilation not sucessful!\nError ID: %s\nError message:\n%s', ...
        me.identifier, me.message)
end

%% add folders and subfolders to path, but ignore .git folder
cd(rootDirectory)
addpath(genpath(rootDirectory), '-begin');
rmpath(genpath(fullfile(rootDirectory, '.git')));