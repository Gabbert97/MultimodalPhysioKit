%RUN_MATLAB_TESTS Run the complete MultimodalPhysioKit MATLAB test suite.

scriptPath = mfilename('fullpath');
testsDirectory = fileparts(scriptPath);
matlabRoot = fileparts(testsDirectory);

addpath(genpath(matlabRoot))

results = runtests( ...
    testsDirectory, ...
    'IncludeSubfolders',true);

disp(table(results))
assertSuccess(results)
