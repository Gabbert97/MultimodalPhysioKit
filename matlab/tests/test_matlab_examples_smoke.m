function tests = test_matlab_examples_smoke
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
testCase.TestData.RepositoryRoot = repositoryRoot;
testCase.TestData.ExampleDirectory = fullfile( ...
    repositoryRoot,'examples','matlab');
testCase.TestData.OriginalPath = path;
addpath(genpath(fullfile(repositoryRoot,'matlab')))
end

function teardownOnce(testCase)
path(testCase.TestData.OriginalPath)
end

function setup(testCase)
testCase.TestData.DirectoryContents = localDirectoryContents(pwd);
end

function testTemperatureExample(testCase)
plotCleanup = localPrepareForPlots(); %#ok<NASGU>
existingFigures = findall(groot,'Type','figure');
run(fullfile(testCase.TestData.ExampleDirectory,'temperature_example.m'))

verifyClass(testCase,temperatureResults,'containers.Map')
verifyTrue(testCase,all(isKey(temperatureResults, ...
    {'segment_1','segment_2'})))
verifyClass(testCase,temperatureResults('segment_1'), ...
    'mpk.ProcessingResult')
localVerifyPlottingAndFiles(testCase,existingFigures)
end

function testRespirationExample(testCase)
plotCleanup = localPrepareForPlots(); %#ok<NASGU>
existingFigures = findall(groot,'Type','figure');
run(fullfile(testCase.TestData.ExampleDirectory,'respiration_example.m'))

verifyClass(testCase,respirationResults,'containers.Map')
verifyTrue(testCase,all(isKey(respirationResults, ...
    {'segment_1','segment_2'})))
verifyClass(testCase,respirationResults('segment_1'), ...
    'mpk.ProcessingResult')
localVerifyPlottingAndFiles(testCase,existingFigures)
end

function testEdaExampleWhenDependenciesAreAvailable(testCase)
assumeTrue(testCase,exist('cvxEDA','file') == 2, ...
    'The optional real cvxEDA dependency is unavailable.')
assumeTrue(testCase,exist('quadprog','file') == 2, ...
    'The optimization dependency required by cvxEDA is unavailable.')

plotCleanup = localPrepareForPlots(); %#ok<NASGU>
existingFigures = findall(groot,'Type','figure');
run(fullfile(testCase.TestData.ExampleDirectory,'eda_example.m'))

verifyClass(testCase,edaResults,'containers.Map')
verifyTrue(testCase,all(isKey(edaResults, ...
    {'segment_1','segment_2'})))
verifyClass(testCase,edaResults('segment_1'),'mpk.ProcessingResult')
localVerifyPlottingAndFiles(testCase,existingFigures)
end

function testEcgExample(testCase)
plotCleanup = localPrepareForPlots(); %#ok<NASGU>
existingFigures = findall(groot,'Type','figure');
run(fullfile(testCase.TestData.ExampleDirectory,'ecg_example.m'))

verifyClass(testCase,ecgResults,'containers.Map')
verifyTrue(testCase,all(isKey(ecgResults, ...
    {'segment_1','segment_2'})))
verifyClass(testCase,ecgResults('segment_1'),'mpk.ProcessingResult')
localVerifyPlottingAndFiles(testCase,existingFigures)
end

function testFnirsExample(testCase)
plotCleanup = localPrepareForPlots(); %#ok<NASGU>
existingFigures = findall(groot,'Type','figure');
run(fullfile(testCase.TestData.ExampleDirectory,'fnirs_example.m'))

verifyClass(testCase,fnirsResults,'containers.Map')
verifyTrue(testCase,all(isKey(fnirsResults,{'task_1','task_2'})))
taskResult = fnirsResults('task_1');
verifyClass(testCase,taskResult,'mpk.ProcessingResult')
verifyEqual(testCase,taskResult.Metadata.baseline_mode,'provided')
localVerifyPlottingAndFiles(testCase,existingFigures)
end

function testMultimodalExample(testCase)
plotCleanup = localPrepareForPlots(); %#ok<NASGU>
existingFigures = findall(groot,'Type','figure');
run(fullfile(testCase.TestData.ExampleDirectory,'multimodal_example.m'))

verifyTrue(testCase,isstruct(results))
verifyClass(testCase,results.Temperature,'containers.Map')
verifyClass(testCase,results.Respiration,'containers.Map')
verifyClass(testCase,results.ECG,'containers.Map')
verifyClass(testCase,results.FNIRS,'containers.Map')
verifyClass(testCase,results.Temperature('baseline'), ...
    'mpk.ProcessingResult')
verifyClass(testCase,results.Respiration('task_1'), ...
    'mpk.ProcessingResult')
verifyClass(testCase,results.ECG('task_2'), ...
    'mpk.ProcessingResult')
verifyTrue(testCase,all(isKey(results.FNIRS,{'task_1','task_2'})))
verifyEqual(testCase, ...
    results.FNIRS('task_1').Metadata.baseline_mode, ...
    'provided')
localVerifyPlottingAndFiles(testCase,existingFigures)
end

function contents = localDirectoryContents(directory)
listing = dir(directory);
contents = sort({listing.name});
end

function cleanup = localPrepareForPlots()
previousVisibility = get(groot,'DefaultFigureVisible');
set(groot,'DefaultFigureVisible','off')
existingFigures = findall(groot,'Type','figure');
cleanup = onCleanup(@() localRestorePlots( ...
    previousVisibility,existingFigures));
end

function localVerifyPlottingAndFiles(testCase,existingFigures)
newFigures = setdiff( ...
    findall(groot,'Type','figure'), ...
    existingFigures);
verifyNotEmpty(testCase,newFigures)
verifyEqual(testCase,localDirectoryContents(pwd), ...
    testCase.TestData.DirectoryContents)
end

function localRestorePlots(previousVisibility,existingFigures)
set(groot,'DefaultFigureVisible',previousVisibility)
currentFigures = findall(groot,'Type','figure');
newFigures = setdiff(currentFigures,existingFigures);
if ~isempty(newFigures)
    close(newFigures)
end
end
