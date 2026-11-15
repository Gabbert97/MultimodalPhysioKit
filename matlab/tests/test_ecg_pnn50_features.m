function tests = test_ecg_pnn50_features
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
testCase.TestData.OriginalPath = path;
addpath(genpath(fullfile(repositoryRoot,'matlab')))
end

function teardownOnce(testCase)
path(testCase.TestData.OriginalPath)
end

function testComparesOnlyAdjacentIbisAfterRejectedInterval(testCase)
% Interval 4 was rejected: runs [0.8 0.8 0.8] and [1.0 1.0 1.1 1.2].
% Evaluated differences: [0 0] and [0 0.1 0.1]; 2 of 5 exceed 50 ms.
flags = [1 1 1 0 1 1 1 1];
ibi = [0.8 0.8 0.8 1.0 1.0 1.1 1.2];
[pnn50,pnn50Max] = ecg_pnn50_features(ibi,flags);
verifyEqual(testCase,pnn50,2/5,'AbsTol',1e-12)
verifyEqual(testCase,pnn50Max,2/(120*numel(ibi)/sum(ibi)),'AbsTol',1e-12)
end

function testCountsIsolatedDifferencesAbove50Ms(testCase)
% Task Force (1996) pNN50: alternating 0.8/0.9 s gives 5 of 9 differences.
ibi = [0.8 0.9 0.9 0.8 0.8 0.9 0.9 0.8 0.8 0.9];
[pnn50,pnn50Max] = ecg_pnn50_features(ibi,ones(size(ibi)));
verifyEqual(testCase,pnn50,5/9,'AbsTol',1e-12)
verifyEqual(testCase,pnn50Max,1/(120*numel(ibi)/sum(ibi)),'AbsTol',1e-12)
end

function testExcludesDifferencesOf50MsOrLess(testCase)
ibi = [0.80 0.84 0.80 0.84 0.80 0.84];
[pnn50,pnn50Max] = ecg_pnn50_features(ibi,ones(size(ibi)));
verifyEqual(testCase,pnn50,0)
verifyEqual(testCase,pnn50Max,0)
end
