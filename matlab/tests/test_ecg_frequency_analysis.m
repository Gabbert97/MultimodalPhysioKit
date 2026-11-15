function tests = test_ecg_frequency_analysis
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

function testPsdKeepsRealTimeAcrossRejectedIntervals(testCase)
% Respiratory sinus arrhythmia at 0.3125 Hz (a PSD grid bin). Rejected
% intervals must leave a gap on the real timeline; removing every 10th
% interval from time compresses the tachogram by ~10% and moves the RSA
% peak two bins up (0.3438 Hz with the old cumulative timeline).
[ibi,starts] = localRsaBeats(300,0.3125);
keep = true(size(ibi));
keep(16:10:numel(ibi)-15) = false;
[~,~,~,~,~,~,frequencyData] = frequency_analysis_ecg( ...
    ibi(keep),sum(keep),1,2,false,0,'test',false,[],{starts(keep)});
[~,peakIndex] = max(frequencyData.PSD{1});
verifyEqual(testCase,frequencyData.Frequencies{1}(peakIndex),0.3125, ...
    'AbsTol',0.0078)
end

function testPsdWithoutRejectionsMatchesCumulativeTimeline(testCase)
[ibi,starts] = localRsaBeats(120,0.25);
[plfA,~,phfA,~,lfhfA] = frequency_analysis_ecg( ...
    ibi,numel(ibi),1,2,false,0,'test',false,[]);
[plfB,~,phfB,~,lfhfB] = frequency_analysis_ecg( ...
    ibi,numel(ibi),1,2,false,0,'test',false,[],{starts + 7});
verifyEqual(testCase,[plfB phfB lfhfB],[plfA phfA lfhfA],'AbsTol',1e-9)
end

function [ibi,starts] = localRsaBeats(durationSeconds,rsaFrequency)
beatTimes = 0;
while beatTimes(end) < durationSeconds
    last = beatTimes(end);
    beatTimes(end+1,1) = last + 0.9 + 0.05*sin(2*pi*rsaFrequency*last); %#ok<AGROW>
end
ibi = diff(beatTimes);
starts = beatTimes(1:end-1);
end
