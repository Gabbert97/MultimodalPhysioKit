function tests = test_utils_create_sliding_phase_intervals
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

function testCompleteNonoverlappingWindows(testCase)
[intervals,labels] = mpk.utils.createSlidingPhaseIntervals(0,8,4,0);
verifyEqual(testCase,intervals,[0 4;4 8])
verifyEqual(testCase,labels,{'window_1','window_2'})
verifyClass(testCase,intervals,'double')
end

function testOverlappingWindowsAndFloatingPointStarts(testCase)
[intervals,labels] = mpk.utils.createSlidingPhaseIntervals( ...
    1.25,7.25,2.5,0.5,'LabelPrefix','segment');
verifyEqual(testCase,intervals,[1.25 3.25;3.75 5.75])
verifyEqual(testCase,labels,{'segment_1','segment_2'})
end

function testTrailingWarningAndDisable(testCase)
verifyWarning(testCase,@() mpk.utils.createSlidingPhaseIntervals( ...
    0,9.5,4,0),'MultimodalPhysioKit:DroppedIncompletePhase')
verifyWarningFree(testCase,@() mpk.utils.createSlidingPhaseIntervals( ...
    0,9.5,4,0,'WarnOnDroppedTrailingWindow',false))
end

function testIntervalsConfigureRecordingDirectly(testCase)
[intervals,labels] = mpk.utils.createSlidingPhaseIntervals(0,12,4,0);
recording = mpk.Recording('sliding');
recording.setPhases(intervals,labels);
verifyEqual(testCase,recording.PhaseIntervalsSeconds,intervals)
verifyEqual(testCase,recording.PhaseLabels,labels)
end

function testAvoidsCumulativeDriftAndDoesNotExceedEnd(testCase)
startTime = 0.1;
endTime = 1.1;
windowDuration = 0.2;
overlap = 0.1;
[intervals,~] = mpk.utils.createSlidingPhaseIntervals( ...
    startTime,endTime,windowDuration,overlap, ...
    'WarnOnDroppedTrailingWindow',false);
expectedStarts = startTime + ...
    (0:(size(intervals,2)-1))*(windowDuration-overlap);
verifyEqual(testCase,intervals(1,:),expectedStarts,'AbsTol',0)
verifyLessThanOrEqual(testCase,max(intervals(2,:)),endTime)
verifyEqual(testCase,intervals(2,:) - intervals(1,:), ...
    repmat(windowDuration,1,size(intervals,2)), ...
    'AbsTol',8*eps)
end

function testValidation(testCase)
verifyError(testCase,@() mpk.utils.createSlidingPhaseIntervals( ...
    -1,10,2,0), ...
    'mpk:utils:createSlidingPhaseIntervals:InvalidStart')
verifyError(testCase,@() mpk.utils.createSlidingPhaseIntervals( ...
    0,10,2,2), ...
    'MultimodalPhysioKit:InvalidPhaseOverlap')
verifyError(testCase,@() mpk.utils.createSlidingPhaseIntervals( ...
    0,10,2,3), ...
    'MultimodalPhysioKit:InvalidPhaseOverlap')
verifyError(testCase,@() mpk.utils.createSlidingPhaseIntervals( ...
    0,10,0,0), ...
    'mpk:utils:createSlidingPhaseIntervals:InvalidDuration')
verifyError(testCase,@() mpk.utils.createSlidingPhaseIntervals( ...
    0,1,2,0), ...
    'MultimodalPhysioKit:SampleTooShort')
verifyError(testCase,@() mpk.utils.createSlidingPhaseIntervals( ...
    0,10,2,0,'MinimumWindowDurationSeconds',3), ...
    'mpk:utils:createSlidingPhaseIntervals:InvalidDuration')
end
