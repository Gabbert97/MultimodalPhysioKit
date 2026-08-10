function tests = test_core_recording
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
testCase.TestData.OriginalPath = path;
addpath(genpath(fullfile(repositoryRoot,'matlab')))
addpath(genpath(fullfile(repositoryRoot,'matlab','tests','helpers')))
end

function teardownOnce(testCase)
path(testCase.TestData.OriginalPath)
end

function testConstructsRecordingWithMetadata(testCase)
metadata = struct('SourceFormat','array');

recording = mpk.Recording('session_001',metadata);

verifyClass(testCase,recording,'mpk.Recording')
verifyEqual(testCase,recording.Id,'session_001')
verifyEqual(testCase,recording.Signals,struct())
verifyEqual(testCase,recording.Metadata,metadata)
verifyEqual(testCase,recording.PhaseIntervalsSeconds,zeros(2,0))
verifyEqual(testCase,recording.PhaseLabels,cell(1,0))
end

function testUsesDefaultConstructorValues(testCase)
recording = mpk.Recording();

verifyEqual(testCase,recording.Id,'')
verifyEqual(testCase,recording.Signals,struct())
verifyEqual(testCase,recording.Metadata,struct())
verifyEqual(testCase,recording.numberOfPhases(),0)
verifyEqual(testCase,recording.phaseDurationsSeconds(),zeros(1,0))
end

function testAddsGetsAndListsSignals(testCase)
recording = mpk.Recording('session');
ecg = localSignal('ecg');
temperature = localSignal('temp');

recording.addSignal('ecg lead 1',ecg);
recording.addSignal('temp',temperature);

verifyTrue(testCase,recording.hasSignal('ecg lead 1'))
verifyTrue(testCase,recording.hasSignal('temp'))
verifySameHandleOrValue(testCase,recording.getSignal('ecg lead 1'),ecg)
verifySameHandleOrValue(testCase,recording.getSignal('temp'),temperature)
verifyEqual(testCase,sort(recording.signalNames()), ...
    {'ecgLead1'; 'temp'})
end

function testAddSignalReplacesExistingName(testCase)
recording = mpk.Recording('session');
first = localSignal('first');
second = localSignal('second');

recording.addSignal('ecg',first);
recording.addSignal('ecg',second);

verifyEqual(testCase,recording.getSignal('ecg').Label,'second')
verifyEqual(testCase,recording.signalNames(),{'ecg'})
end

function testRejectsInvalidMetadata(testCase)
verifyError(testCase,@() mpk.Recording('session','private'), ...
    'mpk:Recording:InvalidMetadata')
end

function testRejectsInvalidSignal(testCase)
recording = mpk.Recording('session');

verifyError(testCase,@() recording.addSignal('ecg',[1; 2; 3]), ...
    'mpk:Recording:InvalidSignal')
end

function testMissingSignalErrors(testCase)
recording = mpk.Recording('session');

verifyFalse(testCase,recording.hasSignal('ecg'))
verifyError(testCase,@() recording.getSignal('ecg'), ...
    'mpk:Recording:MissingSignal')
end

function testSetsFloatingPointPhasesAndDefaultLabels(testCase)
recording = mpk.Recording('session');
intervals = [0.25 10.75 21.125; 5.5 20.625 30.875];

recording.setPhases(intervals);

verifyEqual(testCase,recording.PhaseIntervalsSeconds,intervals)
verifyEqual(testCase,recording.PhaseLabels, ...
    {'phase_1','phase_2','phase_3'})
verifyEqual(testCase,recording.numberOfPhases(),3)
verifyEqual(testCase,recording.phaseDurationsSeconds(), ...
    intervals(2,:) - intervals(1,:))
end

function testNormalizesAndPreservesSuppliedLabels(testCase)
recording = mpk.Recording('session');
recording.setPhases([0 10; 10 20],["second","first"]);

verifyEqual(testCase,recording.PhaseLabels,{'second','first'})
verifySize(testCase,recording.PhaseLabels,[1 2])

recording.setPhases([0 10; 10 20],{'alpha';'beta'});
verifyEqual(testCase,recording.PhaseLabels,{'alpha','beta'})

recording.setPhases([0 5; 10 8],{'outer','inner'});
verifyEqual(testCase,recording.PhaseIntervalsSeconds,[0 5; 10 8])
end

function testRejectsInvalidPhaseLabels(testCase)
recording = mpk.Recording('session');
intervals = [0 10; 10 20];

verifyError(testCase,@() recording.setPhases(intervals,{'same','same'}), ...
    'mpk:Recording:DuplicatePhaseLabels')
verifyError(testCase,@() recording.setPhases(intervals,{'','task'}), ...
    'mpk:Recording:EmptyPhaseLabel')
verifyError(testCase,@() recording.setPhases(intervals,{'only_one'}), ...
    'mpk:Recording:PhaseLabelMismatch')
verifyError(testCase,@() recording.setPhases(intervals,{'ok',2}), ...
    'mpk:Recording:InvalidPhaseLabels')
end

function testRejectsInvalidPhaseIntervals(testCase)
recording = mpk.Recording('session');

verifyError(testCase,@() recording.setPhases([0 1 2]), ...
    'mpk:Recording:InvalidPhaseIntervals')
verifyError(testCase,@() recording.setPhases(reshape(1:8,[2 2 2])), ...
    'mpk:Recording:InvalidPhaseIntervals')
verifyError(testCase,@() recording.setPhases({'0';'1'}), ...
    'mpk:Recording:InvalidPhaseIntervals')
verifyError(testCase,@() recording.setPhases([0;NaN]), ...
    'mpk:Recording:InvalidPhaseIntervals')
verifyError(testCase,@() recording.setPhases([0;Inf]), ...
    'mpk:Recording:InvalidPhaseIntervals')
verifyError(testCase,@() recording.setPhases([-1;1]), ...
    'mpk:Recording:InvalidPhaseIntervals')
verifyError(testCase,@() recording.setPhases([1;1]), ...
    'mpk:Recording:InvalidPhaseIntervals')
verifyError(testCase,@() recording.setPhases([2;1]), ...
    'mpk:Recording:InvalidPhaseIntervals')
end

function testClearPhasesRestoresEmptyState(testCase)
recording = mpk.Recording('session');
recording.setPhases([0 10; 10 20],{'first','second'});

recording.clearPhases();

verifyEqual(testCase,recording.PhaseIntervalsSeconds,zeros(2,0))
verifyEqual(testCase,recording.PhaseLabels,cell(1,0))
verifyEqual(testCase,recording.numberOfPhases(),0)
verifyEqual(testCase,recording.phaseDurationsSeconds(),zeros(1,0))
end

function testPhaseIndicesAcceptsNameAndSignal(testCase)
recording = mpk.Recording('session');
signal = mpk.Signal((1:5)',(5:0.5:7)',2,'ecg','mV');
recording.addSignal('ecg',signal);
recording.setPhases([5 6; 6 7.5],{'first','second'});

expected = [1 3; 3 6];
verifyEqual(testCase,recording.phaseIndices('ecg'),expected)
verifyEqual(testCase,recording.phaseIndices('ecg',[]),expected)
verifyEqual(testCase,recording.phaseIndices(signal),expected)
verifyEqual(testCase,recording.phaseIndices('ecg','first'),expected(:,1))
verifyEqual(testCase,recording.phaseIndices('ecg',"second"),expected(:,2))
verifyEqual(testCase,recording.phaseIndices('ecg',{'second','first'}), ...
    expected(:,[2 1]))
verifyEqual(testCase,recording.phaseIndices('ecg',["second","first"]), ...
    expected(:,[2 1]))
end

function testPhaseIndicesMapsBoundariesBetweenSamples(testCase)
recording = mpk.Recording('session');
signal = mpk.Signal((1:5)',(5:0.5:7)',2,'ecg','mV');
recording.addSignal('ecg',signal);
recording.setPhases([5.1;6.1],{'between'});

verifyEqual(testCase,recording.phaseIndices('ecg'),[2;4])
end

function testPhaseSelectionErrors(testCase)
recording = mpk.Recording('session');
signal = mpk.Signal((1:5)',(0:0.5:2)',2,'ecg','mV');
recording.addSignal('ecg',signal);

verifyError(testCase,@() recording.phaseIndices('ecg'), ...
    'mpk:Recording:NoPhases')

recording.setPhases([0;1],{'known'});
verifyError(testCase,@() recording.phaseIndices('missing'), ...
    'mpk:Recording:MissingSignal')
verifyError(testCase,@() recording.phaseIndices('ecg','missing'), ...
    'mpk:Recording:UnknownPhase')
verifyError(testCase,@() recording.phaseIndices('ecg',cell(1,0)), ...
    'mpk:Recording:InvalidPhaseSelection')
verifyError(testCase,@() recording.phaseIndices('ecg',strings(1,0)), ...
    'mpk:Recording:InvalidPhaseSelection')
verifyError(testCase,@() recording.phaseIndices('ecg',1), ...
    'mpk:Recording:InvalidPhaseSelection')
end

function testPhaseIndicesRejectsOutsideAndEmptyMappings(testCase)
recording = mpk.Recording('session');
signal = mpk.Signal((1:5)',(5:0.5:7)',2,'ecg','mV');
recording.addSignal('ecg',signal);

recording.setPhases([4.9;5.5],{'before'});
verifyError(testCase,@() recording.phaseIndices('ecg'), ...
    'mpk:Recording:PhaseOutsideSignal')

recording.setPhases([7;7.6],{'after'});
verifyError(testCase,@() recording.phaseIndices('ecg'), ...
    'mpk:Recording:PhaseOutsideSignal')

recording.setPhases([5.1;5.4],{'empty'});
verifyError(testCase,@() recording.phaseIndices('ecg'), ...
    'mpk:Recording:EmptyPhase')
end

function testPhaseIndicesValidatesSignalTime(testCase)
recording = mpk.Recording('session');
recording.setPhases([0;1],{'phase'});
signal = mpk.Signal((1:3)',(0:2)',1,'ecg','mV');

signal.Samples = [];
signal.Time = [];
verifyError(testCase,@() recording.phaseIndices(signal), ...
    'mpk:Recording:EmptySignal')

signal = mpk.Signal((1:3)',(0:2)',1,'ecg','mV');
signal.Time(2) = NaN;
verifyError(testCase,@() recording.phaseIndices(signal), ...
    'mpk:Recording:InvalidSignalTime')

signal.Time = [0;2;1];
verifyError(testCase,@() recording.phaseIndices(signal), ...
    'mpk:Recording:InvalidSignalTime')
end

function testExtractPhaseResetsTimeAndPreservesSignalProperties(testCase)
recording = mpk.Recording('session_42');
sourceMetadata = struct('Sensor','ECG','Gain',4);
source = mpk.Signal((10:10:50)',(5:0.5:7)',2, ...
    'lead_1','mV',sourceMetadata);
recording.addSignal('ecg',source);
recording.setPhases([5.5;6.5],{'segment'});
originalIntervals = recording.PhaseIntervalsSeconds;

extracted = recording.extractPhase('ecg','segment');

verifyClass(testCase,extracted,'mpk.Signal')
verifyEqual(testCase,extracted.Samples,[20;30])
verifyEqual(testCase,extracted.Time,[0;0.5])
verifyEqual(testCase,extracted.SamplingFrequency,2)
verifyEqual(testCase,extracted.Label,'lead_1')
verifyEqual(testCase,extracted.Units,'mV')
verifyEqual(testCase,extracted.Metadata.Sensor,'ECG')
verifyEqual(testCase,extracted.Metadata.Gain,4)
verifyEqual(testCase,extracted.Metadata.SourceRecording,'session_42')
verifyEqual(testCase,extracted.Metadata.PhaseLabel,'segment')
verifyEqual(testCase,extracted.Metadata.PhaseStartSeconds,5.5)
verifyEqual(testCase,extracted.Metadata.PhaseEndSeconds,6.5)
verifyTrue(testCase,extracted.Metadata.TimeReset)
verifyEqual(testCase,source.Time,(5:0.5:7)')
verifyEqual(testCase,source.Metadata,sourceMetadata)
verifyEqual(testCase,recording.PhaseIntervalsSeconds,originalIntervals)
end

function testExtractPhaseCanPreserveTimestamps(testCase)
recording = mpk.Recording('session');
recording.addSignal('ecg',mpk.Signal((1:5)',(5:0.5:7)',2, ...
    'lead_1','mV'));
recording.setPhases([5.5;6.5],{'segment'});

extracted = recording.extractPhase( ...
    'ecg','segment','ResetTime',false);

verifyEqual(testCase,extracted.Time,[5.5;6])
verifyFalse(testCase,extracted.Metadata.TimeReset)
verifyError(testCase,@() recording.extractPhase( ...
    'ecg','segment','ResetTime',1), ...
    'mpk:Recording:InvalidResetTime')
end

function testDifferentSamplingFrequenciesMapIndependently(testCase)
recording = mpk.Recording('multimodal');
ecgFs = 1000;
temperatureFs = 10;
ecgTime = (0:(30*ecgFs-1))' ./ ecgFs;
temperatureTime = (0:(30*temperatureFs-1))' ./ temperatureFs;
recording.addSignal('ecg',mpk.Signal( ...
    zeros(size(ecgTime)),ecgTime,ecgFs,'ecg','mV'));
recording.addSignal('temperature',mpk.Signal( ...
    zeros(size(temperatureTime)),temperatureTime,temperatureFs, ...
    'temperature','degC'));
recording.setPhases([10;20],{'segment'});

ecgIndices = recording.phaseIndices('ecg');
temperatureIndices = recording.phaseIndices('temperature');
ecgPhase = recording.extractPhase( ...
    'ecg','segment','ResetTime',false);
temperaturePhase = recording.extractPhase( ...
    'temperature','segment','ResetTime',false);

verifyEqual(testCase,ecgIndices,[10001;20001])
verifyEqual(testCase,temperatureIndices,[101;201])
verifyEqual(testCase,ecgPhase.count(),10000)
verifyEqual(testCase,temperaturePhase.count(),100)
verifyEqual(testCase,[ecgPhase.Time(1),ecgPhase.Time(end)], ...
    [10,19.999],'AbsTol',1e-12)
verifyEqual(testCase,[temperaturePhase.Time(1),temperaturePhase.Time(end)], ...
    [10,19.9],'AbsTol',1e-12)
end

function testRemoveTimeRangesRequiresClearedPhases(testCase)
recording = mpk.Recording('session');
recording.addSignal('ecg',mpk.Signal((1:5)',(0:4)',1,'ecg','mV'));
recording.setPhases([0;2],{'segment'});

verifyError(testCase,@() recording.removeTimeRanges([1 3]), ...
    'mpk:Recording:PhasesMustBeCleared')

recording.clearPhases();
out = recording.removeTimeRanges([1 3],'ResetTime',false);
verifyEqual(testCase,out.getSignal('ecg').Samples,[1;4;5])
verifyEqual(testCase,out.getSignal('ecg').Time,[0;3;4])
verifyEqual(testCase,out.numberOfPhases(),0)
verifyEqual(testCase,out.PhaseIntervalsSeconds,zeros(2,0))
verifyEqual(testCase,out.PhaseLabels,cell(1,0))
end

function testProcessSignalPassesOneCompleteSignalAndReturnsUnchanged(testCase)
recording = localProcessingRecording();
expected = struct('Features',struct('raw_name',42),'Marker','unchanged');
processor = FakeRecordingProcessor(expected);

result = recording.processSignal(processor,'ecg');

verifyEqual(testCase,result,expected)
verifyEqual(testCase,numel(processor.Calls),1)
verifyEqual(testCase,numel(processor.Calls{1}),1)
verifyClass(testCase,processor.Calls{1}{1},'mpk.Signal')
verifyEqual(testCase,processor.Calls{1}{1}.Samples, ...
    recording.getSignal('ecg').Samples)
end

function testProcessSignalAcceptsNameFormsAndPreservesOrder(testCase)
recording = localProcessingRecording();

processor = FakeRecordingProcessor();
recording.processSignal(processor,"temperature");
verifyEqual(testCase,processor.Calls{1}{1}.Label,'temperature')

processor = FakeRecordingProcessor();
recording.processSignal(processor,{'temperature','ecg'});
verifyEqual(testCase,processor.Calls{1}{1}.Label,'temperature')
verifyEqual(testCase,processor.Calls{1}{2}.Label,'ecg')

processor = FakeRecordingProcessor();
recording.processSignal(processor,["ecg","temperature"]);
verifyEqual(testCase,processor.Calls{1}{1}.Label,'ecg')
verifyEqual(testCase,processor.Calls{1}{2}.Label,'temperature')
end

function testProcessSignalValidatesNamesProcessorAndArguments(testCase)
recording = localProcessingRecording();
processor = FakeRecordingProcessor();

verifyError(testCase,@() recording.processSignal(processor,cell(1,0)), ...
    'mpk:Recording:InvalidSignalNames')
verifyError(testCase,@() recording.processSignal(processor,strings(1,0)), ...
    'mpk:Recording:InvalidSignalNames')
verifyError(testCase,@() recording.processSignal(processor,''), ...
    'mpk:Recording:InvalidSignalNames')
verifyError(testCase,@() recording.processSignal(processor,{'ecg',''}), ...
    'mpk:Recording:InvalidSignalNames')
verifyError(testCase,@() recording.processSignal(processor,{'ecg',2}), ...
    'mpk:Recording:InvalidSignalNames')
verifyError(testCase,@() recording.processSignal(processor,42), ...
    'mpk:Recording:InvalidSignalNames')
verifyError(testCase,@() recording.processSignal( ...
    processor,{'ecg','ecg'}), ...
    'mpk:Recording:DuplicateSignalNames')
verifyError(testCase,@() recording.processSignal(processor,'missing'), ...
    'mpk:Recording:MissingSignal')
verifyError(testCase,@() recording.processSignal(42,'ecg'), ...
    'mpk:Recording:InvalidProcessor')
verifyError(testCase,@() recording.processSignal( ...
    struct('process',@sin),'ecg'), ...
    'mpk:Recording:InvalidProcessor')
verifyError(testCase,@() recording.processSignal( ...
    processor,'ecg','ProcessArguments','not-a-cell'), ...
    'mpk:Recording:InvalidProcessArguments')
end

function testProcessSignalForwardsProcessArgumentsUnchanged(testCase)
recording = localProcessingRecording();
processor = FakeRecordingProcessor();
processArguments = {'BaselineRedValue',3.5, ...
    'BaselineInfraredValue',7.25,17};
originalArguments = processArguments;

recording.processSignal( ...
    processor,{'ecg','temperature'}, ...
    'ProcessArguments',processArguments);

received = processor.Calls{1};
verifyEqual(testCase,numel(received),2 + numel(processArguments))
verifyEqual(testCase,received{1}.Label,'ecg')
verifyEqual(testCase,received{2}.Label,'temperature')
verifyEqual(testCase,received(3:end),processArguments)
verifyEqual(testCase,processArguments,originalArguments)
end

function testProcessSignalPropagatesProcessorErrors(testCase)
recording = localProcessingRecording();
processor = FakeRecordingProcessor();
processor.ErrorOnCall = 1;

verifyError(testCase,@() recording.processSignal(processor,'ecg'), ...
    'test:FakeRecordingProcessor:ProcessFailed')
end

function testProcessPhasesRequiresTemporalPhases(testCase)
recording = localProcessingRecording();
processor = FakeRecordingProcessor();

verifyError(testCase,@() recording.processPhases(processor,'ecg'), ...
    'mpk:Recording:NoPhases')
end

function testProcessPhasesProcessesAllAndNumericEmptyInOrder(testCase)
recording = localProcessingRecording();
recording.setPhases([0 2; 2 4],{'first phase','second-phase'});
processor = FakeRecordingProcessor();
processor.ReturnCallNumber = true;

results = recording.processPhases(processor,'ecg');

verifyClass(testCase,results,'containers.Map')
verifyEqual(testCase,results.KeyType,'char')
verifyEqual(testCase,results.ValueType,'any')
verifyTrue(testCase,isKey(results,'first phase'))
verifyTrue(testCase,isKey(results,'second-phase'))
firstResult = results('first phase');
secondResult = results('second-phase');
verifyEqual(testCase,firstResult.CallNumber,1)
verifyEqual(testCase,secondResult.CallNumber,2)
verifyEqual(testCase,firstResult.Features,struct('unchanged',1))
verifyEqual(testCase,numel(processor.Calls),2)
verifyEqual(testCase,processor.Calls{1}{1}.Metadata.PhaseLabel, ...
    'first phase')
verifyEqual(testCase,processor.Calls{2}{1}.Metadata.PhaseLabel, ...
    'second-phase')
verifyEqual(testCase,processor.Calls{1}{1}.Time(1),0)
verifyEqual(testCase,processor.Calls{2}{1}.Time(1),0)

processor = FakeRecordingProcessor();
recording.processPhases(processor,'ecg','Phases',[]);
verifyEqual(testCase,numel(processor.Calls),2)
end

function testProcessPhasesAcceptsSelectorsAndPreservesRequestedOrder(testCase)
recording = localProcessingRecording();
recording.setPhases([0 2; 2 4],{'first','second'});

processor = FakeRecordingProcessor();
recording.processPhases(processor,'ecg','Phases','second');
verifyEqual(testCase,numel(processor.Calls),1)
verifyEqual(testCase,processor.Calls{1}{1}.Metadata.PhaseLabel,'second')

processor = FakeRecordingProcessor();
recording.processPhases(processor,'ecg','Phases',"first");
verifyEqual(testCase,processor.Calls{1}{1}.Metadata.PhaseLabel,'first')

processor = FakeRecordingProcessor();
recording.processPhases( ...
    processor,'ecg','Phases',{'second','first'});
verifyEqual(testCase,processor.Calls{1}{1}.Metadata.PhaseLabel,'second')
verifyEqual(testCase,processor.Calls{2}{1}.Metadata.PhaseLabel,'first')

processor = FakeRecordingProcessor();
recording.processPhases( ...
    processor,'ecg','Phases',["second","first"]);
verifyEqual(testCase,processor.Calls{1}{1}.Metadata.PhaseLabel,'second')
verifyEqual(testCase,processor.Calls{2}{1}.Metadata.PhaseLabel,'first')
end

function testProcessPhasesRejectsInvalidSelectionsAndResetTime(testCase)
recording = localProcessingRecording();
recording.setPhases([0 2; 2 4],{'first','second'});
processor = FakeRecordingProcessor();

verifyError(testCase,@() recording.processPhases( ...
    processor,'ecg','Phases',cell(1,0)), ...
    'mpk:Recording:InvalidPhaseSelection')
verifyError(testCase,@() recording.processPhases( ...
    processor,'ecg','Phases',strings(1,0)), ...
    'mpk:Recording:InvalidPhaseSelection')
verifyError(testCase,@() recording.processPhases( ...
    processor,'ecg','Phases',''), ...
    'mpk:Recording:InvalidPhaseSelection')
verifyError(testCase,@() recording.processPhases( ...
    processor,'ecg','Phases',{'missing'}), ...
    'mpk:Recording:UnknownPhase')
verifyError(testCase,@() recording.processPhases( ...
    processor,'ecg','ResetTime',1), ...
    'mpk:Recording:InvalidResetTime')
verifyError(testCase,@() recording.processPhases( ...
    processor,'ecg','ProcessArguments',42), ...
    'mpk:Recording:InvalidProcessArguments')
verifyError(testCase,@() recording.processPhases(42,'ecg'), ...
    'mpk:Recording:InvalidProcessor')
end

function testProcessPhasesCanPreserveTimeAndForwardArguments(testCase)
recording = localProcessingRecording();
recording.setPhases([1 2; 2 3],{'one','two'});
processor = FakeRecordingProcessor();
processArguments = {'OptionName',17,@sin};

recording.processPhases( ...
    processor,'ecg', ...
    'ResetTime',false, ...
    'ProcessArguments',processArguments);

verifyEqual(testCase,numel(processor.Calls),2)
verifyEqual(testCase,processor.Calls{1}{1}.Time(1),1)
verifyEqual(testCase,processor.Calls{2}{1}.Time(1),2)
verifyEqual(testCase,processor.Calls{1}(2:end),processArguments)
verifyEqual(testCase,processor.Calls{2}(2:end),processArguments)
end

function testProcessPhasesExtractsMultipleSignalsIndependently(testCase)
recording = mpk.Recording('multimodal');
ecgFs = 1000;
temperatureFs = 10;
ecgTime = (0:(30*ecgFs-1))' ./ ecgFs;
temperatureTime = (0:(30*temperatureFs-1))' ./ temperatureFs;
ecg = mpk.Signal(zeros(size(ecgTime)),ecgTime,ecgFs,'ecg','mV');
temperature = mpk.Signal(zeros(size(temperatureTime)), ...
    temperatureTime,temperatureFs,'temperature','degC');
recording.addSignal('ecg',ecg);
recording.addSignal('temperature',temperature);
recording.setPhases([10;20],{'shared phase'});
originalIntervals = recording.PhaseIntervalsSeconds;
originalLabels = recording.PhaseLabels;
processor = FakeRecordingProcessor();

results = recording.processPhases( ...
    processor,{'temperature','ecg'}, ...
    'ResetTime',false);

verifyTrue(testCase,isKey(results,'shared phase'))
verifyEqual(testCase,numel(processor.Calls),1)
received = processor.Calls{1};
verifyEqual(testCase,received{1}.Label,'temperature')
verifyEqual(testCase,received{2}.Label,'ecg')
verifyEqual(testCase,received{1}.count(),100)
verifyEqual(testCase,received{2}.count(),10000)
verifyEqual(testCase,[received{1}.Time(1),received{1}.Time(end)], ...
    [10,19.9],'AbsTol',1e-12)
verifyEqual(testCase,[received{2}.Time(1),received{2}.Time(end)], ...
    [10,19.999],'AbsTol',1e-12)
verifyEqual(testCase,temperature.Time,temperatureTime)
verifyEqual(testCase,ecg.Time,ecgTime)
verifyEqual(testCase,recording.PhaseIntervalsSeconds,originalIntervals)
verifyEqual(testCase,recording.PhaseLabels,originalLabels)
end

function testProcessPhasesPropagatesErrorsAndStops(testCase)
recording = localProcessingRecording();
recording.setPhases([0 1 2; 1 2 3],{'one','two','three'});
processor = FakeRecordingProcessor();
processor.ErrorOnCall = 2;

verifyError(testCase,@() recording.processPhases(processor,'ecg'), ...
    'test:FakeRecordingProcessor:ProcessFailed')
verifyEqual(testCase,numel(processor.Calls),2)
end

function signal = localSignal(label)
signal = mpk.Signal([1; 2; 3],[0; 1; 2],1,label,'a.u.');
end

function recording = localProcessingRecording()
recording = mpk.Recording('processing');
time = (0:0.5:4.5)';
recording.addSignal('ecg',mpk.Signal( ...
    (1:numel(time))',time,2,'ecg','mV'));
recording.addSignal('temperature',mpk.Signal( ...
    (101:(100 + numel(time)))',time,2,'temperature','degC'));
end

function verifySameHandleOrValue(testCase,actual,expected)
verifyEqual(testCase,actual.Samples,expected.Samples)
verifyEqual(testCase,actual.Time,expected.Time)
verifyEqual(testCase,actual.SamplingFrequency,expected.SamplingFrequency)
verifyEqual(testCase,actual.Label,expected.Label)
verifyEqual(testCase,actual.Units,expected.Units)
verifyEqual(testCase,actual.Metadata,expected.Metadata)
end
