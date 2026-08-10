function tests = test_fnirs_processor_validation
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

function testAgeIsRequiredAndConfigurationIsValidated(testCase)
verifyError(testCase,@() mpk.FNIRSProcessor(), ...
    'mpk:FNIRSProcessor:InvalidAge')

invalidAges = {[],true,[20 30],0,-1,NaN,Inf,1+1i,'25'};
for index = 1:numel(invalidAges)
    verifyError(testCase,@() mpk.FNIRSProcessor( ...
        'AgeYears',invalidAges{index}), ...
        'mpk:FNIRSProcessor:InvalidAge')
end
verifyError(testCase,@() mpk.FNIRSProcessor( ...
    'AgeYears',25,'ReturnIntermediates',1), ...
    'mpk:FNIRSProcessor:InvalidReturnIntermediates')

processor = mpk.FNIRSProcessor( ...
    'AgeYears',25,'ReturnIntermediates',true);
verifyEqual(testCase,processor.AgeYears,25)
verifyTrue(testCase,processor.ReturnIntermediates)
end

function testRejectsInvalidSignalInputs(testCase)
processor = localProcessor();
[red,infrared] = localSignals(30,100);

verifyError(testCase,@() processor.process(1,infrared), ...
    'mpk:FNIRSProcessor:InvalidRedSignal')
verifyError(testCase,@() processor.process(red,1), ...
    'mpk:FNIRSProcessor:InvalidInfraredSignal')
verifyError(testCase,@() processor.process({red,infrared},infrared), ...
    'mpk:FNIRSProcessor:InvalidRedSignal')
verifyError(testCase,@() processor.process( ...
    mpk.Recording('invalid'),infrared), ...
    'mpk:FNIRSProcessor:InvalidRedSignal')
verifyError(testCase,@() processor.process(red,struct()), ...
    'mpk:FNIRSProcessor:InvalidInfraredSignal')
end

function testRejectsMismatchedSignalPairs(testCase)
processor = localProcessor();
[red,infrared] = localSignals(30,100);
shortInfrared = mpk.Signal( ...
    infrared.Samples(1:end-1),infrared.Time(1:end-1), ...
    100,'infrared','uA');
verifyError(testCase,@() processor.process(red,shortInfrared), ...
    'mpk:FNIRSProcessor:LengthMismatch')

differentFrequency = mpk.Signal( ...
    infrared.Samples,infrared.Time,200,'infrared','uA');
verifyError(testCase,@() processor.process(red,differentFrequency), ...
    'mpk:FNIRSProcessor:SamplingFrequencyMismatch')

differentTime = mpk.Signal( ...
    infrared.Samples,infrared.Time + 0.01,100,'infrared','uA');
verifyError(testCase,@() processor.process(red,differentTime), ...
    'mpk:FNIRSProcessor:TimeMismatch')
end

function testRejectsInvalidCurrentsAndSamplingFrequency(testCase)
processor = localProcessor();
[red,infrared] = localSignals(30,100);
invalidSamples = red.Samples;
invalidSamples(5) = 0;
invalidRed = mpk.Signal( ...
    invalidSamples,red.Time,red.SamplingFrequency,'red','uA');
verifyError(testCase,@() processor.process(invalidRed,infrared), ...
    'mpk:FNIRSProcessor:InvalidCurrentSamples')

nanSamples = infrared.Samples;
nanSamples(7) = NaN;
invalidInfrared = mpk.Signal( ...
    nanSamples,infrared.Time,infrared.SamplingFrequency, ...
    'infrared','uA');
verifyError(testCase,@() processor.process(red,invalidInfrared), ...
    'mpk:FNIRSProcessor:InvalidCurrentSamples')

[red125,infrared125] = localSignals(30,125);
verifyError(testCase,@() processor.process(red125,infrared125), ...
    'mpk:FNIRSProcessor:InvalidSamplingFrequency')
end

function testBaselinePairValidation(testCase)
processor = localProcessor();
[red,infrared] = localSignals(30,100);
verifyError(testCase,@() processor.process( ...
    red,infrared,'BaselineRedValue',10), ...
    'mpk:FNIRSProcessor:IncompleteBaselinePair')
verifyError(testCase,@() processor.process( ...
    red,infrared,'BaselineInfraredValue',11), ...
    'mpk:FNIRSProcessor:IncompleteBaselinePair')

invalidValues = {true,[1 2],1+1i,NaN,Inf,0,-1};
for index = 1:numel(invalidValues)
    invalid = invalidValues{index};
    verifyError(testCase,@() processor.process( ...
        red,infrared,'BaselineRedValue',invalid, ...
        'BaselineInfraredValue',11), ...
        'mpk:FNIRSProcessor:InvalidBaselineRedValue')
    verifyError(testCase,@() processor.process( ...
        red,infrared,'BaselineRedValue',10, ...
        'BaselineInfraredValue',invalid), ...
        'mpk:FNIRSProcessor:InvalidBaselineInfraredValue')
end
end

function processor = localProcessor()
processor = mpk.FNIRSProcessor('AgeYears',25);
end

function [redSignal,infraredSignal] = localSignals(durationSeconds,fs)
elapsed = (0:1/fs:durationSeconds - 1/fs)';
time = 4.2 + elapsed;
red = 10 + 0.04*sin(2*pi*0.03*elapsed);
infrared = 11 + 0.04*cos(2*pi*0.025*elapsed);
redSignal = mpk.Signal(red,time,fs,'red_current','uA');
infraredSignal = mpk.Signal( ...
    infrared,time,fs,'infrared_current','uA');
end
