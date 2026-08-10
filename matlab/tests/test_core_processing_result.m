function tests = test_core_processing_result
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

function testConstructsProcessingResult(testCase)
features = {'mean_hr',72; 'sdnn',42};
options = struct('ShowPlots',false);
metadata = struct('LegacyFunction','ecg_signal_processing');

result = mpk.ProcessingResult('ecg','ECGProcessor',features, ...
    'ecg',options,metadata);

verifyClass(testCase,result,'mpk.ProcessingResult')
verifyEqual(testCase,result.Modality,'ecg')
verifyEqual(testCase,result.ProcessorName,'ECGProcessor')
verifyEqual(testCase,result.Features,features)
verifyEqual(testCase,result.SourceSignalLabel,'ecg')
verifyEqual(testCase,result.Options,options)
verifyEqual(testCase,result.Metadata,metadata)
verifyEqual(testCase,result.Intermediates,struct())
verifyTrue(testCase,ischar(result.CreatedOn))
verifyNotEmpty(testCase,result.CreatedOn)
end

function testUsesDefaultOptionalFields(testCase)
features = {'feature',1};

result = mpk.ProcessingResult('eda','EDAProcessor',features);

verifyEqual(testCase,result.SourceSignalLabel,'')
verifyEqual(testCase,result.Options,struct())
verifyEqual(testCase,result.Metadata,struct())
verifyEqual(testCase,result.Intermediates,struct())
end

function testSixArgumentConstructionRemainsCompatible(testCase)
features = {'feature',1};
options = struct('ShowPlots',false);
metadata = struct('source','legacy_processor');

result = mpk.ProcessingResult( ...
    'modality','Processor',features,'signal',options,metadata);

verifyEqual(testCase,result.Features,features)
verifyEqual(testCase,result.SourceSignalLabel,'signal')
verifyEqual(testCase,result.Options,options)
verifyEqual(testCase,result.Metadata,metadata)
verifyEqual(testCase,result.Intermediates,struct())
end

function testExplicitEmptyIntermediatesUseScalarStruct(testCase)
result = mpk.ProcessingResult( ...
    'ecg','ECGProcessor',{'feature',1},'',struct(),struct(),[]);

verifyTrue(testCase,isstruct(result.Intermediates))
verifyTrue(testCase,isscalar(result.Intermediates))
verifyEmpty(testCase,fieldnames(result.Intermediates))
end

function testPreservesStructuredIntermediateValues(testCase)
intermediates = struct( ...
    'filtered_signal',[1;2;3], ...
    'r_peak_indices',[2;5;9], ...
    'ibi',[0.8;0.9], ...
    'stft_power',ones(3,4), ...
    'diagnostics',struct('accepted_peaks',3));
intermediates.notes = {'accepted',[1 2 3]};

result = mpk.ProcessingResult( ...
    'ecg','ECGProcessor',{'feature',1}, ...
    'ecg',struct(),struct(),intermediates);

verifyEqual(testCase,fieldnames(result.Intermediates), ...
    fieldnames(intermediates))
verifyEqual(testCase,result.Intermediates.filtered_signal,[1;2;3])
verifyEqual(testCase,result.Intermediates.r_peak_indices,[2;5;9])
verifyEqual(testCase,result.Intermediates.ibi,[0.8;0.9])
verifyEqual(testCase,result.Intermediates.stft_power,ones(3,4))
verifyEqual(testCase,result.Intermediates.diagnostics, ...
    struct('accepted_peaks',3))
verifyEqual(testCase,result.Intermediates.notes, ...
    {'accepted',[1 2 3]})
verifyEqual(testCase,result.Intermediates,intermediates)
end

function testRejectsInvalidIntermediateContainers(testCase)
features = {'feature',1};
invalidValues = {
    repmat(struct('value',1),1,2)
    [1 2 3]
    true
    'filtered_signal'
    "filtered_signal"
    {[1;2;3]}
    table([1;2])
    containers.Map({'filtered_signal'},{[1;2;3]})
    };

for index = 1:numel(invalidValues)
    verifyError(testCase,@() mpk.ProcessingResult( ...
        'ecg','ECGProcessor',features,'ecg',struct(),struct(), ...
        invalidValues{index}), ...
        'mpk:ProcessingResult:InvalidIntermediates')
end
end

function testToCellReturnsFeatureValueUnchanged(testCase)
features = {'feature',1};
intermediates = struct('filtered_signal',[1;2;3]);
result = mpk.ProcessingResult( ...
    'temp','TemperatureProcessor',features,'',struct(),struct(), ...
    intermediates);

verifyEqual(testCase,result.toCell(),features)
end

function testFeatureNamesFromCellFeatures(testCase)
features = {'mean_hr',72; 'sdnn',42};
intermediates = struct( ...
    'filtered_signal',[1;2;3], ...
    'r_peak_indices',[2;5;9]);
result = mpk.ProcessingResult( ...
    'ecg','ECGProcessor',features,'',struct(),struct(), ...
    intermediates);

verifyEqual(testCase,result.featureNames(),{'mean_hr'; 'sdnn'})
verifyFalse(testCase,any(ismember( ...
    fieldnames(intermediates),result.featureNames())))
end

function testFeatureNamesFromStringFeatures(testCase)
features = ["mean_temp" "36.5"; "std_temp" "0.2"];
result = mpk.ProcessingResult('temp','TemperatureProcessor',features);

verifyEqual(testCase,result.featureNames(),{'mean_temp'; 'std_temp'})
end

function testFeatureNamesFromStructFeatureGroups(testCase)
features = struct();
features.Oxy = {'mean_oxy',1; 'std_oxy',2};
features.Deoxy = ["mean_deoxy" "3"; "std_deoxy" "4"];
features.Empty = [1 2 3];
result = mpk.ProcessingResult('fnirs','FNIRSProcessor',features);

names = result.featureNames();

verifyEqual(testCase,names.Oxy,{'mean_oxy'; 'std_oxy'})
verifyEqual(testCase,names.Deoxy,{'mean_deoxy'; 'std_deoxy'})
verifyEqual(testCase,names.Empty,{})
end

function testFeatureNamesReturnsEmptyForUnsupportedFeatureType(testCase)
result = mpk.ProcessingResult('custom','Processor',[1 2 3]);

verifyEqual(testCase,result.featureNames(),{})
end
