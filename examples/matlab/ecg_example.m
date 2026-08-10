% Process two independent ECG windows with phase-local STFT analysis.

showPlots = true;

exampleDirectory = fileparts(mfilename('fullpath'));
repositoryRoot = fileparts(fileparts(exampleDirectory));
addpath(genpath(fullfile(repositoryRoot,'matlab')))

dataFile = fullfile(repositoryRoot,'examples','sample_data', ...
    'synthetic_multimodal_180s.hdf5');
recording = mpk.io.readBioSignalsPluxHdf5(dataFile);

phaseIntervals = [0 60;60 120];
phaseLabels = {'segment_1','segment_2'};
recording.setPhases(phaseIntervals,phaseLabels);

if showPlots
    mpk.utils.plotSignalWithPhases( ...
        recording.getSignal('ecg'), ...
        recording.PhaseIntervalsSeconds, ...
        'PhaseLabels',recording.PhaseLabels, ...
        'Title','Raw ECG with analysis phases');
end

processor = mpk.ECGProcessor( ...
    'LabelFrequencyAnalysis',1, ...
    'ReturnIntermediates',true);
ecgResults = recording.processPhases(processor,'ecg');

ecgResults('segment_1').Features
