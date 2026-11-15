% Process two independent PPG windows with phase-local STFT analysis.

showPlots = true;

exampleDirectory = fileparts(mfilename('fullpath'));
repositoryRoot = fileparts(fileparts(exampleDirectory));
addpath(genpath(fullfile(repositoryRoot,'matlab')))

dataFile = fullfile(repositoryRoot,'examples','sample_data', ...
    'synthetic_multimodal_180s.hdf5');
recording = mpk.io.readBioSignalsPluxHdf5(dataFile);

phaseIntervals = [0 90;90 180];
phaseLabels = {'segment_1','segment_2'};
recording.setPhases(phaseIntervals,phaseLabels);

if showPlots
    mpk.utils.plotSignalWithPhases( ...
        recording.getSignal('ppg'), ...
        recording.PhaseIntervalsSeconds, ...
        'PhaseLabels',recording.PhaseLabels, ...
        'Title','Raw PPG with analysis phases');
end

processor = mpk.PPGProcessor( ...
    'LabelFrequencyAnalysis',1, ...
    'ReturnIntermediates',true);
processor.ShowPlots = showPlots;
ppgResults = recording.processPhases(processor,'ppg');

ppgResults('segment_1').Features
