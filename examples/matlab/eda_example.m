% Process two EDA segments from the bundled synthetic fixture.
% EDA requires external cvxEDA and its required optimization functionality.

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

edaResults = recording.processPhases( ...
    mpk.EDAProcessor(), ...
    'eda', ...
    'ProcessArguments',{'ShowPlots',showPlots});

edaResults('segment_1').Features
