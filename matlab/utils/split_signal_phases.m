function signals = split_signal_phases(i_phases)
%SPLIT_SIGNAL_PHASES Expand phase boundaries into per-phase index vectors.
%
% I_PHASES is a 2-by-N matrix of one-based inclusive phase limits. The output
% is an N-by-1 cell array where each cell contains the sample indexes for one
% phase.

[~,column] = size(i_phases);
signals = {};

for i = 1 : column
    signals{i,1} = (i_phases(1,i) : i_phases(2,i))';
end
end
