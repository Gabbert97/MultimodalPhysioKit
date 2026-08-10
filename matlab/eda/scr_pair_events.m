function [amplitudes,riseTimes,pairedMaxima,pairedMinima] = ...
        scr_pair_events(filteredPhasic,maxima,minima,fs)
%SCR_PAIR_EVENTS Pair maxima with preceding minima and calculate SCR values.
%
% AMPLITUDES are filtered-phasic trough-to-peak changes in microSiemens.

maxima = maxima(:).';
minima = minima(:).';
pairedMaxima = zeros(1,0);
pairedMinima = zeros(1,0);
for index = 1:numel(maxima)
    preceding = minima(minima < maxima(index));
    if isempty(preceding)
        continue
    end
    pairedMaxima(end+1) = maxima(index); %#ok<AGROW>
    pairedMinima(end+1) = preceding(end); %#ok<AGROW>
end
amplitudes = filteredPhasic(pairedMaxima) - ...
    filteredPhasic(pairedMinima);
riseTimes = (pairedMaxima - pairedMinima) / fs;
end
