function [amplitudesChecked,riseTimesChecked,count,indexesCheckedMax,indexesCheckedMin] = ...
        scr_processing(scr,time,fs,minimumScrAmplitude)
%SCR_PROCESSING Detect SCRs using filtered-phasic trough-to-peak amplitude.
%
% Amplitude is the filtered phasic value at a local maximum minus the value
% at its immediately preceding local minimum, in microSiemens. A candidate
% is accepted when this amplitude is at least MINIMUMSCRAMPLITUDE (default
% 0.01 microSiemens). Rise time is reported but is not used for acceptance.

if nargin < 4
    minimumScrAmplitude = 0.01;
end
if ~isnumeric(minimumScrAmplitude) || ~isreal(minimumScrAmplitude) || ...
        ~isscalar(minimumScrAmplitude) || ~isfinite(minimumScrAmplitude) || ...
        minimumScrAmplitude < 0
    error('mpk:EDAProcessor:InvalidMinimumScrAmplitude', ...
        ['minimumScrAmplitude must be a finite numeric scalar greater ' ...
         'than or equal to zero.'])
end

scrFiltered = filtering(scr,fs,0,2,0,2,0);
[maxima,minima] = find_max_min(scrFiltered,0,0);
[amplitudes,riseTimes,pairedMaxima,pairedMinima] = ...
    scr_pair_events(scrFiltered,maxima,minima,fs);

if isempty(amplitudes)
    amplitudesChecked = zeros(1,0);
    riseTimesChecked = zeros(1,0);
    indexesCheckedMax = zeros(1,0);
    indexesCheckedMin = zeros(1,0);
    count = 0;
    return
end

keep = scr_accept_events(amplitudes,minimumScrAmplitude);
amplitudesChecked = amplitudes(keep);
riseTimesChecked = riseTimes(keep);
indexesCheckedMax = pairedMaxima(keep);
indexesCheckedMin = pairedMinima(keep);
count = numel(amplitudesChecked);

if any(~isfinite(amplitudesChecked)) || any(~isfinite(riseTimesChecked))
    error('mpk:EDAProcessor:NonfiniteSCRFeatures', ...
        'SCR event features must contain finite values.')
end
end
