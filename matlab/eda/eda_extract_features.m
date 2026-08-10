function features = eda_extract_features(decomposition,minimumScrAmplitude)
%EDA_EXTRACT_FEATURES Extract the canonical eight EDA features.

if nargin < 2
    minimumScrAmplitude = 0.01;
end

scl = decomposition.Scl;
scr = decomposition.Scr;
time = decomposition.Time;
if numel(scl) < 3 || numel(scr) < 3 || numel(time) < 3
    error('mpk:EDAProcessor:SignalTooShortForFeatures', ...
        'EDA signal is too short for feature extraction.')
end
if any(~isfinite(scl)) || any(~isfinite(scr)) || ...
        any(~isfinite(time))
    error('mpk:EDAProcessor:NonfiniteFeatureInput', ...
        'EDA feature inputs must contain finite values.')
end
duration = time(end) - time(1);
if ~isfinite(duration) || duration <= 0
    error('mpk:EDAProcessor:NonpositiveDuration', ...
        'EDA signal duration must be positive.')
end

try
    [amplitudeScr,riseTimeScr,numberScr,~,~] = ...
        scr_processing( ...
            scr,time,decomposition.SamplingFrequency,minimumScrAmplitude);
catch exception
    rethrow(exception)
end
amplitudeScr = amplitudeScr(:);
riseTimeScr = riseTimeScr(:);
if any(~isfinite(amplitudeScr)) || any(~isfinite(riseTimeScr))
    error('mpk:EDAProcessor:NonfiniteSCRFeatures', ...
        'SCR event features must contain finite values.')
end

slope = polyfit(1:numel(scl),scl.',1);
if numberScr == 0
    scrFeatureValues = zeros(4,1);
else
    scrFeatureValues = [
        mean(amplitudeScr)
        std(amplitudeScr)
        mean(riseTimeScr)
        std(riseTimeScr)
        ];
end
featureValues = [
    mean(scl)
    std(scl)
    slope(1)
    scrFeatureValues
    numberScr / duration * 60
    ];
if any(~isfinite(featureValues))
    error('mpk:EDAProcessor:NonfiniteFeatures', ...
        'EDA feature extraction produced a nonfinite value.')
end
featureNames = {
    'eda_scl_mean'
    'eda_scl_std'
    'eda_scl_slope'
    'eda_scr_amplitude_mean'
    'eda_scr_amplitude_std'
    'eda_scr_rise_time_mean'
    'eda_scr_rise_time_std'
    'eda_scr_per_min'
    };
features = [featureNames,num2cell(featureValues)];
end
