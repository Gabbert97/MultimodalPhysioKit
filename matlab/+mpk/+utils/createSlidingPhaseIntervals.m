function [intervalsSeconds,labels] = createSlidingPhaseIntervals( ...
        startTimeSeconds,endTimeSeconds,windowDurationSeconds, ...
        overlapSeconds,varargin)
%CREATESLIDINGPHASEINTERVALS Create complete half-open temporal windows.
%
% [INTERVALS,LABELS] = mpk.utils.createSlidingPhaseIntervals(START,STOP,
% WINDOW,OVERLAP) returns a 2-by-N floating-point matrix in seconds. Each
% column is [start,end), and only complete windows are returned.

if ~isFiniteRealScalar(startTimeSeconds) || startTimeSeconds < 0
    error('mpk:utils:createSlidingPhaseIntervals:InvalidStart', ...
        'startTimeSeconds must be a nonnegative finite real scalar.')
end
if ~isFiniteRealScalar(endTimeSeconds) || ...
        endTimeSeconds <= startTimeSeconds
    error('mpk:utils:createSlidingPhaseIntervals:InvalidEnd', ...
        'endTimeSeconds must be finite and greater than startTimeSeconds.')
end
if ~isFiniteRealScalar(windowDurationSeconds) || ...
        windowDurationSeconds <= 0
    error('mpk:utils:createSlidingPhaseIntervals:InvalidDuration', ...
        'windowDurationSeconds must be a positive finite real scalar.')
end
if ~isFiniteRealScalar(overlapSeconds) || overlapSeconds < 0 || ...
        overlapSeconds >= windowDurationSeconds
    error('MultimodalPhysioKit:InvalidPhaseOverlap', ...
        ['overlapSeconds must be nonnegative and shorter than ' ...
         'windowDurationSeconds.'])
end

parser = inputParser;
parser.FunctionName = 'mpk.utils.createSlidingPhaseIntervals';
addParameter(parser,'LabelPrefix','window');
addParameter(parser,'WarnOnDroppedTrailingWindow',true);
addParameter(parser,'MinimumWindowDurationSeconds',0);
parse(parser,varargin{:});

labelPrefix = parser.Results.LabelPrefix;
if ~(ischar(labelPrefix) && isrow(labelPrefix)) && ...
        ~(isstring(labelPrefix) && isscalar(labelPrefix))
    error('mpk:utils:createSlidingPhaseIntervals:InvalidLabelPrefix', ...
        'LabelPrefix must be a character vector or string scalar.')
end
labelPrefix = char(labelPrefix);
if isempty(labelPrefix)
    error('mpk:utils:createSlidingPhaseIntervals:InvalidLabelPrefix', ...
        'LabelPrefix must not be empty.')
end
warnOnTrailing = parser.Results.WarnOnDroppedTrailingWindow;
if ~islogical(warnOnTrailing) || ~isscalar(warnOnTrailing)
    error('mpk:utils:createSlidingPhaseIntervals:InvalidWarningFlag', ...
        'WarnOnDroppedTrailingWindow must be a logical scalar.')
end
minimumDuration = parser.Results.MinimumWindowDurationSeconds;
if ~isFiniteRealScalar(minimumDuration) || minimumDuration < 0
    error('mpk:utils:createSlidingPhaseIntervals:InvalidMinimumDuration', ...
        ['MinimumWindowDurationSeconds must be a nonnegative finite ' ...
         'real scalar.'])
end
if windowDurationSeconds < minimumDuration
    error('mpk:utils:createSlidingPhaseIntervals:InvalidDuration', ...
        ['windowDurationSeconds must be at least ' ...
         'MinimumWindowDurationSeconds.'])
end

stepSeconds = windowDurationSeconds - overlapSeconds;
availableDuration = endTimeSeconds - startTimeSeconds;
windowCount = floor((availableDuration - windowDurationSeconds) / ...
    stepSeconds + 16*eps(max(1,availableDuration))) + 1;
if windowCount < 1
    error('MultimodalPhysioKit:SampleTooShort', ...
        'No complete temporal window fits inside the requested range.')
end

windowIndexes = 0:(windowCount - 1);
starts = startTimeSeconds + windowIndexes .* stepSeconds;
ends = starts + windowDurationSeconds;
tolerance = 16*eps(max(1,endTimeSeconds));
complete = ends <= endTimeSeconds + tolerance;
starts = starts(complete);
ends = min(ends(complete),endTimeSeconds);
intervalsSeconds = double([starts;ends]);
labels = arrayfun( ...
    @(index) sprintf('%s_%d',labelPrefix,index), ...
    1:numel(starts),'UniformOutput',false);

nextStart = startTimeSeconds + numel(starts) * stepSeconds;
remaining = endTimeSeconds - nextStart;
if warnOnTrailing && nextStart < endTimeSeconds - tolerance && ...
        remaining < windowDurationSeconds - tolerance
    warning('MultimodalPhysioKit:DroppedIncompletePhase', ...
        ['Trailing phase is dropped because only %.3f s are available; ' ...
         '%.3f s are required.'],remaining,windowDurationSeconds)
end
end

function tf = isFiniteRealScalar(value)
tf = isnumeric(value) && ~islogical(value) && isreal(value) && ...
    isscalar(value) && isfinite(value);
end
