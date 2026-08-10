classdef Signal
    %SIGNAL Sampled physiological signal and acquisition metadata.
    %
    % SIGNAL = mpk.Signal(SAMPLES,TIME,FS) stores numeric samples, matching
    % timestamps in seconds, and a positive sampling frequency in hertz.
    % Optional LABEL, UNITS, and METADATA arguments describe the signal
    % without retaining a source file path. Sample and time vectors are
    % normalized to columns.

    properties
        Samples
        Time
        SamplingFrequency
        Label
        Units
        Metadata
    end

    methods
        function obj = Signal(samples,time,fs,label,units,metadata)
            if nargin < 4 || isempty(label)
                label = '';
            end
            if nargin < 5 || isempty(units)
                units = '';
            end
            if nargin < 6 || isempty(metadata)
                metadata = struct();
            end

            if ~isnumeric(samples) || ~isvector(samples)
                error('mpk:Signal:InvalidSamples', ...
                    'Samples must be a numeric vector.')
            end
            if ~isnumeric(time) || ~isvector(time)
                error('mpk:Signal:InvalidTime', ...
                    'Time must be a numeric vector.')
            end
            if numel(samples) ~= numel(time)
                error('mpk:Signal:LengthMismatch', ...
                    'Samples and Time must have the same length.')
            end
            if ~isscalar(fs) || ~isnumeric(fs) || ~isfinite(fs) || fs <= 0
                error('mpk:Signal:InvalidSamplingFrequency', ...
                    'SamplingFrequency must be a positive numeric scalar.')
            end
            if ~isstruct(metadata)
                error('mpk:Signal:InvalidMetadata', ...
                    'Metadata must be a struct.')
            end

            obj.Samples = samples(:);
            obj.Time = time(:);
            obj.SamplingFrequency = fs;
            obj.Label = char(label);
            obj.Units = char(units);
            obj.Metadata = metadata;
        end

        function n = count(obj)
            %COUNT Return the number of signal samples.
            n = numel(obj.Samples);
        end

        function seconds = durationSeconds(obj)
            %DURATIONSECONDS Return elapsed time between first and last sample.
            if isempty(obj.Time)
                seconds = 0;
            else
                seconds = obj.Time(end) - obj.Time(1);
            end
        end

        function out = removeTimeRanges(obj,intervalsSeconds,varargin)
            %REMOVETIMERANGES Return a copy with half-open time ranges removed.
            %
            % OUT = removeTimeRanges(OBJ,INTERVALSSECONDS) removes samples
            % whose timestamps fall in any [start,end) interval. Use Inf as
            % the end value for open-ended removal.

            parser = inputParser;
            parser.FunctionName = 'mpk.Signal.removeTimeRanges';
            addParameter(parser,'ResetTime',true, ...
                @(value) islogical(value) && isscalar(value));
            parse(parser,varargin{:});
            resetTime = logical(parser.Results.ResetTime);

            mpk.Signal.validateTimeRanges(intervalsSeconds);

            keep = true(size(obj.Time));
            for index = 1:size(intervalsSeconds,1)
                startSeconds = intervalsSeconds(index,1);
                endSeconds = intervalsSeconds(index,2);
                if isinf(endSeconds)
                    remove = obj.Time >= startSeconds;
                else
                    remove = obj.Time >= startSeconds & obj.Time < endSeconds;
                end
                keep = keep & ~remove;
            end

            if ~any(keep)
                error('mpk:Signal:EmptyAfterRemoval', ...
                    'Removing time ranges would remove all samples.')
            end

            samples = obj.Samples(keep);
            if resetTime
                time = (0:(numel(samples)-1)).' ./ obj.SamplingFrequency;
            else
                time = obj.Time(keep);
            end

            metadata = obj.Metadata;
            metadata.TimeRangesRemovedSeconds = intervalsSeconds;
            metadata.ResetTimeAfterRemoval = resetTime;
            out = mpk.Signal(samples,time,obj.SamplingFrequency, ...
                obj.Label,obj.Units,metadata);
        end
    end

    methods (Static, Access = private)
        function validateTimeRanges(intervalsSeconds)
            if ~isnumeric(intervalsSeconds) || size(intervalsSeconds,2) ~= 2 || ...
                    isempty(intervalsSeconds)
                error('mpk:Signal:InvalidTimeRanges', ...
                    'intervalsSeconds must be a nonempty Nx2 numeric matrix.')
            end
            starts = intervalsSeconds(:,1);
            ends = intervalsSeconds(:,2);
            if any(~isfinite(starts)) || any(starts < 0)
                error('mpk:Signal:InvalidTimeRanges', ...
                    'Interval start times must be finite and nonnegative.')
            end
            if any(~(isfinite(ends) | isinf(ends))) || any(ends <= 0)
                error('mpk:Signal:InvalidTimeRanges', ...
                    'Interval end times must be finite positive values or Inf.')
            end
            if any(ends <= starts)
                error('mpk:Signal:InvalidTimeRanges', ...
                    'Interval end times must be greater than start times.')
            end
        end
    end
end
