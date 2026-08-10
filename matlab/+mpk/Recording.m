classdef Recording < handle
    %RECORDING Named collection of physiological signals.
    %
    % RECORDING = mpk.Recording(ID,METADATA) creates a handle object for one
    % logical recording or session. Signals are added and retrieved by name.
    % The class does not read files or retain acquisition paths.

    properties
        Id
        Signals
        Metadata
        PhaseIntervalsSeconds
        PhaseLabels
    end

    methods
        function obj = Recording(id,metadata)
            if nargin < 1 || isempty(id)
                id = '';
            end
            if nargin < 2 || isempty(metadata)
                metadata = struct();
            end
            if ~isstruct(metadata)
                error('mpk:Recording:InvalidMetadata', ...
                    'Metadata must be a struct.')
            end

            obj.Id = char(id);
            obj.Signals = struct();
            obj.Metadata = metadata;
            obj.PhaseIntervalsSeconds = zeros(2,0);
            obj.PhaseLabels = cell(1,0);
        end

        function addSignal(obj,name,signal)
            %ADDSIGNAL Add or replace a named mpk.Signal.
            if ~isa(signal,'mpk.Signal')
                error('mpk:Recording:InvalidSignal', ...
                    'signal must be an mpk.Signal instance.')
            end

            fieldName = matlab.lang.makeValidName(char(name));
            obj.Signals.(fieldName) = signal;
        end

        function signal = getSignal(obj,name)
            %GETSIGNAL Return a named signal or raise MissingSignal.
            fieldName = matlab.lang.makeValidName(char(name));
            if ~isfield(obj.Signals,fieldName)
                error('mpk:Recording:MissingSignal', ...
                    'Recording does not contain signal "%s".',char(name))
            end
            signal = obj.Signals.(fieldName);
        end

        function tf = hasSignal(obj,name)
            %HASSIGNAL Return true when a named signal is present.
            fieldName = matlab.lang.makeValidName(char(name));
            tf = isfield(obj.Signals,fieldName);
        end

        function names = signalNames(obj)
            %SIGNALNAMES Return stored signal names as a cell array.
            names = fieldnames(obj.Signals);
        end

        function setPhases(obj,intervalsSeconds,labels)
            %SETPHASES Store half-open temporal phase intervals in seconds.
            if nargin < 3
                labels = [];
                labelsSupplied = false;
            else
                labelsSupplied = true;
            end

            mpk.Recording.validatePhaseIntervals(intervalsSeconds);
            nPhases = size(intervalsSeconds,2);
            if labelsSupplied
                normalizedLabels = mpk.Recording.normalizePhaseLabels( ...
                    labels,nPhases);
            else
                normalizedLabels = cell(1,nPhases);
                for index = 1:nPhases
                    normalizedLabels{index} = sprintf('phase_%d',index);
                end
            end

            obj.PhaseIntervalsSeconds = double(intervalsSeconds);
            obj.PhaseLabels = normalizedLabels;
        end

        function clearPhases(obj)
            %CLEARPHASES Remove all temporal phase definitions.
            obj.PhaseIntervalsSeconds = zeros(2,0);
            obj.PhaseLabels = cell(1,0);
        end

        function n = numberOfPhases(obj)
            %NUMBEROFPHASES Return the number of temporal phase columns.
            n = size(obj.PhaseIntervalsSeconds,2);
        end

        function durations = phaseDurationsSeconds(obj)
            %PHASEDURATIONSSECONDS Return phase durations in supplied order.
            durations = obj.PhaseIntervalsSeconds(2,:) ...
                - obj.PhaseIntervalsSeconds(1,:);
        end

        function indices = phaseIndices(obj,signalOrName,phases)
            %PHASEINDICES Map temporal phases through one Signal.Time vector.
            %
            % Row 1 contains one-based inclusive start indexes. Row 2
            % contains one-based exclusive stop indexes. Slice a phase with
            % startIndex:(stopIndexExclusive - 1). The exclusive stop may be
            % signal.count() + 1.
            if obj.numberOfPhases() == 0
                error('mpk:Recording:NoPhases', ...
                    'Recording has no temporal phases.')
            end

            signal = mpk.Recording.resolveSignal(obj,signalOrName);
            if nargin < 3
                phases = [];
            end
            [~,selectedIndexes] = mpk.Recording.selectPhaseLabels(obj,phases);

            mpk.Recording.validateSignalTime(signal);
            signalStart = signal.Time(1);
            signalStopExclusive = signal.Time(end) ...
                + 1 / signal.SamplingFrequency;
            intervals = obj.PhaseIntervalsSeconds(:,selectedIndexes);
            indices = zeros(2,numel(selectedIndexes));

            for index = 1:numel(selectedIndexes)
                startSeconds = intervals(1,index);
                endSeconds = intervals(2,index);
                if startSeconds < signalStart || endSeconds > signalStopExclusive
                    error('mpk:Recording:PhaseOutsideSignal', ...
                        ['Phase "%s" must be fully represented by the ' ...
                        'Signal temporal range [%g,%g].'], ...
                        obj.PhaseLabels{selectedIndexes(index)}, ...
                        signalStart,signalStopExclusive)
                end

                startIndex = find(signal.Time >= startSeconds,1,'first');
                stopIndexExclusive = find(signal.Time >= endSeconds,1,'first');
                if isempty(stopIndexExclusive) && endSeconds <= signalStopExclusive
                    stopIndexExclusive = signal.count() + 1;
                end
                if isempty(startIndex) || isempty(stopIndexExclusive) || ...
                        startIndex >= stopIndexExclusive
                    error('mpk:Recording:EmptyPhase', ...
                        'Phase "%s" maps to no Signal samples.', ...
                        obj.PhaseLabels{selectedIndexes(index)})
                end

                indices(:,index) = [startIndex; stopIndexExclusive];
            end
        end

        function signal = extractPhase(obj,signalName,phaseLabel,varargin)
            %EXTRACTPHASE Return one stored Signal restricted to one phase.
            parser = inputParser;
            parser.FunctionName = 'mpk.Recording.extractPhase';
            addParameter(parser,'ResetTime',true);
            parse(parser,varargin{:});
            resetTime = parser.Results.ResetTime;
            if ~islogical(resetTime) || ~isscalar(resetTime)
                error('mpk:Recording:InvalidResetTime', ...
                    'ResetTime must be a logical scalar.')
            end

            selectedLabel = mpk.Recording.normalizeSingleLabel( ...
                phaseLabel,'phaseLabel');
            normalizedSignalName = mpk.Recording.normalizeSingleLabel( ...
                signalName,'signalName');
            sourceSignal = obj.getSignal(normalizedSignalName);
            mapped = obj.phaseIndices(normalizedSignalName,selectedLabel);
            phaseColumn = find(strcmp(obj.PhaseLabels,selectedLabel),1);
            sampleIndexes = mapped(1,1):(mapped(2,1)-1);
            selectedSamples = sourceSignal.Samples(sampleIndexes);
            selectedTime = sourceSignal.Time(sampleIndexes);
            if resetTime
                selectedTime = selectedTime - selectedTime(1);
            end

            metadata = sourceSignal.Metadata;
            metadata.SourceRecording = obj.Id;
            metadata.PhaseLabel = selectedLabel;
            metadata.PhaseStartSeconds = ...
                obj.PhaseIntervalsSeconds(1,phaseColumn);
            metadata.PhaseEndSeconds = ...
                obj.PhaseIntervalsSeconds(2,phaseColumn);
            metadata.TimeReset = resetTime;
            signal = mpk.Signal(selectedSamples,selectedTime, ...
                sourceSignal.SamplingFrequency,sourceSignal.Label, ...
                sourceSignal.Units,metadata);
        end

        function result = processSignal(obj,processor,signalNames,varargin)
            %PROCESSSIGNAL Process complete stored Signals through PROCESSOR.
            parser = inputParser;
            parser.FunctionName = 'mpk.Recording.processSignal';
            addParameter(parser,'ProcessArguments',{});
            parse(parser,varargin{:});

            mpk.Recording.validateProcessor(processor);
            processArguments = mpk.Recording.normalizeProcessArguments( ...
                parser.Results.ProcessArguments);
            signals = mpk.Recording.resolveNamedSignals(obj,signalNames);
            result = processor.process(signals{:},processArguments{:});
        end

        function results = processPhases(obj,processor,signalNames,varargin)
            %PROCESSPHASES Extract and process each selected temporal phase.
            parser = inputParser;
            parser.FunctionName = 'mpk.Recording.processPhases';
            addParameter(parser,'Phases',[]);
            addParameter(parser,'ResetTime',true);
            addParameter(parser,'ProcessArguments',{});
            parse(parser,varargin{:});

            selectedLabels = mpk.Recording.selectPhaseLabels( ...
                obj,parser.Results.Phases);
            mpk.Recording.validateProcessor(processor);
            resetTime = parser.Results.ResetTime;
            if ~islogical(resetTime) || ~isscalar(resetTime)
                error('mpk:Recording:InvalidResetTime', ...
                    'ResetTime must be a logical scalar.')
            end
            processArguments = mpk.Recording.normalizeProcessArguments( ...
                parser.Results.ProcessArguments);
            normalizedNames = mpk.Recording.normalizeSignalNames(signalNames);
            mpk.Recording.resolveNamedSignals(obj,normalizedNames);

            results = containers.Map('KeyType','char','ValueType','any');
            for phaseIndex = 1:numel(selectedLabels)
                phaseLabel = selectedLabels{phaseIndex};
                phaseSignals = cell(1,numel(normalizedNames));
                for signalIndex = 1:numel(normalizedNames)
                    phaseSignals{signalIndex} = obj.extractPhase( ...
                        normalizedNames{signalIndex},phaseLabel, ...
                        'ResetTime',resetTime);
                end
                result = processor.process( ...
                    phaseSignals{:},processArguments{:});
                results(phaseLabel) = result;
            end
        end

        function out = removeTimeRanges(obj,intervalsSeconds,varargin)
            %REMOVETIMERANGES Return a new recording with time ranges removed.
            if obj.numberOfPhases() > 0
                error('mpk:Recording:PhasesMustBeCleared', ...
                    ['Phases must be cleared or redefined before removing ' ...
                    'time ranges.'])
            end

            parser = inputParser;
            parser.FunctionName = 'mpk.Recording.removeTimeRanges';
            addParameter(parser,'ResetTime',true, ...
                @(value) islogical(value) && isscalar(value));
            parse(parser,varargin{:});
            resetTime = logical(parser.Results.ResetTime);

            mpk.Recording.validateTimeRanges(intervalsSeconds);

            metadata = obj.Metadata;
            metadata.TimeRangesRemovedSeconds = intervalsSeconds;
            metadata.ResetTimeAfterRemoval = resetTime;
            out = mpk.Recording(obj.Id,metadata);

            names = fieldnames(obj.Signals);
            for index = 1:numel(names)
                signalName = names{index};
                out.addSignal(signalName,obj.Signals.(signalName).removeTimeRanges( ...
                    intervalsSeconds,'ResetTime',resetTime));
            end
        end
    end

    methods (Static, Access = private)
        function validatePhaseIntervals(intervalsSeconds)
            if ~isnumeric(intervalsSeconds) || ~isreal(intervalsSeconds) || ...
                    ndims(intervalsSeconds) ~= 2 || ...
                    size(intervalsSeconds,1) ~= 2
                error('mpk:Recording:InvalidPhaseIntervals', ...
                    'intervalsSeconds must be a 2-by-N numeric matrix.')
            end
            if any(~isfinite(intervalsSeconds(:)))
                error('mpk:Recording:InvalidPhaseIntervals', ...
                    'Phase interval values must be finite.')
            end
            if any(intervalsSeconds(:) < 0)
                error('mpk:Recording:InvalidPhaseIntervals', ...
                    'Phase interval values must be nonnegative.')
            end
            if any(intervalsSeconds(1,:) >= intervalsSeconds(2,:))
                error('mpk:Recording:InvalidPhaseIntervals', ...
                    'Every phase start must be strictly smaller than its end.')
            end
        end

        function labels = normalizePhaseLabels(rawLabels,nPhases)
            if isstring(rawLabels)
                if any(ismissing(rawLabels(:)) | strlength(rawLabels(:)) == 0)
                    error('mpk:Recording:EmptyPhaseLabel', ...
                        'Phase labels must be nonempty.')
                end
                labels = cellstr(reshape(rawLabels,1,[]));
            elseif iscell(rawLabels)
                if ~all(cellfun(@(value) ischar(value) && ...
                        (isrow(value) || isempty(value)),rawLabels))
                    error('mpk:Recording:InvalidPhaseLabels', ...
                        ['labels must be a cell array of character vectors ' ...
                        'or a string array.'])
                end
                labels = reshape(rawLabels,1,[]);
            else
                error('mpk:Recording:InvalidPhaseLabels', ...
                    ['labels must be a cell array of character vectors ' ...
                    'or a string array.'])
            end

            if numel(labels) ~= nPhases
                error('mpk:Recording:PhaseLabelMismatch', ...
                    'labels must contain exactly one label per phase.')
            end
            if any(cellfun(@isempty,labels))
                error('mpk:Recording:EmptyPhaseLabel', ...
                    'Phase labels must be nonempty.')
            end
            if numel(unique(labels,'stable')) ~= numel(labels)
                error('mpk:Recording:DuplicatePhaseLabels', ...
                    'Phase labels must be unique.')
            end
        end

        function labels = normalizePhaseSelection(phases)
            if ischar(phases)
                if ~isrow(phases) || isempty(phases)
                    error('mpk:Recording:InvalidPhaseSelection', ...
                        'Phase selection labels must be nonempty.')
                end
                labels = {phases};
            elseif isstring(phases)
                if isempty(phases)
                    error('mpk:Recording:InvalidPhaseSelection', ...
                        'Explicit phase selections must be nonempty.')
                end
                if any(ismissing(phases(:)) | strlength(phases(:)) == 0)
                    error('mpk:Recording:InvalidPhaseSelection', ...
                        'Phase selection labels must be nonempty.')
                end
                labels = cellstr(reshape(phases,1,[]));
            elseif iscell(phases)
                if isempty(phases)
                    error('mpk:Recording:InvalidPhaseSelection', ...
                        'Explicit phase selections must be nonempty.')
                end
                if ~all(cellfun(@(value) ischar(value) && ...
                        isrow(value) && ~isempty(value),phases))
                    error('mpk:Recording:InvalidPhaseSelection', ...
                        'Phase selection must contain character-vector labels.')
                end
                labels = reshape(phases,1,[]);
            else
                error('mpk:Recording:InvalidPhaseSelection', ...
                    ['phases must be a label, a cell array of labels, ' ...
                    'or a string array.'])
            end
            if any(cellfun(@isempty,labels))
                error('mpk:Recording:InvalidPhaseSelection', ...
                    'Phase selection labels must be nonempty.')
            end
        end

        function [labels,indexes] = selectPhaseLabels(obj,phases)
            if obj.numberOfPhases() == 0
                error('mpk:Recording:NoPhases', ...
                    'Recording has no temporal phases.')
            end
            if isnumeric(phases) && isempty(phases)
                labels = obj.PhaseLabels;
                indexes = 1:obj.numberOfPhases();
                return
            end

            labels = mpk.Recording.normalizePhaseSelection(phases);
            indexes = zeros(1,numel(labels));
            for index = 1:numel(labels)
                match = find(strcmp(obj.PhaseLabels,labels{index}),1);
                if isempty(match)
                    error('mpk:Recording:UnknownPhase', ...
                        'Recording does not contain phase "%s".',labels{index})
                end
                indexes(index) = match;
            end
        end

        function validateProcessor(processor)
            try
                valid = ismethod(processor,'process');
            catch
                valid = false;
            end
            if ~valid
                error('mpk:Recording:InvalidProcessor', ...
                    'processor must expose a callable public process method.')
            end
        end

        function names = normalizeSignalNames(signalNames)
            if ischar(signalNames)
                if ~isrow(signalNames) || isempty(signalNames)
                    error('mpk:Recording:InvalidSignalNames', ...
                        'Signal names must be nonempty character vectors.')
                end
                names = {signalNames};
            elseif isstring(signalNames)
                if isempty(signalNames) || any( ...
                        ismissing(signalNames(:)) | ...
                        strlength(signalNames(:)) == 0)
                    error('mpk:Recording:InvalidSignalNames', ...
                        'Signal name collections must be nonempty.')
                end
                names = cellstr(reshape(signalNames,1,[]));
            elseif iscell(signalNames)
                if isempty(signalNames)
                    error('mpk:Recording:InvalidSignalNames', ...
                        'Signal name collections must be nonempty.')
                end
                if ~all(cellfun(@(value) ischar(value) && ...
                        isrow(value) && ~isempty(value),signalNames))
                    error('mpk:Recording:InvalidSignalNames', ...
                        'Signal names must be nonempty character vectors.')
                end
                names = reshape(signalNames,1,[]);
            else
                error('mpk:Recording:InvalidSignalNames', ...
                    ['signalNames must be a character vector, string array, ' ...
                    'or cell array of character vectors.'])
            end
            if numel(unique(names,'stable')) ~= numel(names)
                error('mpk:Recording:DuplicateSignalNames', ...
                    'Signal names must be unique.')
            end
        end

        function signals = resolveNamedSignals(obj,signalNames)
            names = mpk.Recording.normalizeSignalNames(signalNames);
            signals = cell(1,numel(names));
            for index = 1:numel(names)
                signals{index} = obj.getSignal(names{index});
            end
        end

        function arguments = normalizeProcessArguments(processArguments)
            if ~iscell(processArguments)
                error('mpk:Recording:InvalidProcessArguments', ...
                    'ProcessArguments must be a cell array.')
            end
            arguments = processArguments;
        end

        function label = normalizeSingleLabel(rawLabel,argumentName)
            if ischar(rawLabel) && isrow(rawLabel) && ~isempty(rawLabel)
                label = rawLabel;
            elseif isstring(rawLabel) && isscalar(rawLabel) && ...
                    ~ismissing(rawLabel) && strlength(rawLabel) > 0
                label = char(rawLabel);
            else
                error('mpk:Recording:InvalidPhaseSelection', ...
                    '%s must identify exactly one nonempty phase label.', ...
                    argumentName)
            end
        end

        function signal = resolveSignal(obj,signalOrName)
            if isa(signalOrName,'mpk.Signal')
                signal = signalOrName;
                return
            end
            signalName = mpk.Recording.normalizeSingleLabel( ...
                signalOrName,'signalOrName');
            signal = obj.getSignal(signalName);
        end

        function validateSignalTime(signal)
            if signal.count() == 0 || isempty(signal.Time)
                error('mpk:Recording:EmptySignal', ...
                    'Signal must contain at least one sample.')
            end
            if any(~isfinite(signal.Time))
                error('mpk:Recording:InvalidSignalTime', ...
                    'Signal time values must be finite.')
            end
            if any(diff(signal.Time) <= 0)
                error('mpk:Recording:InvalidSignalTime', ...
                    'Signal time values must be strictly increasing.')
            end
        end

        function validateTimeRanges(intervalsSeconds)
            if ~isnumeric(intervalsSeconds) || size(intervalsSeconds,2) ~= 2 || ...
                    isempty(intervalsSeconds)
                error('mpk:Recording:InvalidTimeRanges', ...
                    'intervalsSeconds must be a nonempty Nx2 numeric matrix.')
            end
            starts = intervalsSeconds(:,1);
            ends = intervalsSeconds(:,2);
            if any(~isfinite(starts)) || any(starts < 0)
                error('mpk:Recording:InvalidTimeRanges', ...
                    'Interval start times must be finite and nonnegative.')
            end
            if any(~(isfinite(ends) | isinf(ends))) || any(ends <= 0)
                error('mpk:Recording:InvalidTimeRanges', ...
                    'Interval end times must be finite positive values or Inf.')
            end
            if any(ends <= starts)
                error('mpk:Recording:InvalidTimeRanges', ...
                    'Interval end times must be greater than start times.')
            end
        end
    end
end
