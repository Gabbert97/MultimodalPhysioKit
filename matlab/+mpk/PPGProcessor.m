classdef PPGProcessor
    %PPGPROCESSOR Detect pulses and extract features from one PPG Signal.
    %
    % PROCESSOR = mpk.PPGProcessor('LabelFrequencyAnalysis',METHOD, ...
    %     'ReturnIntermediates',FLAG,'ShowPlots',FLAG) configures STFT
    % (METHOD = 1, default), PSD (METHOD = 2), or no HRV frequency analysis
    % (METHOD = [], for signals shorter than 60 s such as real-time windows).
    %
    % RESULT = PROCESSOR.PROCESS(SIGNAL) analyses exactly the supplied
    % mpk.Signal as one piece and returns one mpk.ProcessingResult.
    % Intermediate indexes are one-based.

    properties
        LabelFrequencyAnalysis
        ReturnIntermediates
        ShowPlots
    end

    methods
        function obj = PPGProcessor(varargin)
            parser = inputParser;
            parser.FunctionName = 'mpk.PPGProcessor';
            addParameter(parser,'LabelFrequencyAnalysis',1);
            addParameter(parser,'ReturnIntermediates',false);
            addParameter(parser,'ShowPlots',false);
            parse(parser,varargin{:});
            obj.LabelFrequencyAnalysis = parser.Results.LabelFrequencyAnalysis;
            obj.ReturnIntermediates = parser.Results.ReturnIntermediates;
            obj.ShowPlots = parser.Results.ShowPlots;
            obj.validate();
        end

        function result = process(obj,signal)
            %PROCESS Process exactly one supplied mpk.Signal.
            if ~isa(signal,'mpk.Signal')
                error('mpk:PPGProcessor:InvalidSignal', ...
                    'signal must be an mpk.Signal instance.')
            end
            obj.validate();

            [features,intermediates] = ppg_signal_processing( ...
                signal.Samples,signal.Time,signal.SamplingFrequency, ...
                obj.LabelFrequencyAnalysis,obj.ShowPlots);

            metadata = struct( ...
                'source','process_ppg', ...
                'sampling_frequency',signal.SamplingFrequency, ...
                'n_samples',signal.count(), ...
                'start_time_seconds',signal.Time(1), ...
                'end_time_seconds',signal.Time(end), ...
                'frequency_analysis_method',obj.LabelFrequencyAnalysis, ...
                'n_pulses',numel(intermediates.pulse_valid), ...
                'n_valid_pulses',nnz(intermediates.pulse_valid), ...
                'n_accepted_ibis',numel(intermediates.ibi));
            if ~obj.ReturnIntermediates
                intermediates = struct();
            end
            label = signal.Label;
            if isempty(label)
                label = 'ppg';
            end
            options = struct( ...
                'LabelFrequencyAnalysis',obj.LabelFrequencyAnalysis, ...
                'ReturnIntermediates',obj.ReturnIntermediates, ...
                'ShowPlots',obj.ShowPlots);
            result = mpk.ProcessingResult('ppg','process_ppg',features, ...
                label,options,metadata,intermediates);
        end
    end

    methods (Access = private)
        function validate(obj)
            value = obj.LabelFrequencyAnalysis;
            if ~(isempty(value) || (isnumeric(value) && isscalar(value) && ...
                    any(value == [1 2])))
                error('mpk:PPGProcessor:InvalidFrequencyAnalysis', ...
                    'LabelFrequencyAnalysis must be 1, 2, or [].')
            end
            if ~islogical(obj.ReturnIntermediates) || ~isscalar(obj.ReturnIntermediates)
                error('mpk:PPGProcessor:InvalidReturnIntermediates', ...
                    'ReturnIntermediates must be a logical scalar.')
            end
            if ~islogical(obj.ShowPlots) || ~isscalar(obj.ShowPlots)
                error('mpk:PPGProcessor:InvalidShowPlots', ...
                    'ShowPlots must be a logical scalar.')
            end
        end
    end
end
