classdef ProcessingResult
    %PROCESSINGRESULT Feature output and processing provenance.
    %
    % RESULT = mpk.ProcessingResult(MODALITY,PROCESSORNAME,FEATURES, ...)
    % records processor identity, source signal label, effective options,
    % non-sensitive metadata, and a compact creation timestamp.
    % FEATURES remains in the format returned by the wrapped procedural API.

    properties
        Modality
        ProcessorName
        Features
        SourceSignalLabel
        Options
        Metadata
        Intermediates
        CreatedOn
    end

    methods
        function obj = ProcessingResult( ...
                modality,processorName,features,sourceSignalLabel, ...
                options,metadata,intermediates)
            if nargin < 4 || isempty(sourceSignalLabel)
                sourceSignalLabel = '';
            end
            if nargin < 5 || isempty(options)
                options = struct();
            end
            if nargin < 6 || isempty(metadata)
                metadata = struct();
            end
            if nargin < 7 || isempty(intermediates)
                intermediates = struct();
            elseif ~isstruct(intermediates) || ~isscalar(intermediates)
                error('mpk:ProcessingResult:InvalidIntermediates', ...
                    'intermediates must be a scalar struct.')
            end

            obj.Modality = char(modality);
            obj.ProcessorName = char(processorName);
            obj.Features = features;
            obj.SourceSignalLabel = char(sourceSignalLabel);
            obj.Options = options;
            obj.Metadata = metadata;
            obj.Intermediates = intermediates;
            createdOn = datetime('now','Format',"yyyyMMdd'T'HHmmss");
            obj.CreatedOn = char(createdOn);
        end

        function features = toCell(obj)
            %TOCELL Return the wrapped feature value unchanged.
            features = obj.Features;
        end

        function names = featureNames(obj)
            %FEATURENAMES Return feature names for one or more feature groups.
            if iscell(obj.Features) && size(obj.Features,2) >= 1
                names = obj.Features(:,1);
            elseif isstring(obj.Features) && size(obj.Features,2) >= 1
                names = cellstr(obj.Features(:,1));
            elseif isstruct(obj.Features)
                names = struct();
                groups = fieldnames(obj.Features);
                for index = 1:numel(groups)
                    group = groups{index};
                    values = obj.Features.(group);
                    if iscell(values) && size(values,2) >= 1
                        names.(group) = values(:,1);
                    elseif isstring(values) && size(values,2) >= 1
                        names.(group) = cellstr(values(:,1));
                    else
                        names.(group) = {};
                    end
                end
            else
                names = {};
            end
        end
    end
end
