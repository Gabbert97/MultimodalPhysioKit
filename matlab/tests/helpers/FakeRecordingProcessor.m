classdef FakeRecordingProcessor < handle
    %FAKERECORDINGPROCESSOR Test-only processor that records process calls.

    properties
        Calls = {}
        ReturnValue = struct('Marker','default')
        ReturnCallNumber = false
        ErrorOnCall = Inf
    end

    methods
        function obj = FakeRecordingProcessor(returnValue)
            if nargin >= 1
                obj.ReturnValue = returnValue;
            end
        end

        function result = process(obj,varargin)
            obj.Calls{end + 1} = varargin;
            callNumber = numel(obj.Calls);
            if callNumber == obj.ErrorOnCall
                error('test:FakeRecordingProcessor:ProcessFailed', ...
                    'Fake processor failure on call %d.',callNumber)
            end
            if obj.ReturnCallNumber
                result = struct( ...
                    'CallNumber',callNumber, ...
                    'Features',struct('unchanged',callNumber));
            else
                result = obj.ReturnValue;
            end
        end
    end
end
