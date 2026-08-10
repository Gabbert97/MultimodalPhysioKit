function varargout = fake_cvxeda(varargin)
%FAKE_CVXEDA Deterministic test-only seven-output cvxEDA replacement.
persistent state
if isempty(state)
    state = localDefaultState();
end

if nargin == 1 && (ischar(varargin{1}) || isstring(varargin{1}))
    command = char(varargin{1});
    switch command
        case 'reset'
            state = localDefaultState();
            varargout = {};
        case 'state'
            varargout{1} = state;
        case 'short_output'
            state.OutputLengthAdjustment = -1;
            varargout = {};
        case 'nonfinite_output'
            state.ReturnNonfinite = true;
            varargout = {};
        otherwise
            error('test:FakeCvxEDA:UnknownCommand', ...
                'Unknown fake cvxEDA command.')
    end
    return
end

edaFiltered = varargin{1};
delta = varargin{2};
state.CallCount = state.CallCount + 1;
state.LastInput = edaFiltered(:);
state.LastDelta = delta;

outputLength = numel(edaFiltered) + state.OutputLengthAdjustment;
indexes = (0:(outputLength-1))';
time = indexes .* delta;
scr = 0.2 + 0.1*sin(2*pi*0.2*time);
scl = 1.5 + 0.0002*indexes;
residuals = 0.001*cos(2*pi*0.1*time);
if state.ReturnNonfinite && ~isempty(scr)
    scr(1) = NaN;
end
unused = zeros(outputLength,1);
varargout = {scr,unused,scl,unused,unused,residuals,unused};
end

function state = localDefaultState()
state = struct( ...
    'CallCount',0, ...
    'LastInput',[], ...
    'LastDelta',[], ...
    'OutputLengthAdjustment',0, ...
    'ReturnNonfinite',false);
end
