function [figureHandle,axesHandle] = ...
        plotSignalWithPhases(signal,intervalsSeconds,varargin)
%PLOTSIGNALWITHPHASES Plot a Signal with temporal phase intervals.
%
% [FIG,AX] = mpk.utils.plotSignalWithPhases(SIGNAL,INTERVALS) plots the
% complete Signal and overlays 2-by-N half-open intervals in seconds.

if ~isa(signal,'mpk.Signal')
    error('mpk:utils:plotSignalWithPhases:InvalidSignal', ...
        'signal must be an mpk.Signal instance.')
end
validateIntervals(intervalsSeconds)

parser = inputParser;
parser.FunctionName = 'mpk.utils.plotSignalWithPhases';
defaultTitle = 'Signal with analysis phases';
if ~isempty(signal.Label)
    defaultTitle = sprintf('%s with analysis phases',signal.Label);
end
addParameter(parser,'PhaseLabels',{})
addParameter(parser,'Title',defaultTitle, ...
    @(value) ischar(value) || isstring(value))
addParameter(parser,'ShowPhaseLabels',true, ...
    @(value) islogical(value) && isscalar(value))
addParameter(parser,'PhaseAlpha',0.08, ...
    @(value) isnumeric(value) && isscalar(value) && ...
    isfinite(value) && value >= 0 && value <= 1)
addParameter(parser,'FigureVisible',true, ...
    @(value) islogical(value) && isscalar(value))
addParameter(parser,'ParentAxes',[], ...
    @(value) isempty(value) || ...
    (isscalar(value) && isgraphics(value,'axes')))
addParameter(parser,'SignalDisplayName','Signal', ...
    @(value) ischar(value) || isstring(value))
addParameter(parser,'MaxPhaseLabels',25, ...
    @(value) isnumeric(value) && isscalar(value) && ...
    isfinite(value) && value >= 0)
parse(parser,varargin{:})
options = parser.Results;
phaseLabels = normalizeLabels( ...
    options.PhaseLabels,size(intervalsSeconds,2));

if ~isempty(intervalsSeconds)
    if isempty(signal.Time)
        error('mpk:utils:plotSignalWithPhases:PhaseOutsideSignal', ...
            'Temporal phases cannot be plotted on an empty Signal.')
    end
    signalStart = signal.Time(1);
    signalStopExclusive = signal.Time(end) ...
        + 1 / signal.SamplingFrequency;
    if any(intervalsSeconds(1,:) < signalStart) || ...
            any(intervalsSeconds(2,:) > signalStopExclusive)
        error('mpk:utils:plotSignalWithPhases:PhaseOutsideSignal', ...
            'Every phase interval must be fully represented by the Signal.')
    end
end

if isempty(options.ParentAxes)
    if options.FigureVisible
        figureHandle = figure();
    else
        figureHandle = figure('Visible','off');
    end
    axesHandle = axes('Parent',figureHandle);
else
    axesHandle = options.ParentAxes;
    figureHandle = ancestor(axesHandle,'figure');
end

holdWasOn = ishold(axesHandle);
hold(axesHandle,'on')

plot(axesHandle,signal.Time,signal.Samples, ...
    'DisplayName',char(options.SignalDisplayName))

yLimits = ylim(axesHandle);
if yLimits(1) == yLimits(2)
    yLimits = yLimits + [-1 1];
    ylim(axesHandle,yLimits)
end

nPhases = size(intervalsSeconds,2);
phaseColor = [0 0.4470 0.7410];
for phaseIndex = 1:nPhases
    startTime = intervalsSeconds(1,phaseIndex);
    endTime = intervalsSeconds(2,phaseIndex);

    patch(axesHandle, ...
        [startTime endTime endTime startTime], ...
        [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
        phaseColor, ...
        'FaceAlpha',options.PhaseAlpha, ...
        'EdgeColor','none', ...
        'HandleVisibility','off')
    line(axesHandle,[startTime startTime],yLimits, ...
        'LineStyle','--', ...
        'Color',[0.3 0.3 0.3], ...
        'HandleVisibility','off')
    line(axesHandle,[endTime endTime],yLimits, ...
        'LineStyle',':', ...
        'Color',[0.3 0.3 0.3], ...
        'HandleVisibility','off')

    if options.ShowPhaseLabels && ...
            phaseIndex <= options.MaxPhaseLabels && ...
            numel(phaseLabels) >= phaseIndex
        label = phaseLabels{phaseIndex};
        if ~isempty(label)
            text(axesHandle,mean([startTime endTime]),yLimits(2),label, ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','top', ...
                'Interpreter','none', ...
                'FontSize',8, ...
                'HandleVisibility','off')
        end
    end
end

title(axesHandle,char(options.Title),'Interpreter','none')
xlabel(axesHandle,'Time (s)')
ylabelText = 'Signal';
if ~isempty(signal.Units)
    ylabelText = sprintf('Signal (%s)',signal.Units);
end
ylabel(axesHandle,ylabelText,'Interpreter','none')
grid(axesHandle,'on')

if ~holdWasOn
    hold(axesHandle,'off')
end
end

function validateIntervals(intervalsSeconds)
if ~isnumeric(intervalsSeconds) || ~isreal(intervalsSeconds) || ...
        ~ismatrix(intervalsSeconds) || size(intervalsSeconds,1) ~= 2 || ...
        any(~isfinite(intervalsSeconds(:))) || ...
        any(intervalsSeconds(:) < 0) || ...
        any(intervalsSeconds(1,:) >= intervalsSeconds(2,:))
    error('mpk:utils:plotSignalWithPhases:InvalidPhaseIntervals', ...
        ['intervalsSeconds must be a finite nonnegative real 2-by-N ' ...
         'matrix whose starts are smaller than its ends.'])
end
end

function labels = normalizeLabels(values,nPhases)
if isempty(values)
    labels = cell(1,0);
    return
end
if isstring(values)
    values = cellstr(values(:)');
elseif iscell(values)
    values = reshape(values,1,[]);
else
    error('mpk:utils:plotSignalWithPhases:InvalidPhaseLabels', ...
        'PhaseLabels must be a cell array or string array.')
end
if numel(values) ~= nPhases
    error('mpk:utils:plotSignalWithPhases:PhaseLabelMismatch', ...
        'PhaseLabels must contain one label per phase interval.')
end
labels = cell(1,nPhases);
for index = 1:nPhases
    value = values{index};
    if isstring(value) && isscalar(value)
        value = char(value);
    end
    if ~ischar(value) || isempty(value)
        error('mpk:utils:plotSignalWithPhases:InvalidPhaseLabels', ...
            'Every phase label must be a nonempty character vector.')
    end
    labels{index} = value;
end
end
