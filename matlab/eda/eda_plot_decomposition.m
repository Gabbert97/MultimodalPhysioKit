function eda_plot_decomposition(decomposition,options)
%EDA_PLOT_DECOMPOSITION Optionally plot or save EDA diagnostics.

if ~(options.ShowPlots || options.SavePlots)
    return
end
figureHandle = figure('Name','EDA decomposition');
cleanup = onCleanup(@() localCloseIfHidden( ...
    figureHandle,options.ShowPlots));
axesHandle = axes(figureHandle);
hold(axesHandle,'on')
plot(axesHandle,decomposition.Time, ...
    decomposition.EdaFiltered,'LineWidth',1.5,'DisplayName','eda')
plot(axesHandle,decomposition.Time, ...
    decomposition.Scl,'LineWidth',1.5,'DisplayName','scl')
plot(axesHandle,decomposition.Time, ...
    decomposition.Scr,'LineWidth',1.5,'DisplayName','scr')
plot(axesHandle,decomposition.Time, ...
    decomposition.Residuals,'LineWidth',1.5,'DisplayName','residuals')
if isfield(options,'MinimumScrAmplitude')
    [~,~,~,acceptedMaxima,~] = scr_processing( ...
        decomposition.Scr,decomposition.Time, ...
        decomposition.SamplingFrequency,options.MinimumScrAmplitude);
    if ~isempty(acceptedMaxima)
        scrFiltered = filtering( ...
            decomposition.Scr,decomposition.SamplingFrequency,0,2,0,2,0);
        plot(axesHandle,decomposition.Time(acceptedMaxima), ...
            scrFiltered(acceptedMaxima),'xk', ...
            'DisplayName','accepted SCR')
    end
end
grid(axesHandle,'on')
legend(axesHandle,'show')
xlabel(axesHandle,'Time (s)')
ylabel(axesHandle,'EDA (\muS)')
title(axesHandle,options.OutputLabel,'Interpreter','none')

if options.SavePlots
    if ~isfolder(options.OutputDirectory)
        mkdir(options.OutputDirectory)
    end
    outputFile = fullfile(options.OutputDirectory, ...
        sprintf('%s_eda_decomposition.png',options.OutputLabel));
    print(figureHandle,outputFile,'-dpng','-r150')
end
if options.ShowPlots
    drawnow
end
end

function localCloseIfHidden(figureHandle,showPlots)
if ~showPlots && isgraphics(figureHandle,'figure')
    close(figureHandle)
end
end
