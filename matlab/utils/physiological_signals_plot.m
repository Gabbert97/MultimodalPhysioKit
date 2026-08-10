function [] = physiological_signals_plot(data,i_phases,gr)
%PHYSIOLOGICAL_SIGNALS_PLOT Plot converted physiological signals.
%
% Inputs
% ------
% DATA     : Struct containing converted signals plus time and events fields.
% I_PHASES : 2-by-N matrix of one-based inclusive phase limits. Empty is
%            allowed when no phase markers should be drawn.
% GR       : Logical flag. If true, figures are generated.

%% PPlot section
time = data.time;
events = data.events;

if gr
    % ECG
    if isfield(data, 'ecg')
        figure('Name','ECG')
        fig_ecg(1)=subplot(2,1,1);
        plot(time,data.ecg,'LineWidth',1.5)
        hold on
        if ~isempty(i_phases)
            xline(time(i_phases(1,:)),'--k')
            xline(time(i_phases(2,:)),'--k')
        end
        xlim([0 time(end)])
        ylabel('ECG [mV]')
        grid on
        set(gca,'FontSize',14)
        fig_ecg(2)=subplot(2,1,2);
        plot(time,events,'LineWidth',1.5)
        xlim([0 time(end)])
        grid on
        xlabel('time [s]')
        ylabel('Digital I/O')
        linkaxes(fig_ecg,'x')
        set(gca,'FontSize',14)
    end
    % EDA
    if isfield(data, 'eda')
        figure('Name','EDA')
        fig_eda(1)=subplot(2,1,1);
        plot(time,data.eda,'LineWidth',1.5)
        hold on
        if ~isempty(i_phases)
            xline(time(i_phases(1,:)),'--k')
            xline(time(i_phases(2,:)),'--k')
        end
        xlim([0 time(end)])
        ylabel('EDA [\muS]')
        grid on
        set(gca,'FontSize',14)
        fig_eda(2)=subplot(2,1,2);
        plot(time,events,'LineWidth',1.5)
        xlim([0 time(end)])
        grid on
        xlabel('time [s]')
        ylabel('Digital I/O')
        linkaxes(fig_eda,'x')
        set(gca,'FontSize',14)
    end
    if isfield(data, 'temp')
        % Temperature
        figure('Name','Temperature')
        fig_temp(1)=subplot(2,1,1);
        plot(time,data.temp,'LineWidth',1.5)
        hold on
        if ~isempty(i_phases)
            xline(time(i_phases(1,:)),'--k')
            xline(time(i_phases(2,:)),'--k')
        end
        xlim([0 time(end)])
        ylabel('Temperature [°C]')
        grid on
        set(gca,'FontSize',14)
        fig_temp(2)=subplot(2,1,2);
        plot(time,events,'LineWidth',1.5)
        xlim([0 time(end)])
        grid on
        xlabel('time [s]')
        ylabel('Digital I/O')
        linkaxes(fig_temp,'x')
        set(gca,'FontSize',14)
    end
    % fNIRS
    if isfield(data, 'fnirs')
        for j = 1:size(data.fnirs,2)
            figure('Name',strcat('fNIRS_',string(j)))
            fig_fnirs(1)=subplot(2,1,1);
            plot(time,data.fnirs(:,j),'LineWidth',1.5)
            hold on
            if ~isempty(i_phases)
                xline(time(i_phases(1,:)),'--k')
                xline(time(i_phases(2,:)),'--k')
            end
            xlim([0 time(end)])
            ylabel('fNIRS [\muA]')
            grid on
            set(gca,'FontSize',14)
            fig_fnirs(2)=subplot(2,1,2);
            plot(time,events,'LineWidth',1.5)
            xlim([0 time(end)])
            grid on
            xlabel('time [s]')
            ylabel('Digital I/O')
            linkaxes(fig_fnirs,'x')
            set(gca,'FontSize',14)
        end
    end

    % Respiration
    if isfield(data, 'rip')
        figure('Name','Respiration')
        fig_resp(1)=subplot(2,1,1);
        plot(time,data.rip,'LineWidth',1.5)
        hold on
        if ~isempty(i_phases)
            xline(time(i_phases(1,:)),'--k')
            xline(time(i_phases(2,:)),'--k')
        end
        xlim([0 time(end)])
        ylabel('Respiration [%]')
        grid on
        set(gca,'FontSize',14)
        fig_resp(2)=subplot(2,1,2);
        plot(time,events,'LineWidth',1.5)
        xlim([0 time(end)])
        grid on
        xlabel('time [s]')
        ylabel('Digital I/O')
        linkaxes(fig_resp,'x')
        set(gca,'FontSize',14)
    end
end
end
