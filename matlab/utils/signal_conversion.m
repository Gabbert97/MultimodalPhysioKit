function [p_data] = signal_conversion(data,att,channels)
%SIGNAL_CONVERSION Convert supported raw BioSignalsPlux channels to units.
%
% Inputs
% ------
% DATA     : N-by-C raw integer sample matrix.
% ATT      : 1-by-C cell array of sensor labels.
% CHANNELS : C-by-2 matrix. Column 1 contains channel IDs; column 2 contains
%            ADC resolution in bits.
%
% Output
% ------
% P_DATA   : Struct with converted fields for recognized sensors:
%            ecg, eda, rip, ppg, temp, fnirs, respiration, and/or events.

    VCC = 3;                                                                % Operating Voltage [V]
    n = channels(:,2);                                                      % n of bits of the channel [-]
    l = length(att);    
    count = 1;
    for i = 1:l
        % Temperature
        if strcmp(att{i},'TEMP')
            a0 = 1.12764514e-3;
            a1 = 2.34282709e-4;
            a2 = 8.77303013e-8;
            NTC_V = data(:,i)*VCC/(2^n(i));                                 % NTC output in Volt [V]
            NTC_Ohm = (NTC_V.*1e+4)./(VCC-NTC_V);                             % NTC resistance in Ohm [Ohm]
            TMP_K = 1./(a0 + a1 * log(NTC_Ohm) + a2 * (log(NTC_Ohm)).^3);   % Temperature [°K]
            p_data.temp = TMP_K - 273.15;                                   % Temperature [°C]
        % Respiration
        elseif strcmp(att{i},'RIP')
            p_data.rip = (data(:,i)./(2^n(i)) - 0.5).*100;                  % Displacement value in percentage [%] of full scale
        % PPG (PLUX BVP sensor; no standard physical unit)
        elseif strcmp(att{i},'BVP')
            p_data.ppg = (data(:,i)./(2^n(i)) - 0.5).*100;                  % Percentage [%] of full scale, centred
        % fNIRS
        elseif strcmp(att{i},'hSpO2')
            p_data.fnirs(:,count) = 0.15*data(:,i)./(2^n(i));               % Photodiode current in microamperes [µA]
            count = count +1;
        % EDA
        elseif strcmp(att{i},'EDA')
            p_data.eda = data(:,i).*VCC./(0.12*2^n(i));                     % EDA value in microsiemens [µS]
        % ECG
        elseif strcmp(att{i},'ECG')
            G_ECG = 1019;                                                   % Sensor gain
            ECG_V = (data(:,i)./(2^n(i))-0.5).*VCC./G_ECG;                  % ECG value in Volt [V]
            p_data.ecg = ECG_V*1e+3;                                        % ECG value in mV [mV] 
        elseif strcmp(att{i},'CUSTOM/0.5/1.0/V')
            events = data(:,i) > 4e+4;
            p_data.events = events;
        elseif strcmp(att{i},'RESPIRATION')
            G_PZT = 1;
            p_data.respiration = (data(:,i)./(2^n(i)-1)-1/2)*VCC/G_PZT;
        end        
    end
end
