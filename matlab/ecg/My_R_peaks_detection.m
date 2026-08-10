function [qrs_amp_raw,qrs_i_raw,f_ECG] = My_R_peaks_detection(ECG,time,fs,gr)
%MY_R_PEAKS_DETECTION Detect ECG R-peaks in one signal segment.
%
%   [QRS_AMP_RAW,QRS_I_RAW,F_ECG] = MY_R_PEAKS_DETECTION(ECG,TIME,FS,GR)
%   filters a raw ECG segment and detects candidate R-peaks using a
%   Pan-Tompkins-style adaptive-threshold procedure independently
%   implemented for MultimodalPhysioKit. Methodological reference:
%
%   J. Pan and W. J. Tompkins, "A Real-Time QRS Detection Algorithm,"
%   IEEE Transactions on Biomedical Engineering, vol. BME-32, no. 3,
%   pp. 230-236, 1985. DOI: 10.1109/TBME.1985.325532
%
%   Inputs
%   ------
%   ECG  : Raw ECG samples for one segment.
%   TIME : Time vector associated with ECG, in seconds.
%   FS   : ECG sampling frequency, in Hz.
%   GR   : Logical flag. If true, diagnostic figures are generated.
%
%   Outputs
%   -------
%   QRS_AMP_RAW : Detected R-peak amplitudes on the filtered ECG.
%   QRS_I_RAW   : Detected R-peak sample indices, relative to ECG.
%   F_ECG       : Filtered ECG used by the detector.
%
%   Notes
%   -----
%   - Peak indices are relative to the input segment, not to the original
%     full recording.
%   - The implementation relies on the project filtering helper.
%   - This function is intended for research preprocessing; visually inspect
%     detections when applying it to new data.

    %% Preprocessing
    % Remove the ECG offset, filter the signal, and compute the
    % derivative/squared/integrated signals used by the detector.
    %% Detector state initialization
    delay = 0;
    skip = 0;                                                                  % becomes one when a T wave is detected
    m_selected_RR = 0;
    mean_RR = 0;
    ser_back = 0; 
    ax = zeros(1,6);

    %% Remove baseline offset
    ECG = ECG - mean(ECG);
    %% ECG filtering
    % Apply the project ECG filter before derivative-based QRS detection.
    f_ECG = filtering(ECG,fs,4,16,2,2,0);
    
    %% ECG differentiation and normalization
    diff_f_ECG = diff(f_ECG);
    n_diff_f_ECG = diff_f_ECG/max(diff_f_ECG);
    
    %% ECG rectification
    squared_ECG = n_diff_f_ECG.^2;

    %% Moving-window integration
    window_length = 0.15;                                                   % Moving window length [s]
    nbr_samples_int_window = window_length*fs;                              % # samples in the integration window [-]
    % Smooth the squared derivative using a 150 ms moving integration window.
    int_ECG = conv(squared_ECG ,ones(1 ,round(nbr_samples_int_window))/round(nbr_samples_int_window));
    delay = round(nbr_samples_int_window)/2;
    
    if gr
        figure('Name','Pan-Tompkins Algorithm')
        title('Pan-Tompkins Algorithm')
        fig(1) = subplot(4,1,1);
        hold on
        plot(time(2:end),ECG(2:end),'LineWidth',1.5)
        plot(time(2:end),f_ECG(2:end),'LineWidth',1.5)
        hold off
        grid on
        ylabel('[mv]')
        legend('Raw ECG','Filtered ECG')
        set(gca,'FontSize',14)
        fig(2) = subplot(4,1,2);
        hold on
        plot(time(2:end),n_diff_f_ECG,'LineWidth',1.5)
        grid on
        ylabel('[-]')
        legend('Diff and normalized ECG')
        set(gca,'FontSize',14)
        fig(3) = subplot(4,1,3);
        plot(time(2:end),squared_ECG,'LineWidth',1.5)
        ylabel('[-]')
        grid on
        legend('Squared ECG')
        set(gca,'FontSize',14)
        fig(4) = subplot(4,1,4);
        plot(time(2:end),int_ECG(1:length(squared_ECG)),'LineWidth',1.5)
        grid on
        xlabel('time [s]')
        ylabel('[-]')
        legend('Integrated ECG')
        linkaxes(fig,'x')
        set(gca,'FontSize',14)
    end
    
    %% Candidate peak detection
    % Candidate peaks are constrained by a minimum physiological RR interval.
    % The detector assumes that valid R-peaks cannot occur less than about
    % 200 ms apart.
    [pks,locs] = findpeaks(int_ECG,'MINPEAKDISTANCE',round(0.2*fs));%round(0.3*fs));
    
    % Remove extreme integrated-signal peaks before adaptive thresholding.
    mean_pks = mean(pks);
    locs = locs(pks <= mean_pks*10);
    pks = pks(pks <= mean_pks*10);

    %% Adaptive-threshold buffers
    LLp = length(pks);
    % Store QRS candidates from the integrated and filtered signals.
    qrs_c = zeros(1,LLp);           % amplitude of R
    qrs_i = zeros(1,LLp);           % index
    qrs_i_raw = zeros(1,LLp);       % amplitude of R
    qrs_amp_raw= zeros(1,LLp);      % Index
    % Store rejected noise candidates for threshold adaptation.
    nois_c = zeros(1,LLp);
    nois_i = zeros(1,LLp);
    % Track adaptive signal, noise, and threshold levels over time.
    SIGL_buf = zeros(1,LLp);
    NOISL_buf = zeros(1,LLp);
    SIGL_buf1 = zeros(1,LLp);
    NOISL_buf1 = zeros(1,LLp);
    THRS_buf1 = zeros(1,LLp);
    THRS_buf = zeros(1,LLp);
    
    
    %% Initialize adaptive thresholds from the first two seconds
    THR_SIG = max(int_ECG(1:2*fs))*1/4;                                          % 0.25 of the max amplitude 
    THR_NOISE = mean(int_ECG(1:2*fs))*1/2;                                       % 0.5 of the mean signal is considered to be noise
    SIG_LEV = THR_SIG;
    NOISE_LEV = THR_NOISE;
    
    
    %% Initialize filtered-signal thresholds from the first two seconds
    THR_SIG1 = max(f_ECG(1:2*fs))*1/4;                                          % 0.25 of the max amplitude 
    THR_NOISE1 = mean(f_ECG(1:2*fs))*1/2; 
    SIG_LEV1 = THR_SIG1;                                                        % Signal level in Bandpassed filter
    NOISE_LEV1 = THR_NOISE1;                                                    % Noise level in Bandpassed filter

    %% Adaptive thresholding and decision rules
    Beat_C = 0;                                                                 % Raw Beats
    Beat_C1 = 0;                                                                % Filtered Beats
    Noise_Count = 0;                                                            % Noise Counter

    for i = 1 : LLp  
       %% Locate the corresponding peak in the filtered ECG
        if locs(i)-round(nbr_samples_int_window)>= 1 && locs(i)<= length(f_ECG)
              [y_i,x_i] = max(f_ECG(locs(i)-round(nbr_samples_int_window):locs(i)));
           else
              if i == 1
                [y_i,x_i] = max(f_ECG(1:locs(i)));
                ser_back = 1;
              elseif locs(i)>= length(f_ECG)
                [y_i,x_i] = max(f_ECG(locs(i)-round(nbr_samples_int_window):end));
              end       
        end  

      %% Update RR-interval statistics
        if Beat_C >= 9        
            diffRR = diff(qrs_i(Beat_C-8:Beat_C));                                   % calculate RR interval
            mean_RR = mean(diffRR);                                                  % calculate the mean of 8 previous R waves interval
            comp =qrs_i(Beat_C)-qrs_i(Beat_C-1);                                     % latest RR
        
            if comp <= 0.92*mean_RR || comp >= 1.16*mean_RR
         % ------ lower down thresholds to detect better in MVI -------- %
                    THR_SIG = 0.5*(THR_SIG);
                    THR_SIG1 = 0.5*(THR_SIG1);               
            else
                m_selected_RR = mean_RR;                                       % The latest regular beats mean
            end        
        end
        
     %% Estimate recent RR interval for missed-beat detection
           if m_selected_RR
               test_m = m_selected_RR;                                         %if the regular RR availabe use it   
           elseif mean_RR && m_selected_RR == 0
               test_m = mean_RR;   
           else
               test_m = 0;
           end
            
        if test_m
              if (locs(i) - qrs_i(Beat_C)) >= round(1.66*test_m)                  % it shows a QRS is missed 
                  [pks_temp,locs_temp] = max(int_ECG(qrs_i(Beat_C)+ round(0.200*fs):locs(i)-round(0.200*fs))); % search back and locate the max in this interval
                  locs_temp = qrs_i(Beat_C)+ round(0.200*fs) + locs_temp -1;      % location 
                 
                  if pks_temp > THR_NOISE
                   Beat_C = Beat_C + 1;
                   qrs_c(Beat_C) = pks_temp;
                   qrs_i(Beat_C) = locs_temp;      
                  % ------------- Locate in Filtered Sig ------------- %
                   if locs_temp <= length(f_ECG)
                      [y_i_t,x_i_t] = max(f_ECG(locs_temp-round(nbr_samples_int_window):locs_temp));
                   else
                      [y_i_t,x_i_t] = max(f_ECG(locs_temp-round(nbr_samples_int_window):end));
                   end
                  % ----------- Band pass Sig Threshold ------------------%
                   if y_i_t > THR_NOISE1 
                      Beat_C1 = Beat_C1 + 1;
                      qrs_i_raw(Beat_C1) = locs_temp-round(nbr_samples_int_window)+ (x_i_t - 1);% save index of bandpass 
                      qrs_amp_raw(Beat_C1) = y_i_t;                               % save amplitude of bandpass 
                      SIG_LEV1 = 0.25*y_i_t + 0.75*SIG_LEV1;                      % when found with the second thres 
                   end
                   
                   not_nois = 1;
                   SIG_LEV = 0.25*pks_temp + 0.75*SIG_LEV ;                       % when found with the second threshold             
                 end             
              else
                  not_nois = 0;         
              end
        end
      
        %% Classify candidate peaks as QRS or noise
        if pks(i) >= THR_SIG      
          % Reject likely T-waves occurring too close to the previous QRS.
           if Beat_C >= 3
              if (locs(i)-qrs_i(Beat_C)) <= round(0.3600*fs)
                  Slope1 = mean(diff(int_ECG(locs(i)-round(0.075*fs):locs(i))));       % mean slope of the waveform at that position
                  Slope2 = mean(diff(int_ECG(qrs_i(Beat_C)-round(0.075*fs):qrs_i(Beat_C)))); % mean slope of previous R wave
                  if abs(Slope1) <= abs(0.5*(Slope2))                              % slope less then 0.5 of previous R
                     Noise_Count = Noise_Count + 1;
                     nois_c(Noise_Count) = pks(i);
                     nois_i(Noise_Count) = locs(i);
                     skip = 1;                                                 % T wave identification
                     % ----- adjust noise levels ------ %
                     NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;
                     NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV; 
                  else
                     skip = 0;
                  end
                
               end
            end
            % If skip is true, the candidate is treated as a T-wave.
            if skip == 0    
              Beat_C = Beat_C + 1;
              qrs_c(Beat_C) = pks(i);
              qrs_i(Beat_C) = locs(i);
            
            % Confirm QRS candidate on the filtered ECG.
              if y_i >= THR_SIG1  
                  Beat_C1 = Beat_C1 + 1;
                  if ser_back 
                     qrs_i_raw(Beat_C1) = x_i;                                 % save index of bandpass 
                  else
                     qrs_i_raw(Beat_C1)= locs(i)-round(nbr_samples_int_window)+ (x_i - 1);   % save index of bandpass 
                  end
                  qrs_amp_raw(Beat_C1) =  y_i;                                 % save amplitude of bandpass 
                  SIG_LEV1 = 0.125*y_i + 0.875*SIG_LEV1;                       % adjust threshold for bandpass filtered sig
              end
             SIG_LEV = 0.125*pks(i) + 0.875*SIG_LEV ;                          % adjust Signal level
            end
                  
        elseif (THR_NOISE <= pks(i)) && (pks(i) < THR_SIG)
             NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;                        % adjust Noise level in filtered sig
             NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV;                       % adjust Noise level in MVI       
        elseif pks(i) < THR_NOISE
            Noise_Count = Noise_Count + 1;
            nois_c(Noise_Count) = pks(i);
            nois_i(Noise_Count) = locs(i);    
            NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;                         % noise level in filtered signal    
            NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV;                        % adjust Noise level in MVI     
        end
                   
        %% Update adaptive thresholds from current signal/noise estimates
        if NOISE_LEV ~= 0 || SIG_LEV ~= 0
            THR_SIG = NOISE_LEV + 0.25*(abs(SIG_LEV - NOISE_LEV));
            THR_NOISE = 0.5*(THR_SIG);
        end
        
        % Update adaptive thresholds for the filtered ECG.
        if NOISE_LEV1 ~= 0 || SIG_LEV1 ~= 0
            THR_SIG1 = NOISE_LEV1 + 0.25*(abs(SIG_LEV1 - NOISE_LEV1));
            THR_NOISE1 = 0.5*(THR_SIG1);
        end
        
        
    % Store threshold history for the integrated signal.
    SIGL_buf(i) = SIG_LEV;
    NOISL_buf(i) = NOISE_LEV;
    THRS_buf(i) = THR_SIG;
    
    % Store threshold history for the filtered ECG.
    SIGL_buf1(i) = SIG_LEV1;
    NOISL_buf1(i) = NOISE_LEV1;
    THRS_buf1(i) = THR_SIG1;
    % Reset per-candidate state before the next iteration.
    skip = 0;                                                   
    not_nois = 0; 
    ser_back = 0;    
    end
    %% Trim preallocated peak arrays
    qrs_i_raw = qrs_i_raw(1:Beat_C1);
    qrs_amp_raw = qrs_amp_raw(1:Beat_C1);
    qrs_c = qrs_c(1:Beat_C);
    qrs_i = qrs_i(1:Beat_C);
    %% Optional diagnostic plots
%     if gr
%       hold on,scatter(qrs_i,qrs_c,'m');
%       hold on,plot(locs,NOISL_buf,'--k','LineWidth',2);
%       hold on,plot(locs,SIGL_buf,'--r','LineWidth',2);
%       hold on,plot(locs,THRS_buf,'--g','LineWidth',2);
%      if any(ax)
%       ax(~ax) = []; 
%       linkaxes(ax,'x');
%       zoom on;
%      end
%     end
    
    %% Overlay detections and thresholds on diagnostic figures
     if gr
       figure;
       az(1)=subplot(311);
       plot(f_ECG);
       title('QRS on Filtered Signal');
       axis tight;
       hold on,scatter(qrs_i_raw,qrs_amp_raw,'m');
       hold on,plot(locs,NOISL_buf1,'LineWidth',2,'Linestyle','--','color','k');
       hold on,plot(locs,SIGL_buf1,'LineWidth',2,'Linestyle','-.','color','r');
       hold on,plot(locs,THRS_buf1,'LineWidth',2,'Linestyle','-.','color','g');
       az(2)=subplot(312);plot(int_ECG);
       title('QRS on MVI signal and Noise level(black),Signal Level (red) and Adaptive Threshold(green)');axis tight;
       hold on,scatter(qrs_i,qrs_c,'m');
       hold on,plot(locs,NOISL_buf,'LineWidth',2,'Linestyle','--','color','k');
       hold on,plot(locs,SIGL_buf,'LineWidth',2,'Linestyle','-.','color','r');
       hold on,plot(locs,THRS_buf,'LineWidth',2,'Linestyle','-.','color','g');
       az(3)=subplot(313);
       plot(ECG-mean(ECG));
       title('Pulse train of the found QRS on ECG signal');
       axis tight;
       line(repmat(qrs_i_raw,[2 1]),...
           repmat([min(ECG-mean(ECG))/2; max(ECG-mean(ECG))/2],size(qrs_i_raw)),...
           'LineWidth',2.5,'LineStyle','-.','Color','r');
       linkaxes(az,'x');
       zoom on;
     end
    

end
