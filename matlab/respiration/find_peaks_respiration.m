function [min_rel,max_rel] = find_peaks_respiration(signal,flag_mindistance,percentage_distance) 

up=0;
down=0;
j_min=1;
j_max=1;
min_rel=[];
max_rel=[];
count_min = 1;
count_max = 1;

d = signal(1) - signal(2);

if d > 0
    up = 1;
else
    down = 1;
end

for i = 1 : length(signal) - 1
    if up == 1
        if signal(i+1) < signal(i)
            down=1;
            up=0;
            max_rel(j_max)=i;
            j_max=j_max+1;
        end
    end
    if down==1
        if signal(i+1) > signal(i)
            up=1;
            down=0;
            min_rel(j_min)=i;
            j_min=j_min+1;
        end
    end
end

% %% Minimum distance
% % Minimum distance check algorithm
% if length(min_rel) > 1    
%     if flag_mindistance == 'MinPeakDistance'
%         peaks_distance = diff(min_rel);
%         mean_peak_distance = mean(peaks_distance);
%         for i = 1 : length(peaks_distance)
%             if peaks_distance(i) >= percentage_distance * mean_peak_distance
%                 min_def(count_min) = min_rel(i+1);
%                 count_min = count_min + 1;
%             end
%         end
%     else
%         min_def = min_rel;
%     end
% elseif min_rel ~= 1
%     min_def = min_rel;
% else
%     min_def = [];
% end
% if length(max_rel) > 1    
%     if flag_mindistance == 'MinPeakDistance'
%         peaks_distance = diff(max_rel);
%         mean_peak_distance = mean(peaks_distance);
%         for i = 1 : length(peaks_distance)
%             if peaks_distance(i) >= percentage_distance * mean_peak_distance
%                 max_def(count_max) = max_rel(i+1);
%                 count_max = count_max + 1;
%             end
%         end
%     else
%         max_def = max_rel;
%     end
% elseif max_rel ~= 1
%     max_def = max_rel;
% else
%     max_def = 0;
% end

% Plot figure
% figure;
% hold on
% plot(1:length(signal),signal)
% plot(min_def,signal(min_def),'o')
% plot(max_def,signal(max_def),'*')
% hold off
end
