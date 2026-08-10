function [max_rel,min_rel] = find_max_min(signal,flag_mindistance,percentage_maxdistance) 
%FIND_MAX_MIN Return relative maxima and minima indexes for a vector.
%
% FLAG_MINDISTANCE and PERCENTAGE_MAXDISTANCE are retained for compatibility
% with earlier scripts. The current implementation returns all detected
% turning points.

%% State initialization
up=0;
down=0;
j_min=1;
j_max=1;
min_rel=[];
max_rel=[];

d = signal(2) - signal(1);

if d > 0
    up = 1;
else
    down = 1;
end

for i = 2 : length(signal)-1
    if up == 1
        if signal(i+1) < signal(i)
            down=1;
            up=0;
            max_rel(j_max)=i;
            j_max=j_max+1;
        end
    end
    if down == 1
        if signal(i+1) > signal(i)
            up=1;
            down=0;
            min_rel(j_min)=i;
            j_min=j_min+1;
        end
    end
end

% figure
% hold on
% plot(1:length(signal),signal)
% plot(min_rel,signal(min_rel),'*g')
% plot(max_rel,signal(max_rel),'*r')

% if length(min_rel) > 1    
%     if flag_mindistance == 'MinPeakDistance'
%         peaks_distance = diff(min_rel);
%         mean_peak_distance = mean(peaks_distance);
%         for i = 1 : length(peaks_distance)
%             if peaks_distance(i) >= percentage_maxdistance * mean_peak_distance
%                 min_def(count) = min_rel(i+1);
%                 count = count + 1;
%             end
%         end
%     else
%         min_def = min_rel(2:end);
%     end
% else
%     min_def = 0;
% end

% figure;
% hold on
% plot(1:length(signal),signal)
% plot(min_def,signal(min_def),'o')
% hold off

end
