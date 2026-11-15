function [pnn50,pnn50_max] = ecg_pnn50_features(ibi_def_phase,control_ibi_flag)
%ECG_PNN50_FEATURES pNN50 and longest-run pNN50 from screened IBIs.
%
% Inputs
% ------
% IBI_DEF_PHASE    : Accepted inter-beat intervals [s], in the order of the
%                    nonzero entries of CONTROL_IBI_FLAG.
% CONTROL_IBI_FLAG : Per-interval acceptance flags (nonzero = accepted).
%
% Outputs
% -------
% PNN50     : Share of successive differences between adjacent accepted
%             IBIs that exceed 50 ms (Task Force, 1996). Pairs spanning a
%             rejected interval are not evaluated. Zero if none evaluated.
% PNN50_MAX : Longest run of consecutive differences exceeding 50 ms,
%             divided by the expected number of IBIs in two minutes.

ibi_def_phase = ibi_def_phase(:);
accepted = control_ibi_flag(:) ~= 0;
accepted_count = cumsum(accepted);
n_pairs = max(numel(accepted) - 1,0);
differences = nan(n_pairs,1);
for j = 1:n_pairs
    if accepted(j) && accepted(j+1)
        m = accepted_count(j);
        if m + 1 <= numel(ibi_def_phase)
            differences(j) = abs(ibi_def_phase(m+1) - ibi_def_phase(m));
        end
    end
end

evaluated = ~isnan(differences);
exceeds = evaluated & differences > 0.05;
if any(evaluated)
    pnn50 = sum(exceeds) / sum(evaluated);
else
    pnn50 = 0;
end

longest_run = 0;
current_run = 0;
for j = 1:numel(exceeds)
    if exceeds(j)
        current_run = current_run + 1;
    else
        current_run = 0;
    end
    longest_run = max(longest_run,current_run);
end
expected_ibis_in_2_minutes = 120 * numel(ibi_def_phase) / sum(ibi_def_phase);
pnn50_max = longest_run / expected_ibis_in_2_minutes;
end
