function [deltaC_oxy,deltaC_deoxy,baseline] = LambertBeer( ...
        fnirs_red,fnirs_infrared,time,fs,i_phases,age_years,gr,varargin)
%LAMBERTBEER Convert dual-wavelength fNIRS current to concentration change.
%
% The modified Beer-Lambert implementation, wavelength selection, extinction
% coefficients, source-detector distance, and age-dependent differential
% pathlength factor follow the fNIRS processing method reported in:
%
% Luzzani, G., Pogliano, M., Buraioli, I., Colavincenzo, M., Martorana, S.,
% Guglieri, G., and Demarchi, D. (2025). "ECG, Respiration, fNIRS, and Eye
% Tracking for Stress and Mental Workload Monitoring in Human-Machine
% Interaction." IEEE Access, 13, 122726-122741.
% DOI: 10.1109/ACCESS.2025.3588384
%
% Baseline reference
% ------------------
% This function uses a separate mean current for each wavelength. The
% baseline is computed from the first provided phase, excluding the first
% and last 10 seconds of that phase. Baseline selection is
% protocol-dependent; users should review and, when necessary, adapt the
% indexes below before analysis. See docs/fnirs.md for the full processing
% order and equations.

if gr
    figure
    plot(time,fnirs_red,'r')
    hold on
    plot(time,fnirs_infrared,'b')
    xline(time(i_phases(1,:)),'--')
    xline(time(i_phases(2,:)),'--')
    hold off
    legend('Red emitter - \lambda_1','Infrared emitter - \lambda_2')
    xlabel('Time [s]')
    ylabel('fNIRS [\muA]')
    title('Filtered fNIRS signals')
end

baseline = fnirs_resolve_baseline( ...
    fnirs_red,fnirs_infrared,fs,i_phases,varargin{:});
deltaOD_red = log10(baseline.RedValue./fnirs_red);
deltaOD_infrared = log10(baseline.InfraredValue./fnirs_infrared);

% Extinction coefficients for the two wavelengths and chromophores, source-
% detector distance, and wavelengths used by the cited fNIRS workflow.
epsilon = [3.4408 0.3346 0.7977 1.2071];
source_detector_distance_cm = 2;
lambda1 = 660;
lambda2 = 860;

% Age-dependent differential pathlength factor used in the cited workflow.
DPF_lambda1 = 223.3 + 0.05624*age_years^0.8493 ...
    - 5.723e-7*lambda1^3 + 0.001245*lambda1^2 - 0.9025*lambda1;
DPF_lambda2 = 223.3 + 0.05624*age_years^0.8493 ...
    - 5.723e-7*lambda2^3 + 0.001245*lambda2^2 - 0.9025*lambda2;

deltaC_oxy = (((epsilon(1)*deltaOD_infrared)/DPF_lambda2) ...
    - ((epsilon(3)*deltaOD_red)/DPF_lambda1)) ...
    / (((epsilon(1)*epsilon(4))-(epsilon(3)*epsilon(2))) ...
    * source_detector_distance_cm);

deltaC_deoxy = (((epsilon(4)*deltaOD_red)/DPF_lambda1) ...
    - ((epsilon(2)*deltaOD_infrared)/DPF_lambda2)) ...
    / (((epsilon(1)*epsilon(4))-(epsilon(3)*epsilon(2))) ...
    * source_detector_distance_cm);

if gr
    figure
    plot(time,deltaC_oxy,'r','LineWidth',1.2)
    hold on
    plot(time,deltaC_deoxy,'b','LineWidth',1.2)
    xline(time(i_phases(1,:)),'--')
    xline(time(i_phases(2,:)),'--')
    hold off
    legend('\DeltaHbO_2','\DeltaHHb')
    ylabel('Concentration change [mM]')
    xlabel('Time [s]')
    title('fNIRS concentration change')
end
end
