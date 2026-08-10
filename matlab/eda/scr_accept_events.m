function keep = scr_accept_events(amplitudes,minimumScrAmplitude)
%SCR_ACCEPT_EVENTS Apply the inclusive absolute SCR amplitude threshold.

keep = amplitudes >= minimumScrAmplitude;
end
