function filtered = fnirs_zero_phase_sos(sos,signal,paddingLength)
%FNIRS_ZERO_PHASE_SOS Apply SOS filtering with explicit odd reflection.

signal = signal(:);
if numel(signal) < 2 || any(~isfinite(signal))
    error('MultimodalPhysioKit:FNIRS:InvalidZeroPhaseSignal', ...
        'Zero-phase filtering requires at least two finite samples.')
end
if ~isscalar(paddingLength) || ~isfinite(paddingLength) || ...
        paddingLength < 1 || paddingLength ~= fix(paddingLength)
    error('MultimodalPhysioKit:FNIRS:InvalidPaddingLength', ...
        'paddingLength must be a positive integer.')
end

paddingLength = min(paddingLength,numel(signal)-1);
prePadding = 2*signal(1) - flipud(signal(2:paddingLength+1));
postPadding = 2*signal(end) - ...
    flipud(signal(end-paddingLength:end-1));
padded = [prePadding; signal; postPadding];

forward = localSosCascade(sos,padded);
backward = localSosCascade(sos,flipud(forward));
filteredPadded = flipud(backward);
filtered = filteredPadded( ...
    paddingLength+1:paddingLength+numel(signal));
end

function output = localSosCascade(sos,input)
output = input;
for sectionIndex = 1:size(sos,1)
    b = sos(sectionIndex,1:3);
    a = sos(sectionIndex,4:6);
    sectionInput = output(1);
    sectionOutput = sum(b)/sum(a)*sectionInput;
    initialConditions = [ ...
        sectionOutput-b(1)*sectionInput; ...
        b(3)*sectionInput-a(3)*sectionOutput];
    output = filter(b,a,output,initialConditions);
end
end
