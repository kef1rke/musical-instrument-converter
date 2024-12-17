frequencyThreshold = 5;

fprintf('Segmenting audio2...\n');
audio2 = AudioProcessor('gtr-jazz.wav');
audio2 = audio2.segmentAudio(1024, 512, 0.2, 3);

% Analyze fundamental frequencies of audio2 segments
audio2Fundamentals = zeros(1, length(audio2.Segments));
for j = 1:length(audio2.Segments)
    [audio2Fundamentals(j), ~, ~] = audio2.analyzeSegment(audio2.Segments{j}, 0.05);
end

fprintf('Segmenting and removing harmonics from audio1...\n');
audio1 = AudioProcessor('pno-cs.wav');
audio1 = audio1.segmentAudio(1024, 512, 0.2, 3);
audio1 = audio1.removeHarmonics(); % Remove harmonics

audio1Fundamentals = zeros(1, length(audio1.FilteredSegments));
for i = 1:length(audio1.FilteredSegments)
    [audio1Fundamentals(i), ~, ~] = audio1.analyzeSegment(audio1.FilteredSegments{i}, 0.05);
    fprintf('Audio1 segment %d fundamental frequency: %.2f Hz\n', i, audio1Fundamentals(i));
end

fprintf('\nMatching audio1 segments to audio2 samples...\n');
% Preallocate a large buffer for the combined signal
estimatedLength = sum(cellfun(@length, audio1.FilteredSegments)) + 1000; % Safety buffer
combinedSignal = zeros(estimatedLength, 1);
currentPosition = 1;
fs = audio1.Fs;

for i = 1:length(audio1.FilteredSegments)
    fprintf('Processing audio1 segment %d...\n', i);
    
    fundamentalFreq1 = audio1Fundamentals(i);

    [minDiff, closestIndex] = min(abs(audio2Fundamentals - fundamentalFreq1));

    if minDiff <= frequencyThreshold
        fprintf('Using audio2 segment %d (difference: %.2f Hz)\n', closestIndex, minDiff);
        matchingSample = audio2.Segments{closestIndex};
    else
        fprintf('Pitch shifting audio2 segment %d to match frequency %.2f Hz...\n', closestIndex, fundamentalFreq1);
        matchingSample = audio2.Segments{closestIndex};
        originalFreq = audio2Fundamentals(closestIndex);
        matchingSample = pitchShift(matchingSample, originalFreq, fundamentalFreq1, fs);
    end

    matchingSample = adjustLength(matchingSample, length(audio1.FilteredSegments{i}));

    matchingSample = matchingSample(:);

    noteLength = length(matchingSample);
    combinedSignal(currentPosition:currentPosition + noteLength - 1) = matchingSample;
    currentPosition = currentPosition + noteLength;
end

combinedSignal = combinedSignal(1:currentPosition - 1);

combinedSignal = combinedSignal / max(abs(combinedSignal));
outputFileName = 'converted_audio.wav';
audiowrite(outputFileName, combinedSignal, fs);

fprintf('Converted audio has been exported to %s\n', outputFileName);

% ----------------- Helper Functions ----------------- %

function shiftedSignal = pitchShift(signal, originalFreq, targetFreq, fs)
    % Validate the frequencies
    if originalFreq <= 0 || targetFreq <= 0
        error('Invalid frequency values: originalFreq and targetFreq must be positive and non-zero.');
    end
    
    resampleFactor = targetFreq / originalFreq;

    if resampleFactor > 2
        fprintf('Resampling factor %.2f is too large. Limiting to 2.\n', resampleFactor);
        resampleFactor = 2;
    elseif resampleFactor < 0.5
        fprintf('Resampling factor %.2f is too small. Limiting to 0.5.\n', resampleFactor);
        resampleFactor = 0.5;
    end

    N = length(signal);
    tOriginal = (0:N-1) / fs;
    tNew = linspace(0, tOriginal(end), round(N * resampleFactor)); % New time vector

    shiftedSignal = interp1(tOriginal, signal, tNew, 'linear', 0); % Linear interpolation, zero-padding
end

function adjustedSignal = adjustLength(signal, targetLength)
    signal = signal(:);
    
    currentLength = length(signal);
    if currentLength > targetLength
        adjustedSignal = signal(1:targetLength);
    elseif currentLength < targetLength
        adjustedSignal = [signal; zeros(targetLength - currentLength, 1)];
    else
        adjustedSignal = signal;
    end
end
