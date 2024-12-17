fprintf('Analyzing audio2 for overtone distances...\n');
audio2 = AudioProcessor('gtr-jazz.wav');
audio2 = audio2.segmentAudio(1024, 512, 0.2, 2);

[longestSegment2, ~, ~] = audio2.findLongestSegment();

if ~isempty(longestSegment2)
    % Analyze audio2: Get overtone distances and amplitudes
    thresholdFactor = 0.05;
    [~, overtoneDistances2, sortedPks2] = audio2.analyzeSegment(longestSegment2, thresholdFactor);
    fprintf('Overtone distances (Hz) from audio2: %s\n', num2str(overtoneDistances2));
else
    error('No valid segment found in audio2 for analysis.');
end

fprintf('\nProcessing audio1 and generating combined notes...\n');
audio1 = AudioProcessor('pno-cs.wav');
audio1 = audio1.segmentAudio(1024, 512, 0.2, 3);
audio1 = audio1.removeHarmonics();

estimatedTotalLength = sum(cellfun(@length, audio1.FilteredSegments)) + 1000; 
combinedSignal = zeros(estimatedTotalLength, 1); % Preallocate combined signal
currentPosition = 1; % Track the current write position

for i = 1:length(audio1.FilteredSegments)
    fprintf('Processing note %d...\n', i);
    
    currentSegment = audio1.FilteredSegments{i};
    
    [fundamentalFreq1, ~, ~] = audio1.analyzeSegment(currentSegment, 0.05);
    fprintf('Fundamental frequency of note %d: %.2f Hz\n', i, fundamentalFreq1);

    [~, overtoneDistances2, sortedPks2, sortedPhases2] = audio2.analyzeSegment(longestSegment2, 0.05);

    numOvertones = 10; % Specify how many overtones to use
    generatedNote = audio2.generateSimulatedSignal(...
        fundamentalFreq1, overtoneDistances2, sortedPks2, sortedPhases2, numOvertones, currentSegment);

    generatedNote = generatedNote(:);

    noteLength = length(generatedNote);
    combinedSignal(currentPosition:currentPosition + noteLength - 1) = generatedNote;

    currentPosition = currentPosition + noteLength;
end

combinedSignal = combinedSignal(1:currentPosition - 1);

combinedSignal = combinedSignal / max(abs(combinedSignal));

outputFileName = 'combined_generated_notes.wav';
audiowrite(outputFileName, combinedSignal, audio1.Fs);

fprintf('Combined notes with phase have been exported to %s\n', outputFileName);
