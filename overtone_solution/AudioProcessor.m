classdef AudioProcessor
    properties
        Audio            % Original audio signal
        Fs               % Sampling frequency
        Segments         % Segments of the audio
        FilteredSegments % Harmonic-filtered segments
    end
    
    methods
        % Constructor: Load audio file
        function obj = AudioProcessor(audioFile)
            [audio, fs] = audioread(audioFile);
            if size(audio, 2) > 1
                audio = mean(audio, 2); % Convert stereo to mono
            end
            obj.Audio = audio;
            obj.Fs = fs;
        end

        % Perform segmentation
        function obj = segmentAudio(obj, windowSize, hopSize, thresholdFactor, deleteLastN)
            spectrum = abs(spectrogram(obj.Audio, windowSize, windowSize - hopSize));
            energy = sum(spectrum, 1);
            diffEnergy = [0, diff(energy)];
            threshold = thresholdFactor * max(diffEnergy);
            onsets = find(diffEnergy > threshold);
            onsetSamples = onsets * hopSize;

            % Divide audio into segments
            obj.Segments = {};
            for i = 1:length(onsetSamples)
                if i < length(onsetSamples)
                    obj.Segments{i} = obj.Audio(onsetSamples(i):onsetSamples(i+1)-1);
                else
                    obj.Segments{i} = obj.Audio(onsetSamples(i):end);
                end
            end

            % Optionally delete the last N segments
            if deleteLastN > 0
                obj = obj.removeLastSegments(deleteLastN);
            end
            
            fprintf('Audio segmented into %d parts.\n', length(obj.Segments));
        end

        % Find the longest segment
        function [longestSegment, longestIndex, longestDuration] = findLongestSegment(obj)
            longestDuration = 0;
            longestIndex = 0;
            longestSegment = [];
            for i = 1:length(obj.Segments)
                segmentDuration = length(obj.Segments{i}) / obj.Fs;
                if segmentDuration > longestDuration
                    longestDuration = segmentDuration;
                    longestIndex = i;
                    longestSegment = obj.Segments{i};
                end
            end
            fprintf('Longest segment index: %d, Duration: %.2f seconds.\n', longestIndex, longestDuration);
        end

        % Remove the last N segments
        function obj = removeLastSegments(obj, n)
            if length(obj.Segments) >= n
                obj.Segments(end-n+1:end) = [];
                fprintf('Last %d segments removed. Remaining: %d segments.\n', n, length(obj.Segments));
            else
                fprintf('Not enough segments to remove. Remaining unchanged.\n');
            end
        end

        % Remove harmonics from each segment
        function obj = removeHarmonics(obj)
            obj.FilteredSegments = cell(size(obj.Segments));
            for i = 1:length(obj.Segments)
                segment = obj.Segments{i};
                N = length(segment);
                fftSegment = fft(segment);
                magnitude = abs(fftSegment);
                [~, fundamentalIdx] = max(magnitude(1:round(N/2)));

                % Keep only the fundamental frequency
                filteredFFT = zeros(size(fftSegment));
                filteredFFT(fundamentalIdx) = fftSegment(fundamentalIdx);
                filteredFFT(N-fundamentalIdx+1) = fftSegment(N-fundamentalIdx+1);

                % Inverse FFT to reconstruct the signal
                obj.FilteredSegments{i} = real(ifft(filteredFFT));
            end
            fprintf('Harmonics removed from all segments.\n');
        end

        % Analyze a segment: Detect fundamental frequency and filter overtones
        function [fundamentalFreq, overtoneDistances, sortedPks, sortedPhases] = analyzeSegment(obj, segment, thresholdFactor)
            % Perform FFT
            N = length(segment);
            frequencies = (0:N-1) * (obj.Fs / N);
            fftSegment = fft(segment);
            magnitude = abs(fftSegment(1:floor(N/2)));
            phase = angle(fftSegment(1:floor(N/2))); % Extract phase information
            frequencies = frequencies(1:floor(N/2));
        
            % Detect peaks in the magnitude spectrum
            [pks, locs] = findpeaks(magnitude, frequencies, 'MinPeakHeight', max(magnitude) * thresholdFactor);
        
            % Find the fundamental frequency (highest amplitude peak)
            [~, maxIdx] = max(pks);
            fundamentalFreq = locs(maxIdx);
        
            % Sort frequencies, amplitudes, and phases by peak magnitude
            [sortedPks, sortIdx] = sort(pks, 'descend');
            sortedFreqs = locs(sortIdx);
            sortedPhases = phase(sortIdx); % Sort phases corresponding to peaks
        
            % Calculate overtone distances relative to the fundamental frequency
            overtoneDistances = sortedFreqs(sortedFreqs ~= fundamentalFreq) - fundamentalFreq;
        
            % Remove the fundamental frequency phase from overtones
            sortedPhases = sortedPhases(sortedFreqs ~= fundamentalFreq);
        end



        % Generate a simulated sound using magnitude and phase
        function simulatedSignal = generateSimulatedSignal(obj, fundamentalFreq, overtoneDistances, sortedPks, sortedPhases, numOvertones, segment)
            N = length(segment);
            t = (0:N-1) / obj.Fs;
        
            % Start with the fundamental tone (phase = 0)
            simulatedSignal = sin(2 * pi * fundamentalFreq * t);
        
            % Add specified number of overtones with phase
            numOvertones = min(numOvertones, length(overtoneDistances)); % Limit to available overtones
            for i = 1:numOvertones
                overtoneFreq = fundamentalFreq + overtoneDistances(i);
                relativeAmp = sortedPks(i) / max(sortedPks); % Normalize amplitude
                phase = sortedPhases(i); % Retrieve phase information
        
                % Reconstruct overtone with magnitude and phase
                simulatedSignal = simulatedSignal + relativeAmp * sin(2 * pi * overtoneFreq * t + phase);
            end
        
            % Scale the simulated signal to match the peak amplitude of the original segment
            peakAmplitudeOriginal = max(abs(segment));
            simulatedSignal = simulatedSignal * (peakAmplitudeOriginal / max(abs(simulatedSignal)));
        end

        % Play all filtered segments
        function testFilteredSegments(obj)
            for i = 1:length(obj.FilteredSegments)
                fprintf('Playing filtered segment %d...\n', i);
                sound(obj.FilteredSegments{i}, obj.Fs);
                pause(length(obj.FilteredSegments{i}) / obj.Fs + 1);
            end
            fprintf('All filtered segments played.\n');
        end
    end
end
