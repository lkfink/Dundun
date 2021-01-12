% Quick function to grab peak in AMS spectrum for further analyses
% LF 20201126

function peakFreqs = getAmsPeak(mus_spectra, spe_spectra, allfreqs)

% Input:
% - amplitude modulation spectrum for each music stimulus (15,600 double)
% - amplitude modulation spectrum for each speech stimulus (15,600 double)
% - frequency values corresponding to each amplitude value (1,600 double)

% Output: 
% - Frequency at which highest AMS value is observed

% Get peaks for music stimuli
[~, idx] = max(mus_spectra, [], 2); % give us index of maximum value in amplitude array (each row)
peakFreqs_mus = allfreqs(idx); % gives us frequency corresponding to highest amplitude

% Get peaks for speech stimuli
[~, idx] = max(spe_spectra, [], 2); % give us index of maximum value in amplitude array
peakFreqs_spe = allfreqs(idx); % gives us frequency corresponding to highest amplitude

% Combine music and speech peak freqs into one vector
peakFreqs = vertcat(peakFreqs_mus, peakFreqs_spe);

end

