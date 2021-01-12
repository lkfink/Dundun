% Quick function to grab means for each stimulus for a feature of interest
% LF 20201125

function means = getMeansForFeat(feat_mus, feat_speech)

% Input:
% - vectors of features for music and speech, respectively (custom to
% how Tina has data organized)
% Output:
% - means for each stimulus (ordered 1:15 M, 16:30 S)

mus = splitOnNaN(feat_mus);
mus_means = cellfun(@mean, mus);
spe = splitOnNaN(feat_speech);
spe_means = cellfun(@mean, spe);
means = vertcat(mus_means, spe_means);

end