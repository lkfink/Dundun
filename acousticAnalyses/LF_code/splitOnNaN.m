% Quick function to split up Tina's feature vectors by stimulus
% (stimuli were separated by NaNs)
% LF - 20201125

function C = splitOnNaN(M)
    % Input: 
    % - double containing values for each stimulus, separated by NaNs

    % Output
    % - cell array with values for each stimulus separately
    
    index=find(~isnan(M));
    idx=find(diff(index)~=1);
    A=[idx(1);diff(idx);numel(index)-idx(end)];
    C=mat2cell(M(~isnan(M)),A,1);
    
end