% segmentByClicking
% This code provides a semi-automated segmentation as an alternative to the combination of the functions identifyPeaks_dlg/setThreshold_dlg/setSyllableBoundaries.
% You need to click inside of the amplitude time course plot where you think a good segmentation threshold would sit - then an interpolated line will be
% drawn between the click positions. This will serve as threshold. Finally you either confirm or redo the segmentation. 
% Input variable is an amplitude time course, for instance amp_timeseries. 
% Output variables will be pklocs_AllPeaks, pkheights_AllPeaks,thresh (a vector as long as amplitude_timeseries),SyllStarts,SyllEnds.



%% plotting figure with amplitude time series

%     figure; 
%     plot(amp_timeseries);
%     hold on
%     plot(hpfilter(amp_timeseries,1000),'g')
%     title('Please resize plot now')

%     pause(12) % time to resize the plot
%     hold on;

%     title('Please click into plot for threshold positions, then press Enter')

%% recording positions of up to 150 clicks into the plot

    [x,y] = ginput(150);

%% making a threshold by interpolating the dots.

try    
    yy = interp1(x,y,[1:1:numel(amp_timeseries)]');
catch
    warning('Problem using interpolation between the points you clicked. You have one more shot...');
    [x,y] = ginput(150);
    yy = interp1(x,y,[1:1:numel(amp_timeseries)]');
end
    thresh = yy;

    hold on; plot(thresh,'m')

%% double checking: are you happy with the threshold?

    a = 1;
    while a,
        
        prompt = {'Accept this threshold? [Y/N]'};
                
            dlg_title = 'Input';
            num_lines = 1;
            default = {''};
            w = inputdlg(prompt,dlg_title,num_lines,default);

            if strmatch(w,'Y','exact') % Threshold can be used for segmenting.
                a = 0;
            else % Threshold not accepted, redo.
                [x,y] = ginput(150);
                yy = interp1(x,y,[1:1:numel(amp_timeseries)]');
                thresh = yy;
                hold on; plot(thresh)
            end
            
    end
                
  %% extracting SyllStarts, SyllEnds, pklocs_AllPeaks, and pkheights_AllPeaks               
                        
 ampcut_timeseries = amp_timeseries-thresh; 
 sound_ix = find(ampcut_timeseries>0);
 silent_ix = find(ampcut_timeseries<=0);
 
 diff_sound_ix = diff(sound_ix)-1;
%  figure; plot(diff_sound_ix);
 SyllStarts = [sound_ix(1);sound_ix(find(diff_sound_ix>0)+1)];
 SyllEnds = [sound_ix(find(diff_sound_ix>0));sound_ix(end)];
 for v=1:numel(SyllStarts)
        hold on; plot(SyllStarts(v):SyllEnds(v),10*ones(SyllEnds(v)-SyllStarts(v)+1,1),'.r')
 end
 
 %% extracting peaklocs and peakheights
 
 pklocs_AllPeaks = zeros(numel(SyllStarts),1);
 pkheights_AllPeaks = zeros(numel(SyllStarts),1);
 [prelimpeakheights,prelimpeaklocs] = findpeaks(amp_timeseries,'MinPeakHeight',min(thresh)); 
%  allpks = findpeaks(amp_timeseries,thresh); 
%             prelimpeaklocs = allpks.loc;
%             prelimpeakheights = amp_timeseries(prelimpeaklocs);
            hold on; plot(prelimpeaklocs,prelimpeakheights,'xb')
 for v=1:numel(SyllStarts)
     thisSyll_ix = [SyllStarts(v):SyllEnds(v)]';
     ovrlpPks_ix = intersect(prelimpeaklocs,thisSyll_ix);
     if numel(ovrlpPks_ix)==1
         pklocs_AllPeaks(v) = ovrlpPks_ix;
         pkheights_AllPeaks(v) = amp_timeseries(pklocs_AllPeaks(v));
     elseif numel(ovrlpPks_ix)>1 % if more than 1 peak, use the first one that exceeds half the distance between max peak and threshold.
         maxdist2thresh = max(ampcut_timeseries(ovrlpPks_ix));
         dist_candidates = ampcut_timeseries(ovrlpPks_ix);
         ix_ofOvrlp = find(dist_candidates>(maxdist2thresh/2),1,'first');
         pklocs_AllPeaks(v) = ovrlpPks_ix(ix_ofOvrlp);
         pkheights_AllPeaks(v) = amp_timeseries(pklocs_AllPeaks(v));
     else
     end
 end
 
       hold on; plot(pklocs_AllPeaks,pkheights_AllPeaks,'or')
       txt_h = labelpoints(pklocs_AllPeaks, pkheights_AllPeaks, [1:1:(numel(pklocs_AllPeaks))], 'N', 0.03, 1, 'FontSize', 8, 'Color', 'r');
 
                  