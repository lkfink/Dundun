% analyze_dundun.m

clear all; close all; clc

%% loading the corpus data that contains raw amplitude, frequency, entropy for all phrases.

cache_dir='~/processed'; 
data_dir='~/Dropbox/MATLAB/Dundun/data/wavs';
results_dir='~/Dropbox/MATLAB/Dundun/results';
savefile = 'Dundun_Rhythm.mat';

MAX_READ_FILE=1000;
is_create=false; % set to true if you need the wav files read in anew. Otherwise, the code just loads the existent corpus data.

corpus=load_dundun_corpus(cache_dir,data_dir,MAX_READ_FILE,is_create);

%% segmentation

IS_FILTERED_FEATURES = true; % would you like to use raw entropy and frequency or a filtered version?
IS_SEGMENTBYCLICK = true; % if you want to segment song in a semi-automated way (clicking into plot), set this to true.
IS_SCALED = true; 
% otherwise, segmentation will be based on filtered time series...
IS_CREATENEWRHYTHM = false; % set this to true if you want to re-do the segmentation.

skipped_wavs = cell(1,1);
sk = 1;

TESTRUNNUM = 0; % enter any number if you want to do test runs to compare

%% cycling through corpus

if IS_CREATENEWRHYTHM

        for i = 1:numel(corpus) 

            thisWavName = char(corpus{i}.wavname);  fprintf('Segmenting file %s...\n',thisWavName)      
            ampraw = corpus{i}.rawamplitude;        ampraw(ampraw==-70) = 0;    
            entraw = corpus{i}.rawentropy;          entraw(entraw>=2) = 0;    
            freqraw = corpus{i}.rawfrequency;       %freqraw()    
            
%             % read in yin
%             rdir=pwd;
%             cd(data_dir)
%             YINTEST = yin(thisWavName)
%             cd(rdir)
            if IS_SCALED
                ampscaled = (ampraw-(prctile(ampraw,2))) / (prctile(ampraw,98)-(prctile(ampraw,2))); % scaled versions are basically only needed for plotting
                entscaled = (entraw-(prctile(entraw,2))) / (prctile(entraw,98)-(prctile(entraw,2))); % scaled versions are basically only needed for plotting
                freqscaled = (freqraw-(prctile(freqraw,2))) / (prctile(freqraw,98)-(prctile(freqraw,2))); % scaled versions are basically only needed for plotting
    %             yinscaled = (yinraw-(prctile(yinraw,2))) / (prctile(yinraw,98)-(prctile(yinraw,2))); % scaled versions are basically only needed for plotting
            else
                ampscaled = ampraw;
                entscaled = entraw;
                freqscaled = freqraw;
            end


            if IS_FILTERED_FEATURES
                ampscaled = hpfilter(ampscaled,50);
                entscaled = hpfilter(entscaled,50);
                freqscaled = hpfilter(freqscaled,50);
%                 yinscaled = hpfilter(yinscaled,50);
            end
            corpus{i}.scaledamplitude = ampscaled;
            corpus{i}.scaledfrequency = freqscaled;
            corpus{i}.scaledentropy = entscaled;
            

            %%%%%%%%%%%%%%%%                
            % SEGMENTATION %
            %%%%%%%%%%%%%%%%

            figure(55); clf;
            set(gcf,'Position',[1 528 1680 427]);
            hold on
%             plot(ampscaled,'b','linewidth',1.2);    plot(entscaled,'Color',[1 .8 0]);   plot(freqscaled,'r');       plot(yinscaled,'c')
            axis([0 numel(ampscaled) -.12 1.3])
            title(thisWavName)

            % playing wav 
            cd(data_dir)
            [audio4sound,fs_original] = audioread(thisWavName);
            sound(audio4sound,fs_original)
            pause(length(audio4sound)/fs_original+1)
            fs_slow = fs_original/3; % blackcaps: 5.5, thrush nightingales and blackbirds: 3
            sound(audio4sound,fs_slow)
            pause(length(audio4sound)/fs_slow+1)


                        % here you get the option to SKIP THIS TRILL:
                        prompt = {'Do you want to use this song? [Y/N]'};

                            dlg_title = 'Input';
                            num_lines = 1;
                            default = {''};
                            w = inputdlg(prompt,dlg_title,num_lines,default);

                            if strmatch(w,'Y','exact') % This trill will be used. 


                                if IS_SEGMENTBYCLICK

                                            % 1st way to segment: as an alternative, perform the segmentation by clicking into the plot where you want your threshold to be.
                                            %   this is done by the script "segmentByClicking".

                                                amp_timeseries = ampscaled;
                                                segmentByClicking; % output here is pklocs_AllPeaks, pkheights_AllPeaks, thresh, SyllStarts, SyllEnds.


                                else

                                            % 2nd way to segment: just the way it's done in the nightingale_trill_toussaint code:

                                                [pklocs_AllPeaks,pkheights_AllPeaks] = identifyPeaks_dlg(allpeaklocs,allpeakheights);
                                                [thresh] = setThreshold_dlg(exampleThresholds); %,ampfilt_trill);
                                                [SyllStarts,SyllEnds] = setSyllableBoundaries(m_amplitude_trill,ampfilt_trill,pklocs_AllPeaks,pkheights_AllPeaks,thresh);

                                end

                %% creating Rhythm structure containing all the data


                                        % CALCULATE AND SAVE NOTE LENGTH 
                                        %%%%%%%%%%% WEITERMACHEN: enter these results into a variable that gets saved...
                                        Rhythm(i).wavname = thisWavName;
                                        Rhythm(i).rawamplitude = ampraw;
                                        Rhythm(i).rawfrequency = freqraw;
                                        Rhythm(i).rawentropy = entraw;
                                        Rhythm(i).filteredamp = ampscaled;
                                        Rhythm(i).peaklocs_allpeaks = pklocs_AllPeaks;
                                        Rhythm(i).peakheights_allpeaks = pkheights_AllPeaks;
                                        Rhythm(i).syllablestarts = SyllStarts;
                                        Rhythm(i).syllableends = SyllEnds;
                                        Rhythm(i).sylldur = SyllEnds-SyllStarts;
                                        Rhythm(i).pausedur = [SyllStarts(2:end)-SyllEnds(1:end-1);nan];

                                        meansyllamp=nan(size(pklocs_AllPeaks));  meansyllent=nan(size(pklocs_AllPeaks));   meansyllfreq=nan(size(pklocs_AllPeaks));% mean syllable features

                                        for j = 1:numel(pklocs_AllPeaks)
                                            meansyllamp(j) = mean(ampraw(SyllStarts(j):SyllEnds(j)));
                                            meansyllent(j) = mean(entraw(SyllStarts(j):SyllEnds(j)));
                                            meansyllfreq(j) = mean(freqraw(SyllStarts(j):SyllEnds(j)));
                                        end
                                        Rhythm(i).mean_amplitude_syll = meansyllamp;
                                        Rhythm(i).mean_entropy_syll = meansyllent;
                                        Rhythm(i).mean_frequency_syll = meansyllfreq;


                            else % choosing to skip this trill. THIS DOESN'T WORK YET!!
                                disp(['You chose to skip trill ',thisWavName]);
                                skipped_wavs{sk,1} = thisWavName;
                                skipped_wavs{sk,2} = i;
                                sk = sk+1;
                            end

        end
        
       %% saving Rhythm - CAREFUL! Don't overwrite the current version

        % cd(results_dir)
        % if ~exist('savefile')
        %     save(savefile,'Rhythm')
        % else
        %     save(savefile,'Rhythm','-append')
        % end

else
    cd(results_dir)
    load('Dundun_Rhythm_0920.mat')
end

if ~isfield(Rhythm,'rawpitch')
    [Rhythm] = addPitch2RhythmStruct(Rhythm,data_dir);
end


%% %%------------------------------------------------------------------------------------------------------------------------------------------%%
%%%% ------------------------------------------------------- Visualizing Dundun rhythm --------------------------------------------------------%%
%%%%%------------------------------------------------------------------------------------------------------------------------------------------%%


%% (1) RASTER PLOT: image from scaled amplitude in false colors, aligned to first peak.

    Dundun_M = Rhythm(1:2:end);
    Dundun_s = Rhythm(2:2:end);

    alpoint = 2000;
    plotwidth = 13000;
    raster = nan(numel(Dundun_M)+numel(Dundun_s)+1,plotwidth);

    % aligning to first note peak 
    for i=1:numel(Dundun_M)
        ampfilt = [Dundun_M(i).filteredamp];
            ampfilt(ampfilt<-.1) = -.1;
        pklocs = [Dundun_M(i).peaklocs_allpeaks];
        ampleft = ampfilt(1:pklocs(1));
        ampright = ampfilt(pklocs(1)+1:end);
        raster(i,alpoint-numel(ampleft):alpoint-1) = ampleft;
        raster(i,alpoint:alpoint+numel(ampright)-1) = ampright;
    end
    for i=1:numel(Dundun_s)
        ampfilt = [Dundun_s(i).filteredamp];
            ampfilt(ampfilt<-.1) = -.1;
        pklocs = [Dundun_s(i).peaklocs_allpeaks];
        ampleft = ampfilt(1:pklocs(1));
        ampright = ampfilt(pklocs(1)+1:end);
        raster(i+numel(Dundun_M)+1,alpoint-numel(ampleft):alpoint-1) = ampleft;
        raster(i+numel(Dundun_M)+1,alpoint:alpoint+numel(ampright)-1) = ampright;
    end

    figure(1+100*TESTRUNNUM); clf;
    imagesc(raster);  colormap(jet);     title('Intensity, Dundun_M top, Dundun_s bottom','Interpreter','none')

%% 1b) Dundun autocorrelation #WEITERMACHEN

    figure(100); clf; title('TEST music amp autocorrelation')
    figure(101); clf; title('TEST speech amp autocorrelation')
    for i=1:numel(Dundun_M)
        ampfilt = [Dundun_M(i).filteredamp];
        [r,lags] = xcorr(ampfilt,'coeff');
        figure(100); hold on
        plot(lags,r)
    end
    for i=1:numel(Dundun_s)
        ampfilt = [Dundun_s(i).filteredamp];
            ampfilt(ampfilt<-.1) = -.1;
        [r,lags] = xcorr(ampfilt,'coeff');
        figure(101); hold on
        plot(lags,r)
    end
    figure(100); xlabel('lags'); ylabel('r');  box on; ylim([0 1]); xlim([-10000 10000]);
    figure(101); xlabel('lags'); ylabel('r');  box on; ylim([0 1]); xlim([-10000 10000]);


    
%% (2) note rate scatter (music left, speech right)
    
    noterate_M_then_s = nan(numel(Dundun_M)*2,1);     % for creating table for Pauline later
    filenames_M_then_s = cell(numel(Dundun_M)*2,1);     % for creating table for Pauline later
    
    figure(2+100*TESTRUNNUM); set(gcf,'Position',[594 583 210 282]); clf;
    for i = 1:numel(Dundun_M)
        filenames_M_then_s{i} = Dundun_M(i).wavname;
        ampfilt = [Dundun_M(i).filteredamp];
        pklocs = [Dundun_M(i).peaklocs_allpeaks];
        dur = pklocs(end)-pklocs(1);
        noterate_M_then_s(i) = (numel(pklocs)-1)/(dur/1000);
        hold on;
        plot(1+randn(1,1)/15,noterate_M_then_s(i),'ob')
        
    end
    for i = 1:numel(Dundun_s)
        filenames_M_then_s{i+numel(Dundun_M)} = Dundun_s(i).wavname;
        ampfilt = [Dundun_s(i).filteredamp];
        pklocs = [Dundun_s(i).peaklocs_allpeaks];
        dur = pklocs(end)-pklocs(1);
        noterate_M_then_s(i+numel(Dundun_M)) = (numel(pklocs)-1)/(dur/1000);
        hold on;
        plot(2+randn(1,1)/15,noterate_M_then_s(i+numel(Dundun_M)),'or')
    end
    ylim([0 9]);    xlim([0 3]);    xticks([1 2]);  xticklabels({'music-like','speech-like'});  box on
    ylabel('Note rate [Hz]');   title('Dundun note rates')
    
    
    
%% (3) FEATURES (mean intensity and mean frequency) of notes plotted on time-x-axis. 
    
NoteFeatures_M = cell(numel(Dundun_M),1);
NoteFeatures_s = cell(numel(Dundun_s),1);
  
    % Music-like
    tic
    for i = 1:numel(Dundun_M)
        
        IS_NOTELENGTH30 = true; % if you want to use a standard note length of 30 ms for mean feature calculation, from syllable onset on
        USE_NOTE_PEAK = false; % if you want to use just 19ms around the peak for note intensity calculation
        lineInRhythmStruct = i;
        IS_FILTERED_FEATURES = false;
        IS_PLOT=false;

        [NoteFeaturesSong] = calcNoteFeaturesFromDundunRhythm(Dundun_M,i,IS_NOTELENGTH30,USE_NOTE_PEAK,IS_FILTERED_FEATURES,IS_PLOT); 

        NoteFeatures_M{i,1} = NoteFeaturesSong; % there are scaled and raw mean note features in here.
        
    end
    
    % speech-like
    for i = 1:numel(Dundun_s)
       
        IS_NOTELENGTH30 = true; % if you want to use a standard note length of 30 ms for mean feature calculation, from syllable onset on
        USE_NOTE_PEAK = false; % if you want to use just 19ms around the peak for note intensity calculation
        IS_FILTERED_FEATURES = false;
        IS_PLOT=false; 
        
        [NoteFeaturesSong] = calcNoteFeaturesFromDundunRhythm(Dundun_s,i,IS_NOTELENGTH30,USE_NOTE_PEAK,IS_FILTERED_FEATURES,IS_PLOT); 

        NoteFeatures_s{i,1} = NoteFeaturesSong; % there are scaled and raw mean note features in here.
        
        
    end
    toc
    
    IS_LONG_PLOT = false;
    
    if IS_LONG_PLOT
        
            % now plotting
            figure(11+100*TESTRUNNUM); set(gcf,'Position',[1 111 1664 823]); clf % for amplitude
                [ha, posa] = tight_subplot(15,2,.003,[.1],[.1]); % [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)   
            figure(12+100*TESTRUNNUM); set(gcf,'Position',[1 111 1664 823]); clf % for entropy
                [hp, posp] = tight_subplot(15,2,.003,[.1],[.1]); % [he, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)   

            % music - plotting features of all individual pieces, continuously and for notes. 
            for k = 1:numel(NoteFeatures_M)

                N_ampvect = [NoteFeatures_M{k,1}.AmpScaled_vect];    N_freqvect = [NoteFeatures_M{k,1}.FreqScaled_vect];    
                N_entvect = [NoteFeatures_M{k,1}.EntScaled_vect];    N_ptcvect = [NoteFeatures_M{k,1}.PitchScaled_vect];
                N_noteix = [NoteFeatures_M{k,1}.Note_ix_vect];

                figure(11+100*TESTRUNNUM);   subplot(ha(k*2-1));  hold on;  

                    amp = [Dundun_M(k).rawamplitude];     ampsc = (amp-(prctile(amp,2))) / (prctile(amp,98)-(prctile(amp,2)));          ampsc = hpfilter(ampsc,50);

                    pl=plot(ampsc,'b','linewidth',.05);     pl.Color(4) = 0.25;   box on
                    plot(N_noteix,N_ampvect,'.b');
                    axis([0 12000 0 1.1]); 

               figure(12+100*TESTRUNNUM);   subplot(hp(k*2-1));  hold on;  

                    ptc = [Dundun_M(k).rawpitch];     ptcsc = (ptc-(prctile(ptc,2))) / (prctile(ptc,98)-(prctile(ptc,2)));          ptcsc = hpfilter(ptcsc,10);

                    pl=plot(ptcsc,'b','linewidth',.05);     pl.Color(4) = 0.25;   box on
                    plot(N_noteix,N_ptcvect,'.b');              
                    axis([0 12000 0 1.1]); 

            end

            % speech - plotting features of all individual pieces, continuously and for notes. 
            for k = 1:numel(NoteFeatures_s)

                N_ampvect = [NoteFeatures_s{k,1}.AmpScaled_vect];    N_freqvect = [NoteFeatures_s{k,1}.FreqScaled_vect];    
                N_entvect = [NoteFeatures_s{k,1}.EntScaled_vect];    N_ptcvect = [NoteFeatures_s{k,1}.PitchScaled_vect];
                N_noteix = [NoteFeatures_s{k,1}.Note_ix_vect];


                figure(11+100*TESTRUNNUM);   subplot(ha(k*2));  hold on;   % amplitude

                    amp = [Dundun_s(k).rawamplitude];     ampsc = (amp-(prctile(amp,2))) / (prctile(amp,98)-(prctile(amp,2)));          ampsc = hpfilter(ampsc,50);

                    pl=plot(ampsc,'r','linewidth',.05);     pl.Color(4) = 0.25;   box on
                    plot(N_noteix,N_ampvect,'.r');
                    axis([0 12000 0 1.1]); 

                figure(12+100*TESTRUNNUM);   subplot(hp(k*2));  hold on;    % pitch

                    ptc = [Dundun_s(k).rawpitch];     ptcsc = (ptc-(prctile(ptc,2))) / (prctile(ptc,98)-(prctile(ptc,2)));          ptcsc = hpfilter(ptcsc,10);

                    pl=plot(ptcsc,'r','linewidth',.05);     pl.Color(4) = 0.25;   box on
                    plot(N_noteix,N_ptcvect,'.r');              
                    axis([0 12000 0 1.1]); 

            end

            figure(11+100*TESTRUNNUM); sgtitle('Note amplitude. Music-like left, speech-like right')
            figure(12+100*TESTRUNNUM); sgtitle('Note pitch. Music-like left, speech-like right')
            
    end

%% Now CREATING VECTORS FOR ALL NOTE AMP, ENT, PITCH  

        % first put ALL amp and freq of all songs into one long vector allNoteAmp_M / allNoteFreq_M and allNoteAmp_s / allNoteFreq_s. For backdrop scatter.
        
                allNoteAmp_M = [];          allNotePitch_M = [];        allNoteEnt_M = [];   
                allNoteAmp_s = [];          allNotePitch_s = [];        allNoteEnt_s = [];
                allNoteAmpDiff_M = [];      allNotePitchDiff_M = [];    allNoteEntDiff_M = [];
                allNoteAmpDiff_s = [];      allNotePitchDiff_s = [];    allNoteEntDiff_s = [];     
                ctM=1;  ctMd=1;
                cts=1;  ctsd=1;
                for k = 1:15
                    if IS_SCALED
                        % music
                        amp_M = NoteFeatures_M{k,1}.ampscaled;               pitch_M = NoteFeatures_M{k,1}.ptcscaled;              ent_M = NoteFeatures_M{k,1}.entscaled;
                        % speech
                        amp_s = NoteFeatures_s{k,1}.ampscaled;               pitch_s = NoteFeatures_s{k,1}.ptcscaled;              ent_s = NoteFeatures_s{k,1}.entscaled;  
                    else
                        % music
                        amp_M = NoteFeatures_M{k,1}.amp;                     pitch_M = NoteFeatures_M{k,1}.ptc;                    ent_M = NoteFeatures_M{k,1}.ent;
                        % speech
                        amp_s = NoteFeatures_s{k,1}.amp;                     pitch_s = NoteFeatures_s{k,1}.ptc;                    ent_s = NoteFeatures_s{k,1}.ent;  
                    end
                    allNoteAmp_M(ctM:ctM+numel(amp_M)-1) = amp_M;        allNotePitch_M(ctM:ctM+numel(pitch_M)-1) = pitch_M;   allNoteEnt_M(ctM:ctM+numel(ent_M)-1) = ent_M;
                    allNoteAmp_s(cts:cts+numel(amp_s)-1) = amp_s;        allNotePitch_s(cts:cts+numel(pitch_s)-1) = pitch_s;   allNoteEnt_s(cts:cts+numel(ent_s)-1) = ent_s;
                    
                    ctM = ctM+numel(amp_M);
                    cts = cts+numel(amp_s);
                    % ampdiff, entdiff, & pitchdiff for M and s
                    allNoteAmpDiff_M(ctMd:ctMd+numel(diff(amp_M))-1) = diff(amp_M);   allNotePitchDiff_M(ctMd:ctMd+numel(diff(pitch_M))-1) = diff(pitch_M);     allNoteEntDiff_M(ctMd:ctMd+numel(diff(ent_M))-1) = diff(ent_M);
                    allNoteAmpDiff_s(ctsd:ctsd+numel(diff(amp_s))-1) = diff(amp_s);   allNotePitchDiff_s(ctsd:ctsd+numel(diff(pitch_s))-1) = diff(pitch_s);     allNoteEntDiff_s(ctsd:ctsd+numel(diff(ent_s))-1) = diff(ent_s);
                    ctMd = ctMd+numel(diff(amp_M));
                    ctsd = ctsd+numel(diff(amp_s));
                end
    
%% looks like NOTE AMPLITUDE changes more up and down in music-like stimuli!. 

NoteFeatures_struct_M = NoteFeatures_M;     NoteFeatures_struct_s = NoteFeatures_s;     
if IS_SCALED
    fieldname_string = 'ampscaled';  ylabelstring = 'scaled Note intensity';
else
    fieldname_string = 'amp';  ylabelstring = 'Note intensity';
end
FIGURENUMBER = 22+100*TESTRUNNUM;
plot_M_s_scatter_means_STE(NoteFeatures_struct_M,NoteFeatures_struct_s,fieldname_string,ylabelstring,FIGURENUMBER,IS_SCALED)
title('Dundun note intensity')

yl=ylim;
        
        %%   same data but as a histogram (or hist with pdf)
        
        data_M = allNoteAmp_M;    data_s = allNoteAmp_s;    nbins_M = 34;  nbins_s = 30;  FIGURENUMBER = 221;   
        if IS_SCALED, 
            bandwidth = .02;  xstring = 'intensity (scaled)';   
        else
            bandwidth = .7; xstring = 'intensity';   
        end
        plot_M_s_hist_with_density(data_M,data_s,nbins_M,nbins_s,xstring,bandwidth,FIGURENUMBER);
        xlim(yl)
        set(gca,'YDir','reverse');
        set(gcf,'Position',[398 282 134 392])

        camroll(90)

%         figure(123);    clf;   set(gcf,'Position',[900 607 284 199]);   hold on;
%                 hna_M = histogram(allNoteAmp_M,34,'normalization','pdf'); hna_M.FaceColor = [0 0 1]; hna_M.FaceAlpha = 1; hna_M.LineStyle = 'none';
%                 hna_s = histogram(allNoteAmp_s,30,'normalization','pdf'); hna_s.EdgeColor = [1 0 0]; hna_s.FaceAlpha = 0; hna_s.LineWidth = 2.2; 
%                 xlabel('intensity'); ylabel('pdf'); box on
                    edges_amp = linspace(min([data_M,data_s]),max([data_M,data_s]),44+1)';
                    [Hstcounts_amp_M,edges_amp] = histcounts(data_M,edges_amp,'Normalization','pdf');
                    [Hstcounts_amp_s,edges_amp] = histcounts(data_s,edges_amp,'Normalization','pdf'); % this will be entered in Pauline's table further down

%% is NOTE PITCH also different between speech-like and music-like stimuli? 

NoteFeatures_struct_M = NoteFeatures_M;     NoteFeatures_struct_s = NoteFeatures_s;  
if IS_SCALED
    fieldname_string = 'ptcscaled';  ylabelstring = 'scaled Note pitch';
else
    fieldname_string = 'ptc';  ylabelstring = 'Note pitch [Hz]';
end
FIGURENUMBER = 23+100*TESTRUNNUM;
plot_M_s_scatter_means_STE(NoteFeatures_struct_M,NoteFeatures_struct_s,fieldname_string,ylabelstring,FIGURENUMBER,IS_SCALED)
title('Dundun note pitch')

if IS_SCALED
    yl=ylim;
else
    yl=ylim;
end


        %% same data but as a histogram (or hist & PDF)
        
        data_M = allNotePitch_M;    data_s = allNotePitch_s;    nbins_M = 30;  nbins_s = 30;  FIGURENUMBER = 122;   
        if IS_SCALED 
            bandwidth = .026;  xstring = 'pitch';   
        else
            bandwidth = 4;  xstring = 'pitch [Hz]';   
        end
        plot_M_s_hist_with_density(data_M,data_s,nbins_M,nbins_s,xstring,bandwidth,FIGURENUMBER);
        xlim(yl)
        set(gca,'YDir','reverse');
        set(gcf,'Position',[398 282 134 392])

        camroll(90)
        % real histogram
%         figure(123);    clf;   set(gcf,'Position',[900 607 284 199]);   hold on;
%                 hnp_M = histogram(allNotePitch_M,30,'normalization','pdf'); hnp_M.FaceColor = [0 0 1]; hnp_M.FaceAlpha = 1; hnp_M.LineStyle = 'none';
%                 hnp_s = histogram(allNotePitch_s,30,'normalization','pdf'); hnp_s.EdgeColor = [1 0 0]; hnp_s.FaceAlpha = 0; hnp_s.LineWidth = 2.2; 
%                 xlabel('pitch'); ylabel('pdf'); box on
                    edges_pitch = linspace(min([data_M,data_s]),max([data_M,data_s]),44+1)';
                    [Hstcounts_pitch_M,edges_pitch] = histcounts(data_M,edges_pitch,'Normalization','pdf');
                    [Hstcounts_pitch_s,edges_pitch] = histcounts(data_s,edges_pitch,'Normalization','pdf'); % this will be entered in Pauline's table further down

        
%% NOW TIMBRE/ENTROPY scatter plots 

NoteFeatures_struct_M = NoteFeatures_M;     NoteFeatures_struct_s = NoteFeatures_s;     
if IS_SCALED
    fieldname_string = 'entscaled';  ylabelstring = 'scaled Note entropy';
else
    fieldname_string = 'ent';  ylabelstring = 'Note entropy';
end
FIGURENUMBER = 27+100*TESTRUNNUM;
plot_M_s_scatter_means_STE(NoteFeatures_struct_M,NoteFeatures_struct_s,fieldname_string,ylabelstring,FIGURENUMBER,IS_SCALED)
title('Dundun timbre (note entropy)')

yl = ylim;

        %% same Entropy data but as a histogram
                
        data_M = allNoteEnt_M;      data_s = allNoteEnt_s;    xstring = 'entropy (timbre)';     FIGURENUMBER = 312;   
        if IS_SCALED 
            bandwidth = .026; nbins_M = 36;   nbins_s = 25;   
        else
            bandwidth = 0.06;  nbins_M = 56;   nbins_s = 45;  
        end
        plot_M_s_hist_with_density(data_M,data_s,nbins_M,nbins_s,xstring,bandwidth,FIGURENUMBER);
        if IS_SCALED, xlim(yl); else xlim([-5 0]);  end
        set(gca,'YDir','reverse');
        set(gcf,'Position',[398 282 134 392])

        camroll(90)
                    edges_ent = linspace(min([data_M,data_s]),max([data_M,data_s]),44+1)';
                    [Hstcounts_ent_M,edges_ent] = histcounts(data_M,edges_ent,'Normalization','pdf');
                    [Hstcounts_ent_s,edges_ent] = histcounts(data_s,edges_ent,'Normalization','pdf'); % this will be entered in Pauline's table further down

%         clf;   set(gcf,'Position',[900 607 284 199]);   hold on;
%         
%                 hne_M = histogram(allNoteEnt_M,36,'Normalization','pdf','FaceColor','b','FaceAlpha',.3,'EdgeColor','none');%histogram(allNoteEnt_M,36,'normalization','pdf'); hne_M.FaceColor = [0 0 1]; hne_M.FaceAlpha = 1; hne_M.LineStyle = 'none';
%                 hne_s = histogram(allNoteEnt_s,25,'Normalization','pdf','FaceColor','r','FaceAlpha',.3,'EdgeColor','none');%histcolor1,'EdgeAlpha',.8); ?% histogram(allNoteEnt_s,25,'normalization','pdf'); hne_s.EdgeColor = [1 0 0]; hne_s.FaceAlpha = 0; hne_s.LineWidth = 2.2; 
%                 xlabel('entropy (timbre)'); ylabel('pdf'); box on
%                 % density
%                 [f_M,xi_M] = ksdensity(allNoteEnt_M, 'bandwidth',0.026); %'bandwidth',0.008); %,'support','positive'); 0.008 0.012
%                 [f_s,xi_s] = ksdensity(allNoteEnt_s, 'bandwidth',0.026); %'bandwidth',0.008); %,'support','positive');
%                 pM = plot(xi_M,f_M,'Color','b','linewidth',2); % [.8 .8 .8]);
%                 ps = plot(xi_s,f_s,'Color','r','linewidth',2);
%                 xlim([0 1])

    %% NOW SAME FOR TIMING AS THE FOURTH FEATURE - i.e. individual onset-to-onset times, as histogram and as scatter: 
  
    % TIMING SCATTER PLOTS
    % onset-to-onset durations, scatter
    NoteFeatures_struct_M = NoteFeatures_M;     NoteFeatures_struct_s = NoteFeatures_s;     fieldname_string = 'o2o_dur';  ylabelstring = 'onset2onset duration';
    FIGURENUMBER = 434+100*TESTRUNNUM;
    plot_M_s_scatter_means_STE(NoteFeatures_struct_M,NoteFeatures_struct_s,fieldname_string,ylabelstring,FIGURENUMBER,IS_SCALED)
    title('Dundun onset-to-onset durations')
    ylim([0 700])

    
    % more plots for timing: HISTOGRAMS 
        % onsets
        % Music-like
        ctm1=1;     ctm2=1;
        allDur_M = [];        allRatio_M = [];      allCycledur_M = [];
        for i = 1:numel(NoteFeatures_M)
            
                onsets = NoteFeatures_M{i,1}.onsets;   
                durations = onsets(2:end)-onsets(1:end-1);
                cycledurations = durations(1:end-1)+durations(2:end);
                % calculating ratios
                ratios = durations(1:end-1)./cycledurations;
                allDur_M(ctm1:ctm1+numel(durations)-1) = durations;
                allRatio_M(ctm2:ctm2+numel(ratios)-1) = ratios;
                allCycledur_M(ctm2:ctm2+numel(ratios)-1) = cycledurations;
                ctm1=ctm1+numel(durations);
                ctm2=ctm2+numel(ratios);
        end
    
        % speech-like
        ctm1=1;     ctm2=1;
        allDur_s = [];        allRatio_s = [];      allCycledur_s = [];
        for i = 1:numel(NoteFeatures_s)
            
                onsets = NoteFeatures_s{i,1}.onsets;   
                durations = onsets(2:end)-onsets(1:end-1);
                cycledurations = durations(1:end-1)+durations(2:end);
                % calculating ratios
                ratios = durations(1:end-1)./cycledurations;
                allDur_s(ctm1:ctm1+numel(durations)-1) = durations;
                allRatio_s(ctm2:ctm2+numel(ratios)-1) = ratios;
                allCycledur_s(ctm2:ctm2+numel(ratios)-1) = cycledurations;
                ctm1=ctm1+numel(durations);
                ctm2=ctm2+numel(ratios);
        end
        
        
        % duration
        % histogram and probability density    
        data_M = allDur_M;    data_s = allDur_s;    nbins_M = 30;  nbins_s = 60;  xstring = 'o2o duration';   bandwidth = 8;  FIGURENUMBER = 234;   
        plot_M_s_hist_with_density(data_M,data_s,nbins_M,nbins_s,xstring,bandwidth,FIGURENUMBER);
        xlim([0 700])
        set(gca,'YDir','reverse');
        set(gcf,'Position',[398 282 134 392])

        camroll(90)
% 
%         % histogram only
%         figure(284);  clf;   set(gcf,'Position',[900 607 284 199]); hold on
%                 hd_M = histogram(allDur_M,26,'normalization','pdf'); hd_M.FaceColor = [0 0 1]; hd_M.FaceAlpha = 1; hd_M.LineStyle = 'none';
%                 hd_s = histogram(allDur_s,49,'normalization','pdf'); hd_s.EdgeColor = [1 0 0]; hd_s.FaceAlpha = 0; hd_s.LineWidth = 1.2; hd_s.FaceColor = [1 0 0]; hd_s.FaceAlpha = .3%.4; 
%                 xlim([0 600])
%                 xlabel('onset2onset duration'); ylabel('pdf'); box on
        
        % RATIO
        
        % ratio scatter
        NoteFeatures_struct_M = NoteFeatures_M;     NoteFeatures_struct_s = NoteFeatures_s;     fieldname_string = 'ratio';  ylabelstring = 'ratio interval1/(interval1+interval2)';
        FIGURENUMBER = 435+100*TESTRUNNUM;
        plot_M_s_scatter_means_STE(NoteFeatures_struct_M,NoteFeatures_struct_s,fieldname_string,ylabelstring,FIGURENUMBER,IS_SCALED)
        title('Dundun rhythm: two-interval ratios')
        ylim([0 1])

        
        % hist & pdf 1
        data_M = allRatio_M;    data_s = allRatio_s;    nbins_M = 40;  nbins_s = 40;  xstring = 'ratio';   bandwidth = .01;  FIGURENUMBER = 342;   
        plot_M_s_hist_with_density(data_M,data_s,nbins_M,nbins_s,xstring,bandwidth,FIGURENUMBER);
        xlim([0 1])
        hold on;     yl=ylim;   plot([.5 .5],[0 yl(2)],'k','linewidth',.7)

        
%         % hist & pdf 2
%         histcolor1 = [0 0 1]; histcolor2 = [1 0 0]; 
%         figure(343);  clf;   set(gcf,'Position',[900 607 284 199]); hold on
% %                 hd_M = histogram((allRatio_M),40,'normalization','pdf'); hd_M.FaceColor = [0 0 1]; hd_M.FaceAlpha = 1; hd_M.LineStyle = 'none';
% %                 hd_s = histogram(allRatio_s,40,'normalization','pdf'); hd_s.EdgeColor = [1 0 0]; hd_s.FaceAlpha = 0; hd_s.LineWidth = 1.2; hd_s.FaceColor = [1 0 0]; hd_s.FaceAlpha = .3%.4; 
%                 hd_M = histogram(allRatio_M,40,'Normalization','pdf','FaceColor',histcolor1,'FaceAlpha',.3,'EdgeColor','none');%histcolor1,'EdgeAlpha',.8); % 130 or 153 bins for string quartet?
%                 hd_s = histogram(allRatio_s,40,'Normalization','pdf','FaceColor',histcolor2,'FaceAlpha',.3,'EdgeColor','none');%,'EdgeAlpha',.8); % 130 or 153 bins for string quartet?
%                 xlabel('ratio'); ylabel('pdf'); box on
%                 [f_M,xi_M] = ksdensity(allRatio_M, 'bandwidth',0.01); %'bandwidth',0.008); %,'support','positive'); 0.008 0.012
%                 [f_s,xi_s] = ksdensity(allRatio_s, 'bandwidth',0.01); %'bandwidth',0.008); %,'support','positive');
%         % line in the middle
%         plot([.5 .5],[0 7],'k','linewidth',.76)
%         % plotting distributions 
%         pM = plot(xi_M,f_M,'Color','b','linewidth',2); % [.8 .8 .8]);
%         ps = plot(xi_s,f_s,'Color','r','linewidth',2);
%         xlim([0 1])
                
%% IMMEDIATE (FIRST ORDER) DYNAMICS: feature of note n plotted against feature of note n+1.
% to check whether amplitude changes more up and down in M-stimuli, plotting note-amp as a function of previous-note-amp.     
    
%     figure(24+100*TESTRUNNUM);  clf;   set(gcf,'Position',[1 107 1680 794])
%    
%     % amplitude (ALWAYS USE SCALED TO KEEP COMPARABLE BETWEEN RECORDINGS!!)
%         % Music-like
%         for i = 1:numel(NoteFeatures_M)
%             
%             amp = NoteFeatures_M{i,1}.ampscaled;      ent = NoteFeatures_M{i,1}.entscaled;      ptc = NoteFeatures_M{i,1}.ptcscaled;  
%             
%             subplot(2,3,1); hold on;    plot(amp(1:end-1),amp(2:end),'.b')
%             subplot(2,3,2); hold on;    plot(ent(1:end-1),ent(2:end),'xb')
%             subplot(2,3,3); hold on;    plot(ptc(1:end-1),ptc(2:end),'db')            
%         end
%     
%     
%         % speech-like
%         for i = 1:numel(NoteFeatures_s)
%             
%             amp = NoteFeatures_s{i,1}.ampscaled;      ent = NoteFeatures_s{i,1}.entscaled;      ptc = NoteFeatures_s{i,1}.ptcscaled;  
%             
%             subplot(2,3,4); hold on;    plot(amp(1:end-1),amp(2:end),'.r');    box on;     xlabel('note x');       ylabel('note x+1');      axis([0 1.1 0 1.1]);      %if USE_NOTE_PEAK,   axis([0 1.5 0 1.5]);    end    axis([0 1.1 0 1.1]);      
%             subplot(2,3,5); hold on;    plot(ent(1:end-1),ent(2:end),'xr');    box on;     xlabel('note x');       ylabel('note x+1');      axis([0 1.1 0 1.1]);      %if USE_NOTE_PEAK,   axis([0 1.2 0 1.2]);    end    axis([0 1.1 0 1.1]);    
%             subplot(2,3,6); hold on;    plot(ptc(1:end-1),ptc(2:end),'dr');    box on;     xlabel('note x');       ylabel('note x+1');      axis([0 1.1 0 1.1]);      %if USE_NOTE_PEAK,   axis([0 1.2 0 1.2]);    end    axis([0 1.1 0 1.1]);         
%         
%         end
%         
%         subplot(2,3,1);     hold on;    box on;     xlabel('note x');       ylabel('note x+1');        title('Mean note amplitude');    axis([0 1.1 0 1.1]);       %if USE_NOTE_PEAK,   axis([0 1.5 0 1.5]);    end
%         subplot(2,3,2);     hold on;    box on;     xlabel('note x');       ylabel('note x+1');        title('Mean note entropy');      axis([0 1.1 0 1.1]);       %if USE_NOTE_PEAK,   axis([0 1.2 0 1.2]);    end
%         subplot(2,3,3);     hold on;    box on;     xlabel('note x');       ylabel('note x+1');        title('Mean note pitch');        axis([0 1.1 0 1.1]);       %if USE_NOTE_PEAK,   axis([0 1.2 0 1.2]);    end  
%    
%         sgtitle('Adjacent note features')
%         
%         % RESULT AMPLITUDE: 
%         % Overall amplitude range is higher in M. More diverse, farther from the diagonal!! More movement in amplitude. In s, less movement - closer to diagonal!
%         % TO DO: quantify distance from diagonal, i.e. diff(ampscaled), but in relative, not absolute terms???
%         % RESULT PITCH: 
%         % perhaps there is MORE variability in the s stimuli in pitch: They don't concentrate as heavily around the diagonal. They do seem more concentrated around diagonal
%         % in the M stimuli. TEST THIS!! TO DO: calculate distance from diagonal, and/or diff(pitchscaled), i.e. distribution of pitch intervals... best would be again to use
%         % relative, not absolute, intervals.
%         % % RESULT ENTROPY:
%         % TO DO: Test whether distance-from-diagonal is lower in s than in M stimuli.
  
    % Dynamics of AMPLITUDE, ENTROPY, PITCH, plotting feature of note n vs feature of note n+1.
    figure(28+100*TESTRUNNUM);  clf;   set(gcf,'Position',[351 728 486 220])   
    figure(29+100*TESTRUNNUM);  clf;   set(gcf,'Position',[553 632 486 220])
    figure(30+100*TESTRUNNUM);  clf;   set(gcf,'Position',[700 586 486 220])
   
    % (ALWAYS USE SCALED FEATURES TO KEEP COMPARABLE BETWEEN RECORDINGS!!)
        % Music-like
        for i = 1:numel(NoteFeatures_M)
            
            amp = NoteFeatures_M{i,1}.ampscaled;      ent = NoteFeatures_M{i,1}.entscaled;      ptc = NoteFeatures_M{i,1}.ptcscaled;  
            
            figure(28+100*TESTRUNNUM);  hold on;
                subplot(1,2,1); hold on;     if i==1,    plot([0 1.1], [0 1.1], '-k');   end
                plot(amp(1:end-1),amp(2:end),'.b') 
            figure(29+100*TESTRUNNUM); hold on;    
                subplot(1,2,1); hold on;     if i==1,    plot([0 1.1], [0 1.1], '-k');   end
                plot(ent(1:end-1),ent(2:end),'.b') 
            figure(30+100*TESTRUNNUM); hold on;    
                subplot(1,2,1); hold on;     if i==1,    plot([0 1.1], [0 1.1], '-k');   end
                plot(ptc(1:end-1),ptc(2:end),'.b') 
 
        end
        
        % speech-like
        for i = 1:numel(NoteFeatures_s)
            
            amp = NoteFeatures_s{i,1}.ampscaled;      ent = NoteFeatures_s{i,1}.entscaled;      ptc = NoteFeatures_s{i,1}.ptcscaled;  
            
            figure(28+100*TESTRUNNUM); hold on;    
                subplot(1,2,2); hold on;     if i==1,    plot([0 1.1], [0 1.1], '-k');   end
                plot(amp(1:end-1),amp(2:end),'.r');    box on;     xlabel('note x');       ylabel('note x+1');      axis([0 1.1 0 1.1]);      %if USE_NOTE_PEAK,   axis([0 1.5 0 1.5]);    end    axis([0 1.1 0 1.1]);      
            figure(29+100*TESTRUNNUM); hold on;    
                subplot(1,2,2); hold on;     if i==1,    plot([0 1.1], [0 1.1], '-k');   end
                plot(ent(1:end-1),ent(2:end),'.r');    box on;     xlabel('note x');       ylabel('note x+1');      axis([0 1.1 0 1.1]);      %if USE_NOTE_PEAK,   axis([0 1.2 0 1.2]);    end    axis([0 1.1 0 1.1]);    
            figure(30+100*TESTRUNNUM); hold on;    
                subplot(1,2,2); hold on;     if i==1,    plot([0 1.1], [0 1.1], '-k');   end
                plot(ptc(1:end-1),ptc(2:end),'.r');    box on;     xlabel('note x');       ylabel('note x+1');      axis([0 1.1 0 1.1]);      %if USE_NOTE_PEAK,   axis([0 1.2 0 1.2]);    end    axis([0 1.1 0 1.1]);         
        
        end
        
        figure(28+100*TESTRUNNUM); subplot(1,2,1);  box on;     
            xlabel('note x');   ylabel('note x+1');   axis([0 1.1 0 1.1]);    sgtitle('intensity dynamics')   %if USE_NOTE_PEAK,   axis([0 1.5 0 1.5]);    end
        figure(29+100*TESTRUNNUM); subplot(1,2,1);  box on;     
            xlabel('note x');   ylabel('note x+1');   axis([0 1.1 0 1.1]);    sgtitle('timbre (entropy) dynamics')      %if USE_NOTE_PEAK,   axis([0 1.2 0 1.2]);    end
        figure(30+100*TESTRUNNUM); subplot(1,2,1);  box on;     
            xlabel('note x');   ylabel('note x+1');   axis([0 1.1 0 1.1]);     sgtitle('pitch dynamics')     %if USE_NOTE_PEAK,   axis([0 1.2 0 1.2]);    end  
   
    %% THE FIRST-ORDER DYNAMIC OF RHYTHM: NOTE ONSET-TO-ONSET DURATION. Ratio (see above) is a similar measure, but relative. Also plotting absolute values: 
    
    %     % Dynamics of AMPLITUDE, ENTROPY, PITCH, plotting feature of note n vs feature of note n+1.
    figure(567+100*TESTRUNNUM);  clf;   set(gcf,'Position',[351 728 486 220])   
    figure(675+100*TESTRUNNUM);  clf;   set(gcf,'Position',[351 728 486 220])   
   
        % Music-like
        for i = 1:numel(NoteFeatures_M)
            
            o2o = NoteFeatures_M{i,1}.o2o_dur;  
            ratio = NoteFeatures_M{i,1}.ratio;  
            
            figure(567+100*TESTRUNNUM);  hold on;
                subplot(1,2,1); hold on;     if i==1,    plot([0 750], [0 750], '-k');   end
                scatter(o2o(1:end-2),o2o(2:end-1),5,'markeredgecolor','b','markeredgealpha',.5);    %box on;     xlabel('interval x');       ylabel('interval x+1');      axis([0 1.1 0 1.1]);      %if USE_NOTE_PEAK,   axis([0 1.5 0 1.5]);    end    axis([0 1.1 0 1.1]);       
            figure(675+100*TESTRUNNUM); hold on;    
                subplot(1,2,1); hold on;     if i==1,    plot([0 1.1], [0 1.1], '-k');   end
                scatter(ratio(1:end-3),ratio(2:end-2),5,'markeredgecolor','b','markeredgealpha',.5);%box on;     xlabel('ratio x');       ylabel('ratio x+1');      axis([0 1.1 0 1.1]);      %if USE_NOTE_PEAK,   axis([0 1.5 0 1.5]);    end    axis([0 1.1 0 1.1]);        
        end
        
        % speech-like
        for i = 1:numel(NoteFeatures_s)
            
            o2o = NoteFeatures_s{i,1}.o2o_dur;      
            ratio = NoteFeatures_s{i,1}.ratio;     
            
            figure(567+100*TESTRUNNUM); hold on;    
                subplot(1,2,2); hold on;     if i==1,    plot([0 750], [0 750], '-k');   end
                scatter(o2o(1:end-1),o2o(2:end),5,'markeredgecolor','r','markeredgealpha',.5);      %box on;     xlabel('note x');       ylabel('note x+1');      axis([0 1.1 0 1.1]);      %if USE_NOTE_PEAK,   axis([0 1.5 0 1.5]);    end    axis([0 1.1 0 1.1]);        
            figure(675+100*TESTRUNNUM); hold on;    
                subplot(1,2,2); hold on;     if i==1,    plot([0 1.1], [0 1.1], '-k');   end
                scatter(ratio(1:end-3),ratio(2:end-2),5,'markeredgecolor','r','markeredgealpha',.5);%box on;     xlabel('note x');       ylabel('note x+1');      axis([0 1.1 0 1.1]);      %if USE_NOTE_PEAK,   axis([0 1.5 0 1.5]);    end    axis([0 1.1 0 1.1]);          
        end
        
        figure(567+100*TESTRUNNUM); 
            subplot(1,2,1);  hold on;  box on;   xlabel('interval x');   ylabel('interval x+1');   axis([0 600 0 600]); %axis([0 350 0 350]);    
            subplot(1,2,2);  hold on;  box on;   xlabel('interval x');   ylabel('interval x+1');   axis([0 600 0 600]);    
            sgtitle('interval dynamics')   
        figure(675+100*TESTRUNNUM); 
            subplot(1,2,1);  hold on;  box on;     xlabel('ratio x');   ylabel('ratio x+1');   axis([0 1 0 1]);    
            subplot(1,2,2);  hold on;  box on;     xlabel('ratio x');   ylabel('ratio x+1');   axis([0 1 0 1]);    
            sgtitle('ratio dynamics')     

    
            
    %% DIFF: Dynamics of AMPLITUDE, ENTROPY, PITCH, O2O-TIME, RATIO in terms of difference between adjacent notes.
    IS_FOLD = true;
%     #WEITERMACHEN
%         figure(41+100*TESTRUNNUM);  clf;   set(gcf,'Position',[700 607 284 199])   % for real histogram
%         figure(42+100*TESTRUNNUM);  clf;   set(gcf,'Position',[800 607 284 199])   % for real histogram
%         figure(43+100*TESTRUNNUM);  clf;   set(gcf,'Position',[900 607 284 199])   % for real histogram
        
        % music-like
        all_N_amp_M = nan(numel(NoteFeatures_M)*100,1);       all_N_ent_M = nan(numel(NoteFeatures_M)*100,1);       all_N_ptc_M = nan(numel(NoteFeatures_M)*100,1);
        all_N_O2O_M = nan(numel(NoteFeatures_M)*100,1);       all_N_ratio_M = nan(numel(NoteFeatures_M)*100,1);
        ctm = 1;
        for i = 1:numel(NoteFeatures_M)
            N_amp = NoteFeatures_M{i,1}.ampscaled;      all_N_amp_M(ctm:ctm+numel(N_amp)) = [N_amp;nan]; % adding a nan at the end, so diff stays meaningful
            N_ent = NoteFeatures_M{i,1}.entscaled;      all_N_ent_M(ctm:ctm+numel(N_ent)) = [N_ent;nan]; % adding a nan at the end, so diff stays meaningful
            N_ptc = NoteFeatures_M{i,1}.ptcscaled;      all_N_ptc_M(ctm:ctm+numel(N_ptc)) = [N_ptc;nan]; % adding a nan at the end, so diff stays meaningful
            N_O2O = NoteFeatures_M{i,1}.o2o_dur;        all_N_O2O_M(ctm:ctm+numel(N_O2O)) = [N_O2O;nan]; % adding a nan at the end, so diff stays meaningful
            N_ratio = NoteFeatures_M{i,1}.ratio;        all_N_ratio_M(ctm:ctm+numel(N_ratio)) = [N_ratio;nan]; % adding a nan at the end, so diff stays meaningful
            ctm = ctm+numel(N_amp)+1;
        end
        all_N_amp_M(ctm:end) = [];        all_N_ent_M(ctm:end) = [];        all_N_ptc_M(ctm:end) = [];
        all_N_ratio_M(ctm:end) = [];        all_N_O2O_M(ctm:end) = [];
        
        % speech-like
        all_N_amp_s = nan(numel(NoteFeatures_s)*100,1);       all_N_ent_s = nan(numel(NoteFeatures_s)*100,1);       all_N_ptc_s = nan(numel(NoteFeatures_s)*100,1);
        all_N_O2O_s = nan(numel(NoteFeatures_s)*100,1);       all_N_ratio_s = nan(numel(NoteFeatures_s)*100,1);
        cts = 1;
        for i = 1:numel(NoteFeatures_s)
            N_amp = NoteFeatures_s{i,1}.ampscaled;      all_N_amp_s(cts:cts+numel(N_amp)) = [N_amp;nan]; % adding a nan at the end, so diff stays meaningful
            N_ent = NoteFeatures_s{i,1}.entscaled;      all_N_ent_s(cts:cts+numel(N_ent)) = [N_ent;nan]; % adding a nan at the end, so diff stays meaningful
            N_ptc = NoteFeatures_s{i,1}.ptcscaled;      all_N_ptc_s(cts:cts+numel(N_ptc)) = [N_ptc;nan]; % adding a nan at the end, so diff stays meaningful
            N_O2O = NoteFeatures_s{i,1}.o2o_dur;        all_N_O2O_s(cts:cts+numel(N_O2O)) = [N_O2O;nan]; % adding a nan at the end, so diff stays meaningful
            N_ratio = NoteFeatures_s{i,1}.ratio;        all_N_ratio_s(cts:cts+numel(N_ratio)) = [N_ratio;nan]; % adding a nan at the end, so diff stays meaningful
            cts = cts+numel(N_amp)+1;
        end
        all_N_amp_s(cts:end) = [];        all_N_ent_s(cts:end) = [];        all_N_ptc_s(cts:end) = [];
        all_N_O2O_s(cts:end) = [];        all_N_ratio_s(cts:end) = [];
        
            % figure(43+100*TESTRUNNUM);  clf;   set(gcf,'Position',[900 607 284 199])
            
%             % real histograms
%             figure(41+100*TESTRUNNUM);  hold on;
%                 ha_M = histogram(diff(all_N_amp_M),34,'normalization','pdf'); ha_M.FaceColor = [0 0 1]; ha_M.FaceAlpha = 1; ha_M.LineStyle = 'none';
%                 ha_s = histogram(diff(all_N_amp_s),30,'normalization','pdf'); ha_s.EdgeColor = [1 0 0]; ha_s.FaceAlpha = 0; ha_s.LineWidth = 2.2; 
%                 xlabel('\Delta intensity'); ylabel('pdf'); box on
% %                 title('intensity dynamics, adjacent notes'); box on
%             figure(42+100*TESTRUNNUM); hold on;    
%                 he_M = histogram(diff(all_N_ent_M),40,'normalization','pdf'); he_M.FaceColor = [0 0 1]; he_M.FaceAlpha = 1; he_M.LineStyle = 'none';
%                 he_s = histogram(diff(all_N_ent_s),24,'normalization','pdf'); he_s.EdgeColor = [1 0 0]; he_s.FaceAlpha = 0; he_s.LineWidth = 2.2;
%                 xlabel('\Delta entropy'); ylabel('pdf'); box on
% %                 title('timbre dynamics, adjacent notes'); box on
%             figure(43+100*TESTRUNNUM); hold on;    
%                 hp_M = histogram(diff(all_N_ptc_M),44,'normalization','pdf'); hp_M.FaceColor = [0 0 1]; hp_M.FaceAlpha = 1; hp_M.LineStyle = 'none';
%                 hp_s = histogram(diff(all_N_ptc_s),33,'normalization','pdf'); hp_s.EdgeColor = [1 0 0]; hp_s.FaceAlpha = 0; hp_s.LineWidth = 2.2;
%                 xlabel('\Delta pitch'); ylabel('pdf'); box on
% %                 title('pitch dynamics, adjacent notes')
       
            % histogram with probability density on top
            % amplitude diff
            if IS_FOLD
                data_M = abs(diff(all_N_amp_M)); data_s = abs(diff(all_N_amp_s)); nbins_M = 34;  nbins_s = 30; xstring = '\Delta intensity';  bandwidth = .009;
            else
                data_M = diff(all_N_amp_M); data_s = diff(all_N_amp_s); nbins_M = 34;  nbins_s = 30; xstring = '\Delta intensity';  bandwidth = .02;
            end
            all_N_amp_diff_M = data_M;  all_N_amp_diff_s = data_s;   
            FIGURENUMBER = 141+200*TESTRUNNUM;
            plot_M_s_hist_with_density(data_M,data_s,nbins_M,nbins_s,xstring,bandwidth,FIGURENUMBER)
                    edges_D_amp = linspace(min([data_M;data_s]),max([data_M;data_s]),35+1)';
                    [Hstcounts_D_amp_M,edges_D_amp] = histcounts(data_M,edges_D_amp,'Normalization','pdf');
                    [Hstcounts_D_amp_s,edges_D_amp] = histcounts(data_s,edges_D_amp,'Normalization','pdf');
            if IS_FOLD, xlim([0 .5]); ylim([0 18]); end

            % entropy diff
            if IS_FOLD
                data_M = abs(diff(all_N_ent_M)); data_s = abs(diff(all_N_ent_s)); nbins_M = 40;  nbins_s = 24; xstring = '\Delta entropy';   bandwidth = .016;
            else
                data_M = diff(all_N_ent_M); data_s = diff(all_N_ent_s); nbins_M = 40;  nbins_s = 24; xstring = '\Delta entropy';    bandwidth = .02;
            end
            all_N_ent_diff_M = data_M;  all_N_ent_diff_s = data_s;
            FIGURENUMBER = 142+200*TESTRUNNUM;
            plot_M_s_hist_with_density(data_M,data_s,nbins_M,nbins_s,xstring,bandwidth,FIGURENUMBER)
                    edges_D_ent = linspace(min([data_M;data_s]),max([data_M;data_s]),40+1)'; 
                    [Hstcounts_D_ent_M,edges_D_ent] = histcounts(data_M,edges_D_ent,'Normalization','pdf');
                    [Hstcounts_D_ent_s,edges_D_ent] = histcounts(data_s,edges_D_ent,'Normalization','pdf');
            if IS_FOLD, xlim([0 1]); end

            % pitch diff
            if IS_FOLD
                data_M = abs(diff(all_N_ptc_M)); data_s = abs(diff(all_N_ptc_s)); nbins_M = 44;  nbins_s = 33; xstring = '\Delta pitch';  bandwidth = .016;
            else
                data_M = diff(all_N_ptc_M); data_s = diff(all_N_ptc_s); nbins_M = 44;  nbins_s = 33; xstring = '\Delta pitch';  bandwidth = .02;
            end
            all_N_ptc_diff_M = data_M;  all_N_ptc_diff_s = data_s;
            FIGURENUMBER = 143+200*TESTRUNNUM;
            plot_M_s_hist_with_density(data_M,data_s,nbins_M,nbins_s,xstring,bandwidth,FIGURENUMBER)
                    edges_D_pitch = linspace(min([data_M;data_s]),max([data_M;data_s]),44+1)';
                    [Hstcounts_D_pitch_M,edges_D_pitch] = histcounts(data_M,edges_D_pitch,'Normalization','pdf');
                    [Hstcounts_D_pitch_s,edges_D_pitch] = histcounts(data_s,edges_D_pitch,'Normalization','pdf');
            if IS_FOLD, xlim([0 1]); end
        
            % O2O-times diff
            if IS_SCALED
                if IS_FOLD
                    data_M = abs(diff(all_N_O2O_M)); data_s = abs(diff(all_N_O2O_s)); nbins_M = 40;  nbins_s = 30; xstring = '\Delta O2O intervals [ms]';  bandwidth = 2;   
                else
                    data_M = diff(all_N_O2O_M); data_s = diff(all_N_O2O_s); nbins_M = 44;  nbins_s = 33; xstring = '\Delta O2O intervals [ms]';  bandwidth = 5;   
                end
            else
                if IS_FOLD
                    data_M = abs(diff(all_N_O2O_M)); data_s = abs(diff(all_N_O2O_s)); nbins_M = 40;  nbins_s = 30; xstring = '\Delta O2O intervals [ms]';  bandwidth = 3;   
                else
                    data_M = diff(all_N_O2O_M); data_s = diff(all_N_O2O_s); nbins_M = 44;  nbins_s = 33; xstring = '\Delta O2O intervals [ms]';  bandwidth = 8;   
                end
            end
            % O2O times can assume extreme values that make no sense and bust the pdf. replacing anything beyond 300 ms with NaN....
            data_M(abs(data_M)>300) = nan;  data_s(abs(data_s)>300) = nan;
            all_N_O2O_diff_M = data_M;  all_N_O2O_diff_s = data_s;
            FIGURENUMBER = 144+200*TESTRUNNUM;
            plot_M_s_hist_with_density(data_M,data_s,nbins_M,nbins_s,xstring,bandwidth,FIGURENUMBER)
                    edges_D_O2O = linspace(-300,300,44+1)';
                    [Hstcounts_D_O2O_M,edges_D_O2O] = histcounts(data_M,edges_D_O2O,'Normalization','pdf');
                    [Hstcounts_D_O2O_s,edges_D_O2O] = histcounts(data_s,edges_D_O2O,'Normalization','pdf');
            if IS_FOLD, xlim([0 300]); end

            % Ratios diff
            data_M = diff(all_N_ratio_M); data_s = diff(all_N_ratio_s); nbins_M = 44;  nbins_s = 44; xstring = '\Delta ratio';  
            all_N_ratio_diff_M = data_M;  all_N_amp_ratio_diff_s = data_s;
            bandwidth = .02;   FIGURENUMBER = 145+200*TESTRUNNUM;
            plot_M_s_hist_with_density(data_M,data_s,nbins_M,nbins_s,xstring,bandwidth,FIGURENUMBER)
                    edges_D_ratio = linspace(min([data_M;data_s]),max([data_M;data_s]),44+1)';
                    [Hstcounts_D_ratio_M,edges_D_ratio] = histcounts(data_M,edges_D_ratio,'Normalization','pdf');
                    [Hstcounts_D_ratio_s,edges_D_ratio] = histcounts(data_s,edges_D_ratio,'Normalization','pdf');
            
            % O2O intervals & ratio, not diff        
            bandwidth = .01;   FIGURENUMBER = 146+200*TESTRUNNUM;  nbins_M = 44;  nbins_s = 44; xstring = 'Ratio';  
            plot_M_s_hist_with_density(all_N_ratio_M,all_N_ratio_s,nbins_M,nbins_s,xstring,bandwidth,FIGURENUMBER)
            edges_ratio = linspace(0,1,44+1)';
                    [Hstcounts_ratio_M,edges_ratio] = histcounts(all_N_ratio_M,edges_ratio,'Normalization','pdf');
                    [Hstcounts_ratio_s,edges_ratio] = histcounts(all_N_ratio_s,edges_ratio,'Normalization','pdf');
                    
            if IS_SCALED
                bandwidth = 6;  
            else
                bandwidth = 10;
            end 
            FIGURENUMBER = 147+200*TESTRUNNUM;  nbins_M = 44;  nbins_s = 44; xstring = 'O2O intervals';  
            plot_M_s_hist_with_density(all_N_O2O_M,all_N_O2O_s,nbins_M,nbins_s,xstring,bandwidth,FIGURENUMBER)
            edges_O2O = linspace(0,300,44+1)';%linspace( min([all_N_O2O_M(abs(all_N_O2O_M)<300);all_N_O2O_s(abs(all_N_O2O_s)<300)]) , max([all_N_O2O_M(abs(all_N_O2O_M)<300);all_N_O2O_s(abs(all_N_O2O_s)<300)]),44+1)';
                    [Hstcounts_O2O_M,edges_O2O] = histcounts(all_N_O2O_M(abs(all_N_O2O_M)<300),edges_O2O,'Normalization','pdf');
                    [Hstcounts_O2O_s,edges_O2O] = histcounts(all_N_O2O_s(abs(all_N_O2O_s)<300),edges_O2O,'Normalization','pdf');

                        fieldnames4Pauline = {'BinEdgeLeft_amp','BinEdgeRight_amp','ProbDens_amp_M','ProbDens_amp_s',...
                                              'BinEdgeLeft_pitch','BinEdgeRight_pitch','ProbDens_pitch_M','ProbDens_pitch_s',...
                                              'BinEdgeLeft_ent','BinEdgeRight_ent','ProbDens_ent_M','ProbDens_ent_s',... 
                                              'BinEdgeLeft_O2O','BinEdgeRight_O2O','ProbDens_O2O_M','ProbDens_O2O_s',...
                                              'BinEdgeLeft_ratio','BinEdgeRight_ratio','ProbDens_ratio_M','ProbDens_ratio_s',...
                                                    'BinEdgeLeft_amp_change','BinEdgeRight_amp_change','ProbDens_amp_change_M','ProbDens_amp_change_s',...
                                                    'BinEdgeLeft_pitch_change','BinEdgeRight_pitch_change','ProbDens_pitch_change_M','ProbDens_pitch_change_s',...
                                                    'BinEdgeLeft_ent_change','BinEdgeRight_ent_change','ProbDens_ent_change_M','ProbDens_ent_change_s',...
                                                    'BinEdgeLeft_O2O_change','BinEdgeRight_O2O_change','ProbDens_O2O_change_M','ProbDens_O2O_change_s',...
                                                    'BinEdgeLeft_ratio_change','BinEdgeRight_ratio_change','ProbDens_ratio_change_M','ProbDens_ratio_change_s'};

                        % padding all columns with nan to same length. Otherwise they won't write to xls.
                        targetlength = max([length(Hstcounts_D_amp_M),length(Hstcounts_D_pitch_M),length(Hstcounts_D_ent_M),length(Hstcounts_D_O2O_M),length(Hstcounts_D_ratio_M),length(Hstcounts_amp_M),length(Hstcounts_ent_M),length(Hstcounts_pitch_M),length(Hstcounts_O2O_M),length(Hstcounts_ratio_M)])+1;
                            % diff of note features
                            edges4P_D_amp = [edges_D_amp;nan(targetlength-numel(edges_D_amp)+1,1)];         Hstcounts4P_D_amp_M = [Hstcounts_D_amp_M,nan(1,targetlength-numel(Hstcounts_D_amp_M))]';       Hstcounts4P_D_amp_s = [Hstcounts_D_amp_s,nan(1,targetlength-numel(Hstcounts_D_amp_s))]';
                            edges4P_D_ent = [edges_D_ent;nan(targetlength-numel(edges_D_ent)+1,1)];         Hstcounts4P_D_ent_M = [Hstcounts_D_ent_M,nan(1,targetlength-numel(Hstcounts_D_ent_M))]';       Hstcounts4P_D_ent_s = [Hstcounts_D_ent_s,nan(1,targetlength-numel(Hstcounts_D_ent_s))]';
                            edges4P_D_pitch = [edges_D_pitch;nan(targetlength-numel(edges_D_pitch)+1,1)];   Hstcounts4P_D_pitch_M = [Hstcounts_D_pitch_M,nan(1,targetlength-numel(Hstcounts_D_pitch_M))]'; Hstcounts4P_D_pitch_s = [Hstcounts_D_pitch_s,nan(1,targetlength-numel(Hstcounts_D_pitch_s))]';
                            edges4P_D_O2O = [edges_D_O2O;nan(targetlength-numel(edges_D_O2O)+1,1)];         Hstcounts4P_D_O2O_M = [Hstcounts_D_O2O_M,nan(1,targetlength-numel(Hstcounts_D_O2O_M))]';       Hstcounts4P_D_O2O_s = [Hstcounts_D_O2O_s,nan(1,targetlength-numel(Hstcounts_D_O2O_s))]';
                            edges4P_D_ratio = [edges_D_ratio;nan(targetlength-numel(edges_D_ratio)+1,1)];   Hstcounts4P_D_ratio_M = [Hstcounts_D_ratio_M,nan(1,targetlength-numel(Hstcounts_D_ratio_M))]'; Hstcounts4P_D_ratio_s = [Hstcounts_D_ratio_s,nan(1,targetlength-numel(Hstcounts_D_ratio_s))]';
                            % note features
                            edges4P_amp = [edges_amp;nan(targetlength-numel(edges_amp)+1,1)];         Hstcounts4P_amp_M = [Hstcounts_amp_M,nan(1,targetlength-numel(Hstcounts_amp_M))]';       Hstcounts4P_amp_s = [Hstcounts_amp_s,nan(1,targetlength-numel(Hstcounts_amp_s))]';
                            edges4P_ent = [edges_ent;nan(targetlength-numel(edges_ent)+1,1)];         Hstcounts4P_ent_M = [Hstcounts_ent_M,nan(1,targetlength-numel(Hstcounts_ent_M))]';       Hstcounts4P_ent_s = [Hstcounts_ent_s,nan(1,targetlength-numel(Hstcounts_ent_s))]';
                            edges4P_pitch = [edges_pitch;nan(targetlength-numel(edges_pitch)+1,1)];   Hstcounts4P_pitch_M = [Hstcounts_pitch_M,nan(1,targetlength-numel(Hstcounts_pitch_M))]'; Hstcounts4P_pitch_s = [Hstcounts_pitch_s,nan(1,targetlength-numel(Hstcounts_pitch_s))]';
                            edges4P_O2O = [edges_O2O;nan(targetlength-numel(edges_O2O)+1,1)];         Hstcounts4P_O2O_M = [Hstcounts_O2O_M,nan(1,targetlength-numel(Hstcounts_O2O_M))]'; Hstcounts4P_O2O_s = [Hstcounts_O2O_s,nan(1,targetlength-numel(Hstcounts_O2O_s))]';
                            edges4P_ratio = [edges_ratio;nan(targetlength-numel(edges_ratio)+1,1)];   Hstcounts4P_ratio_M = [Hstcounts_ratio_M,nan(1,targetlength-numel(Hstcounts_ratio_M))]'; Hstcounts4P_ratio_s = [Hstcounts_ratio_s,nan(1,targetlength-numel(Hstcounts_ratio_s))]';

                        T = table(edges4P_amp(1:end-1),edges4P_amp(2:end),Hstcounts4P_amp_M,Hstcounts4P_amp_s,...               % note features, amp
                                  edges4P_pitch(1:end-1),edges4P_pitch(2:end),Hstcounts4P_pitch_M,Hstcounts4P_pitch_s,...       % note features, pitch
                                  edges4P_ent(1:end-1),edges4P_ent(2:end),Hstcounts4P_ent_M,Hstcounts4P_ent_s,...               % note features, ent
                                  edges4P_O2O(1:end-1),edges4P_O2O(2:end),Hstcounts4P_O2O_M,Hstcounts4P_O2O_s,...               % note features, O2O
                                  edges4P_ratio(1:end-1),edges4P_ratio(2:end),Hstcounts4P_ratio_M,Hstcounts4P_ratio_s,...           % note features, ratio
                                  edges4P_D_amp(1:end-1),edges4P_D_amp(2:end),Hstcounts4P_D_amp_M,Hstcounts4P_D_amp_s,...           % diff amp
                                  edges4P_D_pitch(1:end-1),edges4P_D_pitch(2:end),Hstcounts4P_D_pitch_M,Hstcounts4P_D_pitch_s,...   % diff pitch
                                  edges4P_D_ent(1:end-1),edges4P_D_ent(2:end),Hstcounts4P_D_ent_M,Hstcounts4P_D_ent_s,...           % diff ent
                                  edges4P_D_O2O(1:end-1),edges4P_D_O2O(2:end),Hstcounts4P_D_O2O_M,Hstcounts4P_D_O2O_s,...           % diff O2O
                                  edges4P_D_ratio(1:end-1),edges4P_D_ratio(2:end),Hstcounts4P_D_ratio_M,Hstcounts4P_D_ratio_s,...   % diff ratio
                                  'VariableNames',fieldnames4Pauline);
                              
                        % Writing table to excel sheet for Pauline: All histogram bin data for both feature and feature-diff data.
                            cd(results_dir);
                            % transfer to xls
                            writetable(T,'Dundun_Note&ChangeHistograms_4Pauline.xls','WriteRowNames',true); 
                            
                            % diff data, so Pauline can look at the stats of those directly.
                            
%% kstest Kolmogorov Smirnov Test 

            % ampdiff
            x1 = diff(all_N_amp_M);        x2 = diff(all_N_amp_s);
            [h_diffamp,p_diffamp] = kstest2(x1,x2);
                figure(141+200*TESTRUNNUM); hold on
                if h_diffamp==1,             text(.2,6,['p = ',num2str(p_diffamp,'%2.3f')],'fontsize',11)
                elseif h_diffamp==0,         text(.2,6,{'not significant',['p = ',num2str(p_diffamp)]})
                end
            % entdiff
            x1 = diff(all_N_ent_M);        x2 = diff(all_N_ent_s);
            [h_diffent,p_diffent] = kstest2(x1,x2);
                figure(142+200*TESTRUNNUM); hold on
                if h_diffent==1,             text(.2,2.5,['p = ',num2str(p_diffent,'%2.3f')],'fontsize',11)
                elseif h_diffent==0,         text(.2,2.5,{'n.s., ',['p = ',num2str(p_diffent)]})
                end
            % pitchdiff
            x1 = diff(all_N_ptc_M);        x2 = diff(all_N_ptc_s);
            [h_diffpitch,p_diffpitch] = kstest2(x1,x2);
                figure(143+200*TESTRUNNUM); hold on
                if h_diffpitch==1,             text(.2,2.5,['p = ',num2str(p_diffpitch,'%2.3f')],'fontsize',11)
                elseif h_diffpitch==0,         text(.2,2.5,{'n.s., ',['p = ',num2str(p_diffpitch)]})
                end
            % O2O-diff
            x1 = diff(all_N_O2O_M);        x2 = diff(all_N_O2O_s);
            [h_diffO2O,p_diffO2O] = kstest2(x1,x2);
                figure(144+200*TESTRUNNUM); hold on
                if h_diffO2O==1,             text(100,.007,['p = ',num2str(p_diffO2O,'%2.3f')],'fontsize',11)
                elseif h_diffO2O==0,         text(100,.007,{'n.s., ',['p = ',num2str(p_diffO2O)]})
                end
            % ratio diff
            x1 = diff(all_N_ratio_M);        x2 = diff(all_N_ratio_s);
            [h_diffratio,p_diffratio] = kstest2(x1,x2);
                figure(145+200*TESTRUNNUM); hold on
                if h_diffratio==1,             text(.2,2.5,['p = ',num2str(p_diffratio,'%2.3f')],'fontsize',11)
                elseif h_diffratio==0,         text(.2,2.5,{'n.s., ',['p = ',num2str(p_diffratio)]})
                end
            % ratio, no diff
            x1 = all_N_ratio_M;        x2 = all_N_ratio_s;
            [h_ratio,p_ratio] = kstest2(x1,x2);
                figure(146+200*TESTRUNNUM); hold on
                if h_ratio==1,             text(.2,2.5,['p = ',num2str(p_ratio,'%2.3f')],'fontsize',11)
                elseif h_ratio==0,         text(.2,2.5,{'n.s., ',['p = ',num2str(p_ratio)]})
                end
            % o2o, no diff
            x1 = all_N_O2O_M;        x2 = all_N_O2O_s;
            [h_O2O,p_O2O] = kstest2(x1,x2);
                figure(147+200*TESTRUNNUM); hold on
                if h_O2O==1,             text(.2,2.5,['p = ',num2str(p_O2O,'%2.3f')],'fontsize',11)
                elseif h_O2O==0,         text(.2,2.5,{'n.s., ',['p = ',num2str(p_O2O)]})
                end
            % amp, no diff
            x1 = all_N_amp_M;        x2 = all_N_amp_s;
            [h_amp,p_amp] = kstest2(x1,x2);
                figure(148+200*TESTRUNNUM); hold on
                if h_amp==1,             text(.2,2.5,['p = ',num2str(p_amp,'%2.3f')],'fontsize',11)
                elseif h_amp==0,         text(.2,2.5,{'n.s., ',['p = ',num2str(p_amp)]})
                end
            % ent, no diff
            x1 = all_N_ent_M;        x2 = all_N_ent_s;
            [h_ent,p_ent] = kstest2(x1,x2);
                figure(149+200*TESTRUNNUM); hold on
                if h_ent==1,             text(.2,2.5,['p = ',num2str(p_ent,'%2.3f')],'fontsize',11)
                elseif h_ent==0,         text(.2,2.5,{'n.s., ',['p = ',num2str(p_ent)]})
                end
            % pitch, no diff
            x1 = all_N_ptc_M;        x2 = all_N_ptc_s;
            [h_ptc,p_ptc] = kstest2(x1,x2);
                figure(150+200*TESTRUNNUM); hold on
                if h_ptc==1,             text(.2,2.5,['p = ',num2str(p_ptc,'%2.3f')],'fontsize',11)
                elseif h_ptc==0,         text(.2,2.5,{'n.s., ',['p = ',num2str(p_ptc)]})
                end
            
            fprintf('ks test results: \n amplitude: p = %2.3f \n entropy: p = %2.3f \n pitch: p = %2.3f \n onset-onset intervals: p = %2.3f \n ratio: p = %2.3f \n diff_amplitude: p = %2.3f \n diff_entropy: p = %2.3f \n diff_pitch: p = %2.3f \n diff_onset-to-onset interval: p = %2.3f \n diff_ratio: p = %2.3f \n',p_amp,p_ent,p_ptc,p_O2O,p_ratio,p_diffamp,p_diffent,p_diffpitch,p_diffO2O,p_diffratio);
            
            
%% NOW PLOTTING ADJACENT NOTE MOVEMENTS (Intervals)

    figure(25+100*TESTRUNNUM);  clf;   set(gcf,'Position',[1 107 1680 794])
            subplot(2,3,1); hold on;    plot([-1 1],[0 0],'Color',[.6 .6 .6]);      plot([0 0],[-1 1],'Color',[.6 .6 .6])
            subplot(2,3,2); hold on;    plot([-1 1],[0 0],'Color',[.6 .6 .6]);      plot([0 0],[-1 1],'Color',[.6 .6 .6])
            subplot(2,3,3); hold on;    plot([-2 2],[0 0],'Color',[.6 .6 .6]);      plot([0 0],[-2 2],'Color',[.6 .6 .6])
            subplot(2,3,4); hold on;    plot([-1 1],[0 0],'Color',[.6 .6 .6]);      plot([0 0],[-1 1],'Color',[.6 .6 .6])
            subplot(2,3,5); hold on;    plot([-1 1],[0 0],'Color',[.6 .6 .6]);      plot([0 0],[-1 1],'Color',[.6 .6 .6])
            subplot(2,3,6); hold on;    plot([-2 2],[0 0],'Color',[.6 .6 .6]);      plot([0 0],[-2 2],'Color',[.6 .6 .6])
   
    % amplitude, frequency, and entropy
        % Music-like. USE SCALED, TO COMPARE BETWEEN RECORDINGS
        for i = 1:numel(NoteFeatures_M)
            
            amp = NoteFeatures_M{i,1}.ampscaled;      ent = NoteFeatures_M{i,1}.entscaled;      ptc = NoteFeatures_M{i,1}.ptcscaled;  
            ampdiff = diff(amp);                      entdiff = diff(ent);                      ptcdiff = diff(ptc);  
            
            subplot(2,3,1); hold on;    plot(ampdiff(1:end-1),ampdiff(2:end),'.b')
            subplot(2,3,2); hold on;    plot(entdiff(1:end-1),entdiff(2:end),'xb')
            subplot(2,3,3); hold on;    plot(ptcdiff(1:end-1),ptcdiff(2:end),'db')            
        end
    
    
        % speech-like
        for i = 1:numel(NoteFeatures_s)
            
            amp = NoteFeatures_s{i,1}.ampscaled;      ent = NoteFeatures_s{i,1}.entscaled;      ptc = NoteFeatures_s{i,1}.ptcscaled;  
            ampdiff = diff(amp);                      entdiff = diff(ent);                      ptcdiff = diff(ptc);  
            
            subplot(2,3,4); hold on;    plot(ampdiff(1:end-1),ampdiff(2:end),'.r');    box on;      axis([-.5 .5 -.5 .5]);    
                                        xlabel('note x');           ylabel('note x+1');         
            subplot(2,3,5); hold on;    plot(entdiff(1:end-1),entdiff(2:end),'xr');    box on;      axis([-.7 .7 -.7 .7]);    
                                        xlabel('note x');           ylabel('note x+1');         
            subplot(2,3,6); hold on;    plot(ptcdiff(1:end-1),ptcdiff(2:end),'dr');  box on;        axis([-.8 .8 -.8 .8]);    
                                        xlabel('note x');           ylabel('note x+1');                
        
        end
        
        subplot(2,3,1);         hold on;    box on;                         axis([-.5 .5 -.5 .5]);    
        xlabel('interval x');               ylabel('interval x+1');         title('Amplitude intervals')
        subplot(2,3,2);         hold on;    box on;                         axis([-.7 .7 -.7 .7]);    
        xlabel('interval x');               ylabel('interval x+1');         title('Entropy intervals')
        subplot(2,3,3);         hold on;    box on;                         axis([-.8 .8 -.8 .8]);    
        xlabel('interval x');               ylabel('interval x+1');         title('Pitch intervals')     
   
        sgtitle('Adjacent intervals / note movements')

        
%% plot freq intervals as a function of amp intervals. Then check how drummers cycle through them. 
        
figure(31+100*TESTRUNNUM); clf;    %set(gcf,'Position',)

        for i = 1:numel(NoteFeatures_M)            
            amp = NoteFeatures_M{i,1}.ampscaled;      ptc = NoteFeatures_M{i,1}.ptcscaled;  
            ampdiff = diff(amp);                      ptcdiff = diff(ptc);   
            
            subplot(2,2,1); hold on;    
                plot([-1 1],[0 0],'Color',[.6 .6 .6]);      plot([0 0],[-1 1],'Color',[.6 .6 .6])
                plot(amp,ptc,'.b');         xlabel('amplitude');           ylabel('pitch');                title('Music-like amp+pitch');              box on
                axis([0 1.25 0 1.25])
            subplot(2,2,3); hold on;    
                plot([-1 1],[0 0],'Color',[.6 .6 .6]);      plot([0 0],[-1 1],'Color',[.6 .6 .6])
                plot(ampdiff,ptcdiff,'xb'); xlabel('amplitude interval');  ylabel('pitch interval');       title('Music-like amp+pitch intervals');    box on  
                axis([-.8 .8 -.8 .8])
        end
    
    
        % speech-like
        for i = 1:numel(NoteFeatures_s)            
            amp = NoteFeatures_s{i,1}.ampscaled;      ptc = NoteFeatures_s{i,1}.ptcscaled;  
            ampdiff = diff(amp);                 ptcdiff = diff(ptc); 
            
            subplot(2,2,2); hold on;    
                plot([-1 1],[0 0],'Color',[.6 .6 .6]);      plot([0 0],[-1 1],'Color',[.6 .6 .6])
                plot(amp,ptc,'.r');          xlabel('amplitude');           ylabel('pitch');            title('Speech-like amp+pitch');             box on  
                axis([0 1.25 0 1.25])
            subplot(2,2,4); hold on;    
                plot([-1 1],[0 0],'Color',[.6 .6 .6]);      plot([0 0],[-1 1],'Color',[.6 .6 .6])
                plot(ampdiff,ptcdiff,'xr');  xlabel('amplitude interval');  ylabel('pitch interval');   title('Speech-like amp+pitch intervals');   box on  
                axis([-.8 .8 -.8 .8])
        end
                
    
        %% EACH SONG AS A CONNECTED LINE ON GRAY BACKDROP OF AMP/FREQ FEATURES.

%         % first put ALL amp and freq of all songs into one long vector allNoteAmp_M / allNoteFreq_M and allNoteAmp_s / allNoteFreq_s. For backdrop scatter.
%         
%                 allNoteAmp_M = [];          allNotePitch_M = [];        allNoteEnt_M = [];   
%                 allNoteAmp_s = [];          allNotePitch_s = [];        allNoteEnt_s = [];
%                 allNoteAmpDiff_M = [];      allNotePitchDiff_M = [];    allNoteEntDiff_M = [];
%                 allNoteAmpDiff_s = [];      allNotePitchDiff_s = [];    allNoteEntDiff_s = [];     
%                 ctM=1;  ctMd=1;
%                 cts=1;  ctsd=1;
%                 for k = 1:15
%                     % music
%                     amp_M = NoteFeatures_M{k,1}.ampscaled;               pitch_M = NoteFeatures_M{k,1}.ptcscaled;              ent_M = NoteFeatures_M{k,1}.entscaled;  
%                     allNoteAmp_M(ctM:ctM+numel(amp_M)-1) = amp_M;        allNotePitch_M(ctM:ctM+numel(pitch_M)-1) = pitch_M;   allNoteEnt_M(ctM:ctM+numel(ent_M)-1) = ent_M;
%                     % speech
%                     amp_s = NoteFeatures_s{k,1}.ampscaled;               pitch_s = NoteFeatures_s{k,1}.ptcscaled;              ent_s = NoteFeatures_s{k,1}.entscaled;  
%                     allNoteAmp_s(cts:cts+numel(amp_s)-1) = amp_s;        allNotePitch_s(cts:cts+numel(pitch_s)-1) = pitch_s;   allNoteEnt_s(cts:cts+numel(ent_s)-1) = ent_s;
%                     ctM = ctM+numel(amp_M);
%                     cts = cts+numel(amp_s);
%                     % ampdiff, entdiff, & pitchdiff for M and s
%                     allNoteAmpDiff_M(ctMd:ctMd+numel(diff(amp_M))-1) = diff(amp_M);   allNotePitchDiff_M(ctMd:ctMd+numel(diff(pitch_M))-1) = diff(pitch_M);     allNoteEntDiff_M(ctMd:ctMd+numel(diff(ent_M))-1) = diff(ent_M);
%                     allNoteAmpDiff_s(ctsd:ctsd+numel(diff(amp_s))-1) = diff(amp_s);   allNotePitchDiff_s(ctsd:ctsd+numel(diff(pitch_s))-1) = diff(pitch_s);     allNoteEntDiff_s(ctsd:ctsd+numel(diff(ent_s))-1) = diff(ent_s);
%                     ctMd = ctMd+numel(diff(amp_M));
%                     ctsd = ctsd+numel(diff(amp_s));
%                 end

        % now plotting the backdrop scatters, and the individual songs on top.
        
            for k = 1:15 % one plot for each individual song. 15 music-like ones in blue, 15 speech-like ones in red.
                    
                figure(k+240+100*TESTRUNNUM); clf;

                % MUSIC: plotting gray backdrop for Note features and Interval features (amp & pitch) - MUSIC-LIKE:

                    subplot(2,2,1); hold on;    %plot([-2 2],[0 0],'Color',[.5 .5 .5]);  plot([0 0],[-2 2],'Color',[.5 .5 .5]); 
                    plot(allNoteAmp_M,allNotePitch_M,'.','Color',[.75 .75 .75]);    xlabel('amplitude');    ylabel('pitch');   box on;      
                    title('Music-like');    
                                                
                    subplot(2,2,3); hold on;    plot([-2 2],[0 0],'Color',[.5 .5 .5]);  plot([0 0],[-2 2],'Color',[.5 .5 .5]); 
                    plot(allNoteAmpDiff_M,allNotePitchDiff_M,'.','Color',[.75 .75 .75]); xlabel('amplitude interval');  ylabel('pitch interval');   box on 
%                     title('Music-like amp+pitch intervals');     

                    % adding progression lines for three songs, MUSIC-LIKE
                            % notes
                            subplot(2,2,1); hold on;
                            amp = NoteFeatures_M{k,1}.ampscaled;      ptc = NoteFeatures_M{k,1}.ptcscaled;  
                            plot(amp,ptc,'ob','markersize',5);     hold on
                            newxy = fnplt(cscvn([amp';ptc']),'b',1); % plotting interpolated round curve through the acutal points
                            plt = plot(newxy(1,:),newxy(2,:),'b', 'LineWidth',1); hold on
                                %// modified blueish winter-colormap
                                n=size(newxy,2)+20; % this is used to make a good colormap of the correct size
                                cd = [uint8(winter(n)*255) uint8(ones(n,1))].'; %' for speech, use colormap autumn or hot
                                drawnow
                                set(plt.Edge, 'ColorBinding','interpolated', 'ColorData',cd(:,1:size(newxy,2)))
                            axis([0 1.25 0 1.25])

                            %intervals
                            subplot(2,2,3); hold on;    plot([-2 2],[0 0],'Color',[.5 .5 .5]);      plot([0 0],[-2 2],'Color',[.5 .5 .5]);
                            plot(diff(amp),diff(ptc),'ob','markersize',5);     hold on
                            newxy = fnplt(cscvn([diff(amp)';diff(ptc)']),'b',1); % plotting interpolated round curve through the acutal points
                            plt = plot(newxy(1,:),newxy(2,:),'b', 'LineWidth',1); hold on
                                %// using again the modified blueish winter-colormap generated above
                                drawnow
                                set(plt.Edge, 'ColorBinding','interpolated', 'ColorData',cd(:,1:size(newxy,2)))
                             axis([-.8 .8 -.8 .8])
   
                % SPEECH: plotting gray backdrop for Note features and Interval features (amp & pitch) - SPEECH-LIKE:

                    subplot(2,2,2); hold on;    %plot([-2 2],[0 0],'Color',[.5 .5 .5]);  plot([0 0],[-2 2],'Color',[.5 .5 .5]); 
                    plot(allNoteAmp_s,allNotePitch_s,'.','Color',[.75 .75 .75]);    xlabel('amplitude');    ylabel('pitch');    box on;    
                    title('Speech-like');    
                    
                    subplot(2,2,4); hold on;    plot([-2 2],[0 0],'Color',[.5 .5 .5]);  plot([0 0],[-2 2],'Color',[.5 .5 .5]); 
                    plot(allNoteAmpDiff_s,allNotePitchDiff_s,'.','Color',[.75 .75 .75]); xlabel('amplitude interval');  ylabel('pitch interval');   box on  
%                     title('Speech-like amp+pitch intervals');    

                    % adding progression lines for three songs, SPEECH-LIKE
                            % notes
                            subplot(2,2,2); hold on;
                            amp = NoteFeatures_s{k,1}.ampscaled;      ptc = NoteFeatures_s{k,1}.ptcscaled;  
                            plot(amp,ptc,'or','markersize',5);     hold on
                            newxy = fnplt(cscvn([amp';ptc']),'r',1); % plotting interpolated round curve through the acutal points
                            plt = plot(newxy(1,:),newxy(2,:),'r', 'LineWidth',1); hold on
                                %// modified autumn-colormap
                                n=size(newxy,2)+40; % this is used to make a good colormap of the correct size
                                cd = [uint8(autumn(n)*255) uint8(ones(n,1))].'; %' for speech, use colormap autumn or hot
                                drawnow
                                set(plt.Edge, 'ColorBinding','interpolated', 'ColorData',cd(:,1:size(newxy,2)))
                            axis([0 1.25 0 1.25])

                            %intervals
                            subplot(2,2,4); hold on; 
                            plot([-1 1],[0 0],'Color',[.5 .5 .5]);      plot([0 0],[-1 1],'Color',[.5 .5 .5]);
                            plot(diff(amp),diff(ptc),'or','markersize',5);     hold on
                            newxy = fnplt(cscvn([diff(amp)';diff(ptc)']),'r',1); % plotting interpolated round curve through the acutal points
                            plt = plot(newxy(1,:),newxy(2,:),'r', 'LineWidth',1); hold on
                                %// using again the modified autumn-colormap generated above
                                drawnow
                                set(plt.Edge, 'ColorBinding','interpolated', 'ColorData',cd(:,1:size(newxy,2)))                
                            axis([-.8 .8 -.8 .8])
                            
                            wavname_M = NoteFeatures_M{k,1}.filename;   vec = strsplit(wavname_M,'.');  wavname_M = vec{1};
                            wavname_s = NoteFeatures_s{k,1}.filename;   vec = strsplit(wavname_s,'.');  wavname_s = vec{1};
                            sgtitle([wavname_M,' & ',wavname_s])
            end
            
            
            
            
            
%% NOW RHYTHM: start with durations, ratios, etc.

    
    % onsets
        % Music-like
        figure(71+100*TESTRUNNUM);  clf;   set(gcf,'Position',[560 710 431 238])
        ctm1=1;     ctm2=1;
        allDur_M = [];        allRatio_M = [];      allCycledur_M = [];
        for i = 1:numel(NoteFeatures_M)
            
                onsets = NoteFeatures_M{i,1}.onsets;   
                durations = onsets(2:end)-onsets(1:end-1);
                cycledurations = durations(1:end-1)+durations(2:end);
                % calculating ratios
                ratios = durations(1:end-1)./cycledurations;
                allDur_M(ctm1:ctm1+numel(durations)-1) = durations;
                allRatio_M(ctm2:ctm2+numel(ratios)-1) = ratios;
                allCycledur_M(ctm2:ctm2+numel(ratios)-1) = cycledurations;
                ctm1=ctm1+numel(durations);
                ctm2=ctm2+numel(ratios);
        end
        IS_CDSLOPASS = false;   IS_CDFASTPASS = false;  cd_cutoff=1;	IS_OPENNEWFIG=false;    nbins=35;
        RatioHist = plotRatioHist_pickNBins(allRatio_M,allCycledur_M,nbins,'music-like',IS_CDSLOPASS,IS_CDFASTPASS,cd_cutoff,IS_OPENNEWFIG);           

    
        % speech-like
        figure(72+100*TESTRUNNUM);  clf;   set(gcf,'Position',[560 398 431 238])
        ctm1=1;     ctm2=1;
        allDur_s = [];        allRatio_s = [];      allCycledur_s = [];
        for i = 1:numel(NoteFeatures_s)
            
                onsets = NoteFeatures_s{i,1}.onsets;   
                durations = onsets(2:end)-onsets(1:end-1);
                cycledurations = durations(1:end-1)+durations(2:end);
                % calculating ratios
                ratios = durations(1:end-1)./cycledurations;
                allDur_s(ctm1:ctm1+numel(durations)-1) = durations;
                allRatio_s(ctm2:ctm2+numel(ratios)-1) = ratios;
                allCycledur_s(ctm2:ctm2+numel(ratios)-1) = cycledurations;
                ctm1=ctm1+numel(durations);
                ctm2=ctm2+numel(ratios);
        end
        IS_CDSLOPASS = false;   IS_CDFASTPASS = false;  cd_cutoff=1;	IS_OPENNEWFIG=false;    nbins=30;
        RatioHist = plotRatioHist_pickNBins(allRatio_s,allCycledur_s,nbins,'speech-like',IS_CDSLOPASS,IS_CDFASTPASS,cd_cutoff,IS_OPENNEWFIG);           

            
        %% now histograms for durations
        figure(73+100*TESTRUNNUM); clf
        
        % Music, blue pdf
        subplot(2,1,1)
        histcolor1 = [.75 .75 .75];      alphaval = 1;       nbins = 35;
        DurHist_M = histogram(allDur_M,nbins,'Normalization','pdf','FaceColor',histcolor1,'FaceAlpha',alphaval,'EdgeColor',histcolor1,'EdgeAlpha',alphaval); % 130 or 153 bins for string quartet?
        plotylim = ylim;
        hold on;
        % adding pdf 
            [f,xi] = ksdensity(allDur_M,'bandwidth',5); %'bandwidth',0.008); %,'support','positive');
            % plotting distribution 
            py = plot(xi,f,'Color','b','linewidth',2); % [.8 .8 .8]);
            xlim([0 800])
            ylabel('probability density')

        % speech, red pdf
        subplot(2,1,2)
        histcolor1 = [.75 .75 .75];      alphaval = 1;       nbins = 30;
        DurHist_s = histogram(allDur_s,nbins,'Normalization','pdf','FaceColor',histcolor1,'FaceAlpha',alphaval,'EdgeColor',histcolor1,'EdgeAlpha',alphaval); % 130 or 153 bins for string quartet?
        plotylim = ylim;
        hold on;
        % adding pdf 
            [f,xi] = ksdensity(allDur_s,'bandwidth',15); %'bandwidth',0.008); %,'support','positive');
            % plotting distribution 
            py = plot(xi,f,'Color','r','linewidth',2); 
            xlim([0 800])
            xlabel('duration [ms]');    ylabel('probability density')
            
        sgtitle('Onset-to-onset durations')
            
            
%% Now make a table for Pauline with all data she needs for stats. 

% writing data into a table Pauline can use

meanDur_M_then_s = nan(numel(NoteFeatures_M)*2,1);       SEMDur_M_then_s = nan(numel(NoteFeatures_M)*2,1);
meanAmp_M_then_s = nan(numel(NoteFeatures_M)*2,1);       SEMAmp_M_then_s = nan(numel(NoteFeatures_M)*2,1);
meanEnt_M_then_s = nan(numel(NoteFeatures_M)*2,1);       SEMEnt_M_then_s = nan(numel(NoteFeatures_M)*2,1);
meanPtc_M_then_s = nan(numel(NoteFeatures_M)*2,1);       SEMPtc_M_then_s = nan(numel(NoteFeatures_M)*2,1);

meanDurChange_M_then_s = nan(numel(NoteFeatures_M)*2,1);       SEMDurChange_M_then_s = nan(numel(NoteFeatures_M)*2,1);
meanAmpChange_M_then_s = nan(numel(NoteFeatures_M)*2,1);       SEMAmpChange_M_then_s = nan(numel(NoteFeatures_M)*2,1);
meanEntChange_M_then_s = nan(numel(NoteFeatures_M)*2,1);       SEMEntChange_M_then_s = nan(numel(NoteFeatures_M)*2,1);
meanPtcChange_M_then_s = nan(numel(NoteFeatures_M)*2,1);       SEMPtcChange_M_then_s = nan(numel(NoteFeatures_M)*2,1);

for i = 1:numel(NoteFeatures_M)

        onsets = NoteFeatures_M{i,1}.onsets;   
        durations = onsets(2:end)-onsets(1:end-1);
        cycledurations = durations(1:end-1)+durations(2:end);
        % calculating ratios
        ratios = durations(1:end-1)./cycledurations;
        meanDur_M_then_s(i) = mean(durations);                  SEMDur_M_then_s(i) = std(durations)/sqrt(length(durations));   
        
        amp = NoteFeatures_M{i,1}.amp;      meanAmp_M_then_s(i) = mean(amp);    SEMAmp_M_then_s(i) = std(amp)/sqrt(length(amp));
        ent = NoteFeatures_M{i,1}.ent;      meanEnt_M_then_s(i) = mean(ent);    SEMEnt_M_then_s(i) = std(ent)/sqrt(length(ent));
        ptc = NoteFeatures_M{i,1}.ptc;      meanPtc_M_then_s(i) = mean(ptc);    SEMPtc_M_then_s(i) = std(ptc)/sqrt(length(ptc));
        
        ampChange = diff(amp);      meanAmpChange_M_then_s(i) = mean(ampChange);    SEMAmpChange_M_then_s(i) = std(ampChange)/sqrt(length(ampChange));
        entChange = diff(ent);      meanEntChange_M_then_s(i) = mean(entChange);    SEMEntChange_M_then_s(i) = std(entChange)/sqrt(length(entChange));
        ptcChange = diff(ptc);      meanPtcChange_M_then_s(i) = mean(ptcChange);    SEMPtcChange_M_then_s(i) = std(ptcChange)/sqrt(length(ptcChange));


        
end

for i = 1:numel(NoteFeatures_s)
    
        onsets = NoteFeatures_s{i,1}.onsets;   
        durations = onsets(2:end)-onsets(1:end-1);
        cycledurations = durations(1:end-1)+durations(2:end);
        % calculating ratios
        ratios = durations(1:end-1)./cycledurations;
        meanDur_M_then_s(i+numel(NoteFeatures_M)) = mean(durations);        SEMDur_M_then_s(i+numel(NoteFeatures_M)) = std(durations)/sqrt(length(durations));   

        amp = NoteFeatures_s{i,1}.amp;      meanAmp_M_then_s(i+numel(NoteFeatures_M)) = mean(amp);    SEMAmp_M_then_s(i+numel(NoteFeatures_M)) = std(amp)/sqrt(length(amp));
        ent = NoteFeatures_s{i,1}.ent;      meanEnt_M_then_s(i+numel(NoteFeatures_M)) = mean(ent);    SEMEnt_M_then_s(i+numel(NoteFeatures_M)) = std(ent)/sqrt(length(ent));
        ptc = NoteFeatures_s{i,1}.ptc;      meanPtc_M_then_s(i+numel(NoteFeatures_M)) = mean(ptc);    SEMPtc_M_then_s(i+numel(NoteFeatures_M)) = std(ptc)/sqrt(length(ptc));

        ampChange = diff(amp);      meanAmpChange_M_then_s(i+numel(NoteFeatures_M)) = mean(ampChange);    SEMAmpChange_M_then_s(i+numel(NoteFeatures_M)) = std(ampChange)/sqrt(length(ampChange));
        entChange = diff(ent);      meanEntChange_M_then_s(i+numel(NoteFeatures_M)) = mean(entChange);    SEMEntChange_M_then_s(i+numel(NoteFeatures_M)) = std(entChange)/sqrt(length(entChange));
        ptcChange = diff(ptc);      meanPtcChange_M_then_s(i+numel(NoteFeatures_M)) = mean(ptcChange);    SEMPtcChange_M_then_s(i+numel(NoteFeatures_M)) = std(ptcChange)/sqrt(length(ptcChange));

    
end

    fieldnames4Pauline = {'Filename','Mean_NoteRate_Hz',...
                            'Mean_onset2onsetDuration_ms','SEM_onset2onsetDuration'...
                            'Mean_Amplitude','SEM_Amplitude',...
                            'Mean_Entropy','SEM_Entropy',...
                            'Mean_Pitch','SEM_Pitch',...
                            'Mean_AmpChange','SEM_AmpChange',...
                            'Mean_EntChange','SEM_EntChange',...
                            'Mean_PitchChange','SEM_PitchChange'};

    T = table(filenames_M_then_s,noterate_M_then_s,meanDur_M_then_s,SEMDur_M_then_s,...
        meanAmp_M_then_s,SEMAmp_M_then_s,...
        meanEnt_M_then_s,SEMEnt_M_then_s,...
        meanPtc_M_then_s,SEMPtc_M_then_s,...
        meanAmpChange_M_then_s,SEMAmpChange_M_then_s,...
        meanEntChange_M_then_s,SEMEntChange_M_then_s,...
        meanPtcChange_M_then_s,SEMPtcChange_M_then_s,...
        'VariableNames',fieldnames4Pauline);

    % Writing table to excel sheet for Pauline

    cd(results_dir);

    % transfer to xls

    writetable(T,'Dundun_Data4Pauline.xls','WriteRowNames',true);  


            
%% \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
%  \\\\\\ -------------------------------------------------------------------------------------------------------------------------------- \\\\\\\\\
%  \\\\\\ -------------------------------------------------------------------------------------------------------------------------------- \\\\\\\\\
%  \\\\\\ -------------------------------------------------------------------------------------------------------------------------------- \\\\\\\\\            
% \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\        
%%





% %% OLD    
%     
% %% adding global time to structure "Rhythm" - do this if Rhythm is from one performance of one individual.
% % Rhythm = addGlobalTime2Rhythm(Rhythm);
% 
% 
%       %% going through entire corpus' amp/freq/ent to determine running 1st & 99.5th percentiles... 
%       fprintf('preparing standardization of amplitude, entropy, and frequency vectors...\n') 
%       
%           allrawamp = []; allrawfreq = []; allrawent = [];
%           ixs_corpus = zeros(length(corpus),2); % start and end indexes of each corpus entry
%           dur_corpus = zeros(length(corpus),1); % duration of each corpus entry
%           ctc=1;
%           for c = 1:length(corpus)
%               thisrawamp=corpus{c,1}.rawamplitude;
%               thisrawfreq=corpus{c,1}.rawfrequency;
%               thisrawent=corpus{c,1}.rawentropy;
%               allrawamp(ctc:ctc+numel(thisrawamp)-1) = thisrawamp;
%               allrawfreq(ctc:ctc+numel(thisrawfreq)-1) = thisrawfreq;
%               allrawent(ctc:ctc+numel(thisrawent)-1) = thisrawent;
%               
%               ixs_corpus(c,1) = ctc; ixs_corpus(c,2) = ctc+numel(thisrawamp)-1; % this variable just gives you start and end indexes for each corpus entry
%               dur_corpus(c,1) = numel(thisrawamp); 
%               
%               ctc = ctc+numel(thisrawamp);
%           end
%           
%        if IS_LOCALNORM
%           tic
%           
%           normperc_amp = nan(length(corpus),2);        
%           normperc_freq = nan(length(corpus),2);       
%           normperc_ent = nan(length(corpus),2); 
%           % calculating "running" percentiles
%           windowsize = 4000;
%           for c = 1:length(corpus)
%               % which corpus entries to include to calculate percentiles for standardizing?
%               cx_4norm = c;
%               while sum(dur_corpus(cx_4norm))<windowsize
%                   % as long as windowsize isn't reached, include more corpus entries (left and right) into normalization procedure 
%                   cx_4norm = sort(unique([cx_4norm,max(cx_4norm(1)-1,1),min(cx_4norm(end)+1,length(corpus))]));
%               end
% %               test(c).cx_4norm = cx_4norm;
%               ixs4norm_thisc = [ixs_corpus(cx_4norm(1),1) : ixs_corpus(cx_4norm(end),2)];
%               normperc_amp(c,1) = prctile(allrawamp(ixs4norm_thisc),percLo);     normperc_amp(c,2) = prctile(allrawamp(ixs4norm_thisc),percHi);
%               normperc_freq(c,1) = prctile(allrawfreq(ixs4norm_thisc),percLo);   normperc_freq(c,2) = prctile(allrawfreq(ixs4norm_thisc),percHi);
%               normperc_ent(c,1) = prctile(allrawent(ixs4norm_thisc),percLo);     normperc_ent(c,2) = prctile(allrawent(ixs4norm_thisc),percHi);              
% 
%           end
%           
%           toc
%       else % no local, but global normalization
%           
%           globnormperc_amp(1,1) = prctile(allrawamp,percLo);             globnormperc_amp(1,2) = prctile(allrawamp,percHi);
%           globnormperc_ent(1,1) = prctile(allrawent,percLo);             globnormperc_ent(1,2) = prctile(allrawent,percHi);
%           globnormperc_freq(1,1) = prctile(allrawfreq,percLo);           globnormperc_freq(1,2) = prctile(allrawfreq,percHi);
%           
%       end
%       
% %%
% fprintf('standardizing amplitude, entropy, and frequency, then segmenting...\n')
% 
% for i = 1:size(corpus,1)
%     
%     %% for scaling/normalizing, either use a local window (IS_LOCALNORM = true), or global data, and scale between percentiles given above.
%     % Scaling on the basis of individual phrases is not a good idea, some of them may be too short!
%     
%     ampraw = corpus{i,1}.rawamplitude;          
%     freqraw = corpus{i,1}.rawfrequency;         
%     entraw = corpus{i,1}.rawentropy;            
%         
%         % standardization/normalization
%         if IS_LOCALNORM
%             % standardize using the local standardization parameters
%             ampscaled = (ampraw-(normperc_amp(i,1))) / (normperc_amp(i,2)-(normperc_amp(i,1))); % local scaling
%             freqscaled = (freqraw-(normperc_freq(i,1))) / (normperc_freq(i,2)-(normperc_freq(i,1))); % local scaling
%             entscaled = (entraw-(normperc_ent(i,1))) / (normperc_ent(i,2)-(normperc_ent(i,1))); % local scaling
%         else % global normalization
%             ampscaled = (ampraw-(globnormperc_amp(1,1))) / (globnormperc_amp(1,2)-(globnormperc_amp(1,1))); % general, global scaling
%             %ampscaled = (ampraw-(-14)) / (65-(-14)); % good for Ari, global scaling
%             freqscaled = (freqraw-(globnormperc_freq(1,1))) / (globnormperc_freq(1,2)-(globnormperc_freq(1,1)));% general, global scaling
%             %freqscaled = (freqraw-1080) / (6950-1080); % good for Ari, global scaling
%             entscaled = (entraw-(globnormperc_ent(1,1))) / (globnormperc_ent(1,2)-(globnormperc_ent(1,1)));% general, global scaling
%             %entscaled = (entraw-(-7)) / ((-.38)-(-7)); % good for Ari, global scaling
%         end
%                                                 
%                     if IS_FILTERED_FEATURES
%                         freqscaled = hpfilter(freqscaled,40);
%                         entscaled = hpfilter(entscaled,50);
%                     else
%                     end
%                                                 
% 
%     % "stupid segmentation"
%     % These thresholds were determined empirically using the "Rhythm" variable that was created on the first ten 20s-snippets of the AriMS-Song, 
%     % and on amp & ent scaled for range between 1st and 99th percentile, in those snippets.
%     
%     featurevec1 = ampscaled;
%     featurevec2 = entscaled;
%     threshold1 = .5; % threshold for scaled amplitude in AriMS
%     threshold2 = .82; % threshold for scaled entropy in AriMS. was .72 but there are tch-sounds that are higher entropy!
%     t1_higher1lower0 = 1; % for features like amplitude, sound is HIGHER than threshold
%     t2_higher1lower0 = 0; % for entropy, sound is LOWER than threshold
%     silence_mindur = 70; % 70
%     note_mindur = 20; % tentative minimal duration for sounds
%     
%     [note_starts,note_ends,unit_starts,unit_ends,silence_starts,silence_ends,note_ix,silence_ix] = segmentByTwoFeatureThresholds...
%         (featurevec1,featurevec2,threshold1,t1_higher1lower0,threshold2,t2_higher1lower0,silence_mindur,note_mindur);
%     
% %         % outcomment when all is good
% %             figure(100); clf
% %             set(gcf,'Position',[12 282 1650 420]); hold on
% %             plot(ampscaled,'linewidth',1.5)
% %             plot(note_ix,ones(size(note_ix)).*.7,'.k'); 
% %             plot(silence_ix,.75*ones(size(silence_ix)),'.r'); 
% %             ylim([-.1 1.3])
% 
% 
%     % now save stuff in the structure "units"
%     
%                 buffer_ms = 4;               
%                 
%                 for j = 1:numel(unit_starts)
%                     unit{i,j}.filename = char(corpus{i,1}.wavname);
%                     unit{i,j}.buffersize_ms = buffer_ms;
%                     bufferedstartix = max(unit_starts(j)-buffer_ms,1);
%                     bufferedendix = min(unit_ends(j)+buffer_ms,numel(ampscaled));
%                     if ~isempty(unit_starts)
%                         unit{i,j}.startix = bufferedstartix;
%                         unit{i,j}.endix = bufferedendix;
%                         unit{i,j}.duration = unit_ends(j)-unit_starts(j);
%                         unit{i,j}.syllablestarts = intersect(note_starts,[unit_starts(j):unit_ends(j)]);
%                         unit{i,j}.syllableends = intersect(note_ends,[unit_starts(j):unit_ends(j)]);
%                     else
%                         unit{i,j}.startix = [];
%                         unit{i,j}.endix = [];
%                         unit{i,j}.duration = [];
%                         unit{i,j}.syllablestarts = [];
%                         unit{i,j}.syllableends = [];                    
%                     end
%                     unit{i,j}.amp_scaled = ampscaled(bufferedstartix:bufferedendix);
%                     unit{i,j}.ent_scaled = entscaled(bufferedstartix:bufferedendix);
%                     unit{i,j}.freq_scaled = freqscaled(bufferedstartix:bufferedendix);
% 
% 
%                 end
% 
% 
%     
% end
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %--------------------------------------////// JUST LOOKING AT TRANSITIONS //////-------------------------------------%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Now go through all transitions between phrases and compare first and last units. Plot until you know it's working!
% 
% 
% %% make sure between the pre-transition unit and the post-transition unit, there is no longer break than a certain threshold (e.g. 500 ms).
% % WEITERMACHEN VORSCHLAG: Wenn der duration Unterschied zwischen pre- und post-transition syllable sehr lang ist, checken ob er geringer
% % wäre zwischen pre und post-transition SYLLABLES?. 
% 
% thresh_transtime = 3000; % threshold for silence between the two units around the transition.
% 
% transition_ok1_notok0 = nan(size(unit,1),1);
% 
%     u_exists = nan(size(unit,1),size(unit,2)); % this table tells you whether there are units or not (1 means it's empty)
%         for x = 1:size(unit,1)
%             for y = 1:size(unit,2)
%                 if ~isempty(unit{x,y})
%                     u_exists(x,y) = 1;
%                 end
%             end
%         end
% 
% u_beforetrans_ix = nan(size(unit,1),1);
%         
% for i = 1:size(unit,1)-1 % number of phrases -1 is number of transitions
%     
%     % determining the silence after unit 1, before the transition boundary.
%     % First, which unit is the last one before the transition?
%     tu1_unitix = nansum(u_exists(i,:)); 
%     u_beforetrans_ix(i) = tu1_unitix;
%     
%     if tu1_unitix>0 % there must be SOME unit in there for the transition silence length to be calculated
%         tu1_endix = unit{i,tu1_unitix}.endix;
%         tu1_filelength = length(corpus{i,1}.rawamplitude);
%         tu1_pretransgap = tu1_filelength-tu1_endix;
%     elseif tu1_unitix==0  % there is no unit in this phrase. don't use this transition.
%         transition_ok1_notok0(i) = 0;
%     else
%     end
%     
%     % determining the silence before unit 2, after the transition boundary.
%     tu2_exist = nansum(u_exists(i+1,:)); 
%     
%     if tu2_exist && tu1_unitix>0
%         tu2_startix = unit{i+1,1}.startix; % transition unit 2 (after the transition)
%         if (tu2_startix+tu1_pretransgap)<=thresh_transtime
%             transition_ok1_notok0(i) = 1;
%         elseif (tu2_startix+tu1_pretransgap)>thresh_transtime
%             transition_ok1_notok0(i) = 0;
%         else
%         end
%     else
%         transition_ok1_notok0(i) = 0;
%     end
%     
%     
% end
% transition_ok1_notok0(size(unit,1),1) = 0;
