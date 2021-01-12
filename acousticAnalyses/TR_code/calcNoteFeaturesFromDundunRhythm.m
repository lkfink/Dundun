% calcNoteFeaturesFromDundunRhythm. This is only for one line of rhythm!!

function [NoteFeatureStruct] = calcNoteFeaturesFromDundunRhythm(RhythmStruct,lineInRhythmStruct,IS_NOTELENGTH30,USE_NOTE_PEAK,IS_FILTERED_FEATURES,IS_PLOT) 

%     % the NOTELENGTH30 (30ms from note onset) seem to be best for feature calculation. Checked pitch and amplitude, both seem good.
%     IS_NOTELENGTH30 = true; % if you want to use a standard note length of 30 ms for mean feature calculation, from syllable onset on
%     USE_NOTE_PEAK = false; % if you want to use just 19ms around the peak for note intensity calculation
%     RhythmStruct = Dundun_M;
%     lineInRhythmStruct = 10;
%     IS_FILTERED_FEATURES = false;
%     IS_PLOT=true;

% NoteFeatureStruct = cell(1,1);

data_dir = '~/Dropbox/MATLAB/Dundun/data/wavs';
rdir=pwd;


    tic
    i = lineInRhythmStruct;
    
    if ~isfield(RhythmStruct,'rawpitch')
        [RhythmStruct] = addPitch2RhythmStruct(RhythmStruct,data_dir);        
    end

        
        filename = [RhythmStruct(i).wavname];
        amp = [RhythmStruct(i).rawamplitude];     ampsc = (amp-(prctile(amp,2))) / (prctile(amp,98)-(prctile(amp,2)));
        ent = [RhythmStruct(i).rawentropy];       entsc = (ent-(prctile(ent,2))) / (prctile(ent,98)-(prctile(ent,2)));
        freq = [RhythmStruct(i).rawfrequency];    freqsc = (freq-(prctile(freq,2))) / (prctile(freq,98)-(prctile(freq,2)));
        ptc = [RhythmStruct(i).rawpitch];    ptcsc = (ptc-(prctile(ptc,2))) / (prctile(ptc,98)-(prctile(ptc,2)));
        
        % calculate pitch if necessary:
%         if ~isfield(RhythmStruct,'rawpitch')
%                         
%             t_SAP = (0:length(mfreq_SAP)-1)/1000;
%             
%             cd(data_dir)
%             [x,fs] = audioread(filename);
%             
%             [f0mtlb,locs] = pitch(x,fs); % TURNS OUT THIS LOOKS PRETTY GOOD. COMPARABLE TO YIN_ESTIMATOR, AND MEANINGFUL. except for length of outputs which is weird. 100Hz only?!.
%             %t = (0:length(x)-1)/fs;
%             t0mtlb = (locs-1)/fs;
% 
%             % how to make this ms wise?? just resample, stupid workaround but whatever.
%             P=size(t_SAP,2);
%             Q=size(t0mtlb,1);
% 
%             pitchvect = resample(f0mtlb,P,Q);
%             
%             cd(rdir)
%         end
        
               
        if IS_FILTERED_FEATURES
            amp = hpfilter(amp,50);         ampsc = hpfilter(ampsc,50);
            ent = hpfilter(ent,50);         entsc = hpfilter(entsc,50);
            freq = hpfilter(freq,50);       freqsc = hpfilter(freqsc,100);
            ptc = hpfilter(ptc,10);         ptcsc = hpfilter(ptcsc,10);
        end
        
            if IS_NOTELENGTH30
                onsets = [RhythmStruct(i).syllablestarts];
                offsets = onsets+32;
            elseif USE_NOTE_PEAK
                onsets = [RhythmStruct(i).peaklocs_allpeaks]-3;
                offsets = [RhythmStruct(i).peaklocs_allpeaks]+15;
            else
                onsets = [RhythmStruct(i).syllablestarts];
                offsets = [RhythmStruct(i).syllableends];
            end

       NoteAmp = nan(numel(onsets),1);           NoteEnt = nan(numel(onsets),1);           NoteFreq = nan(numel(onsets),1);           NotePitch = nan(numel(onsets),1);     
       NoteAmpScaled = nan(numel(onsets),1);     NoteEntScaled = nan(numel(onsets),1);     NoteFreqScaled = nan(numel(onsets),1);     NotePitchScaled = nan(numel(onsets),1);     
       Note_o2o_dur = nan(numel(onsets),1);   Note_ratio = nan(numel(onsets),1);   
       
       for j = 1:numel(onsets)   
           
                NoteAmpScaled(j) = mean(ampsc(onsets(j):offsets(j))); % Play around! Also check whether assuming 100ms notes gives "cleaner" results
                NoteEntScaled(j) = mean(entsc(onsets(j):offsets(j)));
                NoteFreqScaled(j) = mean(freqsc(onsets(j):offsets(j)));
                NotePitchScaled(j) = mean(ptcsc(onsets(j):offsets(j)));
            
                NoteAmp(j) = mean(amp(onsets(j):offsets(j))); % Play around! Also check whether assuming 100ms notes gives "cleaner" results
                NoteEnt(j) = mean(ent(onsets(j):offsets(j)));
                NoteFreq(j) = mean(freq(onsets(j):offsets(j)));
                NotePitch(j) = mean(ptc(onsets(j):offsets(j)));
                
                % Rhythm
                if j<numel(onsets) % for the last onset, you cannot calculate o2o_dur.
                    o2o_dur = onsets(j+1)-onsets(j);
                    Note_o2o_dur(j) = o2o_dur;
                else 
                    Note_o2o_dur(j) = nan;
                end
                if j<numel(onsets)-1 % for the last two onsets, you cannot calculate ratio.
                    o2o_nextdur = onsets(j+2)-onsets(j+1);
                    Note_ratio(j) = o2o_dur / (o2o_dur+o2o_nextdur);
                else 
                    Note_ratio(j) = nan;
                end
                                    
        end
%         NoteFeatureStruct{i,1}.filename = filename;  NoteFeatureStruct{i,1}.onsets = onsets;      NoteFeatureStruct{i,1}.offsets = offsets;  
%         NoteFeatureStruct{i,1}.amp = NoteAmp;        NoteFeatureStruct{i,1}.ent = NoteEnt;        NoteFeatureStruct{i,1}.freq = NoteFreq;  
%         NoteFeatureStruct{i,1}.ptc = NotePitch;  
        
        NoteFeatureStruct(1).filename = filename;       NoteFeatureStruct(1).onsets = onsets;      NoteFeatureStruct(1).offsets = offsets;  
        NoteFeatureStruct(1).amp = NoteAmp;             NoteFeatureStruct(1).ent = NoteEnt;        NoteFeatureStruct(1).freq = NoteFreq;  
        NoteFeatureStruct(1).ptc = NotePitch;  
        NoteFeatureStruct(1).ampscaled = NoteAmpScaled;        NoteFeatureStruct(1).entscaled = NoteEntScaled;        NoteFeatureStruct(1).freqscaled = NoteFreqScaled;  
        NoteFeatureStruct(1).ptcscaled = NotePitchScaled;  
        NoteFeatureStruct(1).o2o_dur = Note_o2o_dur;    NoteFeatureStruct(1).ratio = Note_ratio;    
    
    
   %% 
%        onsets = [NoteFeaturesStruct{k,1}.onsets];
%             offsets = [NoteFeatures_M{k,1}.offsets];
%             amp = [NoteFeatures_M{k,1}.amp];    freq = [NoteFeatures_M{k,1}.freq];    ent = [NoteFeatures_M{k,1}.ent];
            Note_ix_vect = [];
            AmpSc_vect = [];
            FreqSc_vect = [];
            EntSc_vect = [];
            PitchSc_vect = [];
            ct=1;
            for l = 1:numel(onsets)
                nx = [onsets(l):offsets(l)];
                Note_ix_vect(ct:ct+numel(nx)-1) = nx;
                AmpSc_vect(ct:ct+numel(nx)-1) = ones(size(nx)).*NoteAmpScaled(l); % amp
                EntSc_vect(ct:ct+numel(nx)-1) = ones(size(nx)).*NoteEntScaled(l); % ent
                FreqSc_vect(ct:ct+numel(nx)-1) = ones(size(nx)).*NoteFreqScaled(l); % freq
                PitchSc_vect(ct:ct+numel(nx)-1) = ones(size(nx)).*NotePitchScaled(l); % pitch

                ct=ct+numel(nx);
            end
            
        NoteFeatureStruct(1).Note_ix_vect = Note_ix_vect;
        NoteFeatureStruct(1).AmpScaled_vect = AmpSc_vect;
        NoteFeatureStruct(1).FreqScaled_vect = FreqSc_vect;
        NoteFeatureStruct(1).EntScaled_vect = EntSc_vect;
        NoteFeatureStruct(1).PitchScaled_vect = PitchSc_vect;

            
    toc


    %% Plotting
    
    if IS_PLOT
        figure(997);  clf;   set(gcf,'Position',[3 463 1680 456]);     hold on; 
        
        p1=plot(ampsc,'b'); p1.Color(4)=.25;
        pf=plot(Note_ix_vect,hpfilter(FreqSc_vect,800),'.r');         pp=plot(Note_ix_vect,PitchSc_vect,'.m');         pa=plot(Note_ix_vect,AmpSc_vect,'.b');
        
        ylim([-.1 1.4]);   
        legend([pf,pp,pa],{'frequency','pitch','amplitude'});   legend boxoff
        title(filename,'Interpreter','none')
        
    end
    
    