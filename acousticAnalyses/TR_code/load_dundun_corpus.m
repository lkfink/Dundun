% TEST_load_dundun_corpus
% the idea of this code is to create a "corpus" structure from dundun recordings. The corpus contains amp, freq, ent time series
% but no note onsets yet. 

function corpus=load_dundun_corpus(cache_dir,data_dir,MAX_READ_FILE,is_create)


% cache_dir='~/processed';
% data_dir='~/Dropbox/MATLAB/Dundun/data';
% MAX_READ_FILE=1000;
% is_create=true;
% 
% MSEC_limit=1/1000; % double check what Nori's MSEC_LIMIT was for...


% vec=data_dir;vec=strsplit(vec,'/');
% vec=vec{end};

cfname=sprintf('DUNDUN-FEATURES-MAX-%d.mat',MAX_READ_FILE); % for giving file in the cache (which will contain ms-wise data) a unique name
cd (cache_dir)
adir=dir(cfname); 

if isempty(adir)||is_create
    
%%    
    MAX_FILE=MAX_READ_FILE;
    cd (data_dir)
%     adir=dir('*.mat');
        adir=dir('*.wav');
    
    NF=min(MAX_FILE,length(adir)); % how many files will be read? (in case 1000 generates a too-big corpus file, draw this number down so you end up with different files. )
    corpus=cell(NF,1); % 
    tic
    
    % now read each file's amplitude, entropy, frequency
    for I=1:NF
        
            thisWavName = char(adir(I).name);
            % read in wav
            cd(data_dir);
            fprintf('now in file %d of %d... (%s)\n',I,NF,thisWavName);
            [trill,fs_original] = audioread(thisWavName);
            [~,~,~,m_entropy, m_amplitude, m_frequency,~,~,~,~]=extract_features(trill,44100); % extracting amplitude, ms-wise.
%             ampfilt = hpfilter(m_amplitude,50); % 50 filtering amplitude. This is good for thrush nightingales. seems good in mockingbirds too. No need to save it in "corpus", though, in which no segmentation is stored so far....
%             figure(10); clf; plot(m_amplitude); hold on; plot(ampfilt)
    
            corpus{I}.wavname = thisWavName;
            corpus{I}.rawamplitude = m_amplitude;
            corpus{I}.rawfrequency = m_frequency;
            corpus{I}.rawentropy = m_entropy;

            
    end
    toc
    
    cd (cache_dir);
    save(cfname,'corpus');
else
    fprintf('loading saved corpus from cache %s\n',cfname);
    cd (cache_dir);
    corpus=load(cfname);
    corpus=corpus.corpus;
    
end
    

%% NORI'S OLD CODE

%     mperm=randperm(NF);
%     for I=1:NF
%         %tic()
%         tfname=adir(mperm(I)).name;
% 
%         cnt=cnt+1;
%         cdata{cnt}=load(tfname);
%         
%         Notes=cdata{cnt}.data.Notes;
%         endtime=cdata{cnt}.data.endtime;
%         ofname=cdata{cnt}.data.fname;
%         tfname=tfname;
%         
%         fprintf('now in file %d of %d... (%s)\n',I,NF,ofname);
%         perc_ch=9;
%         pos=Notes(:,2)~=perc_ch;
%         Notes=Notes(pos,:);
%         %length(Notes)
%         
%         
%         %[sonors,places,durations,channels,max_poly]=sonor_me_2018(Notes, MSEC_limit,IS_PLOT);
%         [sonors,places,durations,durations_notes,max_poly]=sonor_me_rhythm(Notes,MSEC_limit);
%         
%         
%         corpus{cnt}.ofname=ofname;
%         corpus{cnt}.tfname=tfname;
%         corpus{cnt}.Notes=Notes;
%         corpus{cnt}.oNotes=cdata{cnt}.data.Notes;
%         corpus{cnt}.sonors=sonors;
%         corpus{cnt}.places=places;
%         corpus{cnt}.durations=durations;
%         corpus{cnt}.durations_notes=durations_notes;
%         corpus{cnt}.max_poly=max_poly;
%         
%         num_poly=nan(size(sonors));
%         for ll=1:length(sonors)
%             num_poly(ll)=length(sonors{ll});
%             assert(num_poly(ll)>=1)
%         end
%         corpus{cnt}.num_poly=num_poly;
%         LIMIT_POLY = 0; % added this, you can probably delete all poly stuff instead to clean up
%         pos=num_poly>=LIMIT_POLY;
%         
%         corpus{cnt}.sonors=corpus{cnt}.sonors(pos);
%         corpus{cnt}.places=corpus{cnt}.places(pos);
%         corpus{cnt}.durations=corpus{cnt}.durations(pos);
%         corpus{cnt}.durations_notes=corpus{cnt}.durations_notes(pos);
%         corpus{cnt}.num_poly=corpus{cnt}.num_poly(pos);
%         corpus{cnt}.cdata=cdata{cnt}.data;
%         %toc
%         
% %                     sonors: {1×691 cell}
% %              places: [1×691 double]
% %           durations: [1×691 double]
% %     durations_notes: [1×691 double]
% %     
%         
%     end
%     toc
%     cd (cache_dir);
%     save(cfname,'corpus');
% else
%     fprintf('loading saved corpus from cache %s\n',cfname);
%     cd (cache_dir);
%     corpus=load(cfname);
%     corpus=corpus.corpus;
%     
% end