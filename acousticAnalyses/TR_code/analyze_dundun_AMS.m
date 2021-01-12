% analyze_dundun.m

clear all; close all; clc

%% loading the corpus data that contains raw amplitude, frequency, entropy for all phrases.

cache_dir='~/processed'; 
data_dir='~/Dropbox/MATLAB/Dundun/data/wavs';
results_dir='~/Dropbox/MATLAB/Dundun/results';
savefile = 'Dundun_AMS.mat';

    
    cd(data_dir)
    
    figure(444); clf
    
    %% speech-like
    
    adir=dir('*s.wav');
    Spectra_s = nan(numel(adir),600);
    Spectra_s_smooth = nan(numel(adir),600);
    
    for I = 1:numel(adir)
        
        filename = char(adir(I).name);
            Filenames_M{I,1} = filename;
            fprintf('now in file %d of %d... (%s)\n',I,numel(adir),filename);
            generateAMS_script; % relevant output (spectrum) is variable ms_rms. the according frequencies are f(1:size(ms,1)), or "freq_on_xaxis".            
            Spectra_s(I,:) = ms_rms;
            Spectra_s_smooth(I,:) = smooth(ms_rms,10);
            
            figure(444); hold on
            p=plot(freq_on_xaxis,smooth(ms_rms,10),'r','linewidth',.2);  p.Color(4)=.2;
            clear ms_rms ms ind
            
    end
     
    %% music-like
    
    adir=dir('*M.wav');
    Spectra_M = nan(numel(adir),600);
    Spectra_M_smooth = nan(numel(adir),600);

    
    for I = 1:numel(adir)
        
        filename = char(adir(I).name);
            Filenames_s{I,1} = filename;
            fprintf('now in file %d of %d... (%s)\n',I,numel(adir),filename);            
            generateAMS_script; % relevant output (spectrum) is variable ms_rms. the according frequencies are f(1:size(ms,1)), or "freq_on_xaxis".
            Spectra_M(I,:) = ms_rms;
            Spectra_M_smooth(I,:) = smooth(ms_rms,10);
            
            figure(444); hold on
            p=plot(freq_on_xaxis,smooth(ms_rms,10),'b','linewidth',.2); p.Color(4)=.2;
            clear ms_rms ms ind
            
    end       
    
   %% Plotting
   
   figure(446); 
   
   meanspect_M = mean(Spectra_M,1);
   meanspect_s = mean(Spectra_s,1);
    
%    plot(freq_on_xaxis,meanspect_M,'b','linewidth',2);  hold on
%    plot(freq_on_xaxis,meanspect_s,'r','linewidth',2);

   % smoothed
   plot(freq_on_xaxis,smooth(meanspect_M,10),'b','linewidth',2);  hold on
   plot(freq_on_xaxis,smooth(meanspect_s,10),'r','linewidth',2);

%    meanspect_M_smooth = mean(Spectra_M_smooth,1);
%    meanspect_s_smooth = mean(Spectra_s_smooth,1);
%    
%    plot(freq_on_xaxis,meanspect_M_smooth,'b','linewidth',2);  hold on
%    plot(freq_on_xaxis,meanspect_s_smooth,'r','linewidth',2);

    xlim([.25 32]);
    set(gca,'xscale','log')
    set(gca,'xtick',[.25 .5 1 2 4 8 16 32])
    xlabel('frequency (Hz)')
    ylabel('normalized amplitude')
    
    
        