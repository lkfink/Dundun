% Quick script to read in dundun stims and calculate other timing measures
% I am interested in 
% LF - 20201204

% add path to mir toolbox
addpath(genpath('/Users/lauren.fink/Documents/MATLAB_utils/MIRtoolbox1.7.2'));

% initialize variables
tf = table;
stimdir = dir(fullfile('/Users/lauren.fink/Documents/Projects/dundun/Stimuli/Dundun-drum-recordings/','*.mp3'));
nfiles = length(stimdir);
nr = 1;
desired_Fs = 100;

% read in all files, compute measures of interest
for ifile = 1:nfiles
    
    % load in audio file
    fname = fullfile(stimdir(ifile).folder, stimdir(ifile).name);
    newStr = split(stimdir(ifile).name,'.');
    stim = newStr{1}; 
    fprintf('\nReading audio file: %s\n', fname)
    [y,fs] = audioread(fname);
    
    % calculate metrical centroid and strength, and pulse clarity
%     [mc,ms] = mirmetroid(fname, 'Frame', .5, 's');
%     msd = mirgetdata(ms);
%     mcd = mirgetdata(mc);
    [r, ac] = mirpulseclarity(fname);%, 'Frame', 5, 's', 10, '%');
    rd = mirgetdata(r);
    rac = mirgetdata(ac);
    
    % calculate amplitude envelope
    fprintf('\nExtracting envelope for audio file: %s\n', fname)
    N = fs*.05; % 50ms - 20 hz. Filter length should be function of lowest freq want to capture. 0.05 - 1 period of stim
    [YUPPER,YLOWER] = envelope(y, N, 'rms');
    
    % lowpass filter envelope at 10Hz or 20Hz
    fprintf('\nLow pass filtering envelope for audio file: %s\n', fname)
    [b,a] = butter(3, 50/(fs/2), 'low');
    filtaud = filtfilt(b,a,YUPPER(:,1));
    
    % downsample amplitude envelope
    fprintf('\nDownsampling envelope for audio file: %s\n', fname)
    [p,q] = rat(desired_Fs / fs); % desired fs/ original fs
    yu = resample(filtaud, p, q);
    
    % compute autocovariance
    maxlags = round(numel(yu)*0.5); % 50% of signal
    [xc,lag] = xcov(yu,maxlags, 'normalized');
    [~,df] = findpeaks(xc,'MinPeakDistance', desired_Fs);
    
%     figure
%     plot(lag/desired_Fs,xc,'k',...
%         lag(df)/desired_Fs, xc(df),'kv','MarkerFaceColor','r')
%     grid on
%     %xlim([-15 15])
%     xlabel('Time (secs)')
%     title('Auto-covariance')   
    

%     
    % Save all vars to table
    fprintf('\nSaving all data for file: %s\n', fname)
    pat = '\d{1,2}';
    pat2 = '\w';

    stimnum = regexp(stim, pat, 'match');
    stimclass = regexp(stim, pat2, 'match');
    tf.stim{nr,1} = stim;
    tf.stimnum(nr,1) = str2double(stimnum);
    tf.class{nr,1} = stimclass{end};
    tf.raw_audio{nr,1} = y;
    tf.ampEnv_downsampled_u{nr, 1} = yu;
    tf.xc{nr,1} = xc; 
    tf.lags{nr,1} = lag;
    tf.peaks{nr,1} = df;
%     tf.metrCent{nr,1} = mcd;
%     tf.metrStrength{nr,1} = msd;
    tf.pulseClarity(nr,1) = rd;
    tf.mirautocorr{nr,1} = rac;
%     tf.pulseClarity_mean(nr,1) = mean(rd);
%     tf.pulseClarity_std(nr,1) = std(rd);
    
    
    nr = nr+1;
end


tf = sortrows(tf,'stimnum','ascend');
tf = sortrows(tf,'class','ascend');
%% Make a few quick plots 
musmask = strcmp(tf.class, 'M');

% plot number of peaks in music vs speech
for irow = 1:length(tf.peaks)
    tf.numPeaks(irow) = numel(tf.peaks{irow});
end
figure()
histogram(tf.numPeaks(musmask))
hold on 
histogram(tf.numPeaks(~musmask))

% plot all acfs for each file, by category
figure()
for i = 1:length(tf.stim)
    subplot(6,5,i)
    plot(tf.lags{i}/desired_Fs,tf.xc{i},'k')%,...
         %tf.lags{i}(tf.peaks{i})/desired_Fs, tf.xc{i}(tf.peaks{i}),'kv','MarkerFaceColor','r')
         % TODO check something seems still off with peak plotting location
    ylim([-0.5, 1])
    xlim([-5, 5])
    grid on
    title(tf.stim{i})
    
    % check if signal is stationary
    tf.stationary{i} = adftest(tf.ampEnv_downsampled_u{i}, 'lags',round(numel(tf.ampEnv_downsampled_u{i})/4),'model','AR');
end

[H,P,STAT,C,REG] = adftest(tf.ampEnv_downsampled_u{i});


%     figure
%     plot(lag/desired_Fs,xc,'k',...
%         lag(df)/desired_Fs, xc(df),'kv','MarkerFaceColor','r')
%     grid on
%     %xlim([-15 15])
%     xlabel('Time (secs)')
%     title('Auto-covariance')   

% plot pulse clarity for music vs. speech
figure()
subplot(2,1,1)
hist(tf.pulseClarity(musmask))
xlim([.6, 1])
ylabel('# stimuli')
xlabel('Pulse clarity')
title('Music')
hold on 
subplot(2,1,2)
hist(tf.pulseClarity(~musmask))
xlim([.6, 1])
title('Speech')
ylabel('# stimuli')
xlabel('Pulse clarity')
suptitle('Pulse clarity of music vs. speech-like dundun')

%%
% plot all mir acfs for each file, by category
figure()
for i = 1:length(tf.stim)
    subplot(6,5,i)
    plot(tf.mirautocorr{i}, 'k')
    %plot(tf.lags{i}/desired_Fs,tf.xc{i},'k')%,...
         %tf.lags{i}(tf.peaks{i})/desired_Fs, tf.xc{i}(tf.peaks{i}),'kv','MarkerFaceColor','r')
         % TODO check something seems still off with peak plotting location
    ylim([-0.2, .4])
    xlim([0, 150])
    
    grid on
    title(tf.stim{i})
    
    % check if signal is stationary
    %tf.stationary{i} = adftest(tf.ampEnv_downsampled_u{i}, 'lags',round(numel(tf.ampEnv_downsampled_u{i})/4),'model','AR');
end

%%
% Save table to .mat file
% fpath = params.paths.MSM_table_path;
t = tf(:,{'stim', 'stimnum', 'class', 'pulseClarity'});
outfname = '/Users/lauren.fink/Documents/Projects/dundun/Data/pulseClarity.csv';
writetable(t, outfname)
% fprintf('\nSaving mat file: %s\n', outfname)
% save(outfname,'ampEnv', '-v7.3');
% fprintf('%s saved', outfname)