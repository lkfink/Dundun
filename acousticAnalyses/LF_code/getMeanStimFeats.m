% quick script to get Tina's data in useable format for modeling in R
% want to grab means of all features 

% load Tina's data
load('/Users/lauren.fink/Documents/Projects/dundun/Data/Dundun_data_for_Lauren.mat')
AMS = load('/Users/lauren.fink/Documents/Projects/dundun/Data/Dundun_AMS_for_Lauren.mat');

%-------------------------------------------------------------------------%
% Create new table with means for each feature
% NOTE: relies on the custom functions splitOnNaN.m and getMeanForFeat.m

nd = table;
nd.stimulus = {'1M';'2M';'3M';'4M';'5M';'6M';'7M';'8M';'9M';'10M';'11M';'12M';'13M';'14M';'15M';'1S';'2S';'3S';'4S';'5S';'6S';'7S';'8S';'9S';'10S';'11S';'12S';'13S';'14S';'15S'};
% NOTE: hard-coding stimulus names ^ because this is the order they are
% returned from the functions below

% Intensity
nd.amp = getMeansForFeat(all_N_amp_M, all_N_amp_s);
% switched to abs value diff 20201201
nd.ampDiff = getMeansForFeat(abs(all_N_amp_diff_M), abs(all_N_amp_diff_s));

% Pitch
nd.pitch = getMeansForFeat(all_N_ptc_M, all_N_ptc_s);
nd.pitchDiff = getMeansForFeat(abs(all_N_ptc_diff_M), abs(all_N_ptc_diff_s));

% Timbre
nd.ent = getMeansForFeat(all_N_ent_M, all_N_ent_s);
nd.entDiff = getMeansForFeat(abs(all_N_ent_diff_M), abs(all_N_ent_diff_s));

% O2O
nd.o2o = getMeansForFeat(abs(all_N_O2O_M), abs(all_N_O2O_s));
%nd.o2oDiff = getMeansForFeat(all_N_O2O_diff_M, all_N_O2O_diff_s);
% ^ ahh there is issue of speech vector containing NaNs within a stimulus,
% so splitting on NaN doesn't work properly and mean vector is too long
% Ahh maybe this doesn't matter because we ended up going with the ratio
% measure, rather than the O2O change in the paper figs

% ratio
nd.ratio = getMeansForFeat(all_N_ratio_M, all_N_ratio_s);
% TODO CHECK with Tina - why is there not all_N_amp_ratio_diff_M?
% Ahh it seems this could be the label of 'all_N_ratio_diff_M', right?
nd.ratioDiff = getMeansForFeat(abs(all_N_ratio_diff_M), abs(all_N_amp_ratio_diff_s));


% Amplitude modulation spectrum
%nd.AMS_peak = getAmsPeak(AMS.Spectra_M_smooth, AMS.Spectra_s_smooth, AMS.freq_on_xaxis);
nd.AMS_peak = getAmsPeak(AMS.Spectra_M, AMS.Spectra_s, AMS.freq_on_xaxis);
% NOTE: will get different results for peak freq depending on whether you
% use raw or smoothed spectrum. Good with smooth because that is what is
% plotted in papers (ours and Ding)

%-------------------------------------------------------------------------%
% Ouput data table to csv
writetable(nd, '/Users/lauren.fink/Documents/Projects/dundun/Data/dundun_acoustFeat_means_absDiffs.csv')


