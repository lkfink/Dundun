% Plots ratio histogram, with the option of only plotting a certain tempo

% INPUTS:
% IS_OPENNEWFIG ("true" or "false"): whether you want a new figure to open or plot in an existing figure.
% corpusname (some string): used for instance in plot title
% IS_CDSLOPASS: set to "true" if you want to plot only ratios of slow rhythms (cycle durations longer than what you set as cd_cutoff).
% IS_CDFASTPASS: set to "true" if you want to plot only ratios of fast rhythms (cycle durations shorter than what you set as cd_cutoff).
% cd_cutoff: only used if either IS_CDFASTPASS or IS_CDSLOPASS are set to "true". Cutoff in ms cycle duration.

% OUTPUTS:
% histogram plot of ratios
% RatioHist object (Matlab's standard histogram object)

% USE LIKE THIS:

    % IS_OPENNEWFIG = true;       corpusname = 'ZebraFinch';      IS_CDSLOPASS = false;
    % IS_CDFASTPASS = false;      cd_cutoff = 0;
    % RatioHist = plotRatioHist(ratios,cycledurations,corpusname,IS_CDSLOPASS,IS_CDFASTPASS,cd_cutoff,IS_OPENNEWFIG)

        % (To first calculate ratios and cycledurations from interval data, load "IntervalData_ZebraFinch_corpus.mat" and then do the following:)
        % cycledurations = ZebraFinch_intervals_ms(1:end-1)+ZebraFinch_intervals_ms(2:end);
        % ratios = ZebraFinch_intervals_ms(1:end-1)./(ZebraFinch_intervals_ms(1:end-1)+ZebraFinch_intervals_ms(2:end));

% IF YOU WANT TO PLOT JUST SLOW RHYTHMS, USE IS_CDSLOPASS (for FAST RHYTHMS, use IS_CDFASTPASS), and indicate the cycle duration cutoff 
% as cd_cutoff (in ms), like this:

%     IS_OPENNEWFIG = true;       corpusname = 'ZebraFinch';      IS_CDSLOPASS = true;
%     IS_CDFASTPASS = false;      cd_cutoff = 100;
%     RatioHist = plotRatioHist(ratios,cycledurations,corpusname,IS_CDSLOPASS,IS_CDFASTPASS,cd_cutoff,IS_OPENNEWFIG)
% 

function RatioHist = plotRatioHist_pickNBins(ratios,cycledurations,nbins,corpusname,IS_CDSLOPASS,IS_CDFASTPASS,cd_cutoff,IS_OPENNEWFIG)


if IS_OPENNEWFIG
    figure('name','RatioHistogram','Position',[360 480 369 218])
end

if IS_CDSLOPASS
    allratio_limitedcd = ratios(cycledurations>=cd_cutoff);
elseif IS_CDFASTPASS
    allratio_limitedcd = ratios(cycledurations<cd_cutoff);
else
    allratio_limitedcd = ratios;
end

histcolor = [.75 .75 .75]; %[.4 .4 .4]; 
alphaval = 1;
RatioHist = histogram(allratio_limitedcd,nbins,'Normalization','pdf','FaceColor',histcolor,'FaceAlpha',alphaval,'EdgeColor',histcolor,'EdgeAlpha',alphaval); % 130 or 153 bins for string quartet?
plotylim = ylim;
hold on;


% adding green dashed lines for the small integer ratios
    hold on;
    plot([1/2 1/2],[0 plotylim(2)],'--g','linewidth',1.5)
    plot([1/3 1/3],[0 plotylim(2)],'--g','linewidth',1.5)
    plot([1/4 1/4],[0 plotylim(2)],'--g','linewidth',1.5)
    % plot([1/5 1/5],[0 plotylim(2)],'--g','linewidth',1.5)
    plot([2/3 2/3],[0 plotylim(2)],'--g','linewidth',1.5)
    plot([3/4 3/4],[0 plotylim(2)],'--g','linewidth',1.5)
    % plot([4/5 4/5],[0 plotylim(2)],'--g','linewidth',1.5)
    xlim([0 1])


% adding pdf 
if ~isempty(allratio_limitedcd)
    [f,xi] = ksdensity(allratio_limitedcd, 'bandwidth',0.012); %'bandwidth',0.008); %,'support','positive');
        % plotting distribution 
        py = plot(xi,f,'Color','k','linewidth',2); % [.8 .8 .8]);
        xlim([0 1])
else
    disp('no ratio data, perhaps you excluded this tempo range?')
end

% adding title
if IS_CDFASTPASS 
    title([corpusname,' ratios pdf, cd<', num2str(cd_cutoff),' ms'],'Interpreter','none','FontSize',12)
elseif IS_CDSLOPASS            
    title([corpusname,' ratios pdf, cd>', num2str(cd_cutoff), ' ms'],'Interpreter','none','FontSize',12)
else
    title([corpusname,' ratios pdf'],'Interpreter','none')
end

xlabel('ratio')
ylabel('probability density')

