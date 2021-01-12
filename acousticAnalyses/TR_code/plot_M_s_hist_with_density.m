function [] = plot_M_s_hist_with_density(data_M,data_s,nbins_M,nbins_s,xstring,bandwidth,FIGURENUMBER)



    figure(FIGURENUMBER);     clf;   set(gcf,'Position',[900 607 284 199]);   hold on;

%         histogram(data_M,nbins_M,'Normalization','pdf','FaceColor','b','FaceAlpha',.3,'EdgeColor','none');%histogram(allNoteEnt_M,36,'normalization','pdf'); hne_M.FaceColor = [0 0 1]; hne_M.FaceAlpha = 1; hne_M.LineStyle = 'none';
%         histogram(data_s,nbins_s,'Normalization','pdf','FaceColor','r','FaceAlpha',.3,'EdgeColor','none');%histcolor1,'EdgeAlpha',.8); ?% histogram(allNoteEnt_s,25,'normalization','pdf'); hne_s.EdgeColor = [1 0 0]; hne_s.FaceAlpha = 0; hne_s.LineWidth = 2.2; 
        
        % attempt with edges on bars:
        histogram(data_M,nbins_M,'Normalization','pdf','FaceColor','b','FaceAlpha',.2,'EdgeColor','b','EdgeAlpha',.3);%histogram(allNoteEnt_M,36,'normalization','pdf'); hne_M.FaceColor = [0 0 1]; hne_M.FaceAlpha = 1; hne_M.LineStyle = 'none';
        histogram(data_s,nbins_s,'Normalization','pdf','FaceColor','r','FaceAlpha',.2,'EdgeColor','r','EdgeAlpha',.3);%histcolor1,'EdgeAlpha',.8); ?% histogram(allNoteEnt_s,25,'normalization','pdf'); hne_s.EdgeColor = [1 0 0]; hne_s.FaceAlpha = 0; hne_s.LineWidth = 2.2; 
        xlabel(xstring); ylabel('probability density'); box on
        % density
        [f_M,xi_M] = ksdensity(data_M, 'bandwidth',bandwidth);%0.026); %'bandwidth',0.008); %,'support','positive'); 0.008 0.012
        [f_s,xi_s] = ksdensity(data_s, 'bandwidth',bandwidth);%0.026); %'bandwidth',0.008); %,'support','positive');
        plot(xi_M,f_M,'Color','b','linewidth',2.6); % [.8 .8 .8]);
        plot(xi_s,f_s,'Color','r','linewidth',2.6);
%         xlim([0 1])

end
