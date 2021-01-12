function [] = plot_M_s_scatter_means_STE(NoteFeatures_struct_M,NoteFeatures_struct_s,fieldname_string,ylabelstring,FIGURENUMBER,IS_SCALED)


    

    figure(FIGURENUMBER);  clf;   set(gcf,'Position',[108 282 289 392])
    meanfeat_M = nan(numel(NoteFeatures_struct_M),1);          STD_feat_M = nan(numel(NoteFeatures_struct_M),1);    STE_feat_M = nan(numel(NoteFeatures_struct_M),1);  
    meanfeat_s = nan(numel(NoteFeatures_struct_s),1);          STD_feat_s = nan(numel(NoteFeatures_struct_s),1);    STE_feat_s = nan(numel(NoteFeatures_struct_s),1); 

        % Music-like
        for i = 1:numel(NoteFeatures_struct_M)    
            
            NF_M = NoteFeatures_struct_M{i};
            val_N_feat_M = getfield(NF_M,fieldname_string);
            
            hold on;
            xval = randn(numel(val_N_feat_M),1)/8 + 1;
            plot(xval,val_N_feat_M,'o','Color',[.75 .75 .75],'markersize',4); %scM.MarkerEdgeAlpha = .2;
            meanfeat_M(i) = nanmean(val_N_feat_M);
            STD_feat_M(i) = std(val_N_feat_M(~isnan(val_N_feat_M)));%/sqrt(length(val_N_feat_M));  
            STE_feat_M(i) = std(val_N_feat_M(~isnan(val_N_feat_M)))/sqrt(length(val_N_feat_M(~isnan(val_N_feat_M))));  
        end    
        
        % speech-like
        for i = 1:numel(NoteFeatures_struct_s)
                                    
            NF_s = NoteFeatures_struct_s{i};
            val_N_feat_s = getfield(NF_s,fieldname_string);
            
            hold on;
            xval = randn(numel(val_N_feat_s),1)/8 + 2;
            plot(xval,val_N_feat_s,'o','Color',[.75 .75 .75],'markersize',4);    
            meanfeat_s(i) = nanmean(val_N_feat_s);
            STD_feat_s(i) = std(val_N_feat_s(~isnan(val_N_feat_s)));%/sqrt(length(val_N_feat_s));  
            STE_feat_s(i) = std(val_N_feat_s(~isnan(val_N_feat_s)))/sqrt(length(val_N_feat_s(~isnan(val_N_feat_s))));  
        end
        % means on top in colors, STD as line
        hold on

            xval = linspace(.7,1.3,(numel(meanfeat_M)));
            [sorted_meanfeat_M,sortix] = sort(meanfeat_M);
            sorted_STD_feat_M = STD_feat_M(sortix);
            sorted_STE_feat_M = STE_feat_M(sortix);
            plot(xval,sorted_meanfeat_M,'ob','markersize',4.5,'linewidth',1.5);
            plot([.69 1.31],[mean(meanfeat_M) mean(meanfeat_M)],'k','linewidth',1.1)
                        
            for s = 1:numel(STE_feat_M)
%                 plot([xval(s) xval(s)],[sorted_meanfeat_M(s)-sorted_STD_feat_M(s) sorted_meanfeat_M(s)+sorted_STD_feat_M(s)],'b','linewidth',1.1)
                % better plot STE, not STD
                plot([xval(s) xval(s)],[sorted_meanfeat_M(s)-sorted_STE_feat_M(s) sorted_meanfeat_M(s)+sorted_STE_feat_M(s)],'b','linewidth',1.1)
            end
        
%             xval = randn(numel(meanfeat_s),1)/11 + 2;            
            xval = linspace(1.7,2.3,(numel(meanfeat_s)));
            [sorted_meanfeat_s,sortix] = sort(meanfeat_s);
            sorted_STD_feat_s = STD_feat_s(sortix);
            sorted_STE_feat_s = STE_feat_s(sortix);
            plot(xval,sorted_meanfeat_s,'or','markersize',4.5,'linewidth',1.5);
            plot([1.69 2.31],[mean(meanfeat_s) mean(meanfeat_s)],'k','linewidth',1.1)
            
            for s = 1:numel(STE_feat_s)
%                 plot([xval(s) , xval(s)],[sorted_meanfeat_s(s)-sorted_STD_feat_s(s) , sorted_meanfeat_s(s)+sorted_STD_feat_s(s)],'r','linewidth',1.1)
                % better plot STE, not STD
                plot([xval(s) , xval(s)],[sorted_meanfeat_s(s)-sorted_STE_feat_s(s) , sorted_meanfeat_s(s)+sorted_STE_feat_s(s)],'r','linewidth',1.1)
            end
%             plot(2.02,mean(meanfeat_s),'+k','markersize',16,'linewidth',3)
            
        if IS_SCALED
                ylim([0 1.1]);    
        else
        end
        
        xlim([0.5 2.5]);    
        xticks([1 2]);  xticklabels({'music-like','speech-like'});  box on
        ylabel(ylabelstring);   
