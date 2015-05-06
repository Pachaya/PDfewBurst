%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    GPM - VL Connnection
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rasterplot of VL  and GPm-VL chance of burst
SAVE_RASTER_FIG =0; 

% GPm VL
CntBurstSpkWT = zeros(ACT_Rec_size);
CntBurstSpkKO = zeros(ACT_Rec_size);
ChanceOfBurstingWT = zeros(ACT_Rec_size);
ChanceOfBurstingKO = zeros(ACT_Rec_size);
                    
p1_ii = 1; p2_ii=1;
for p3_ii = 1 : length(PARAM3)
    fLA1 = figure; set(fLA1,'position',[302          49        1290         948]) ;  set(gcf,'PaperPositionMode','auto')
    fLA2 = figure; set(fLA2,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
    fB1 = figure; set(fB1,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
    fB2 = figure; set(fB2,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
      cnt = 0;
    nR = length(PARAM4); nC = length(PARAM5);
    for p4_ii = 1 : length(PARAM4)
        for p5_ii = 1 : length(PARAM5)

            
            BasalAct = ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii}.VL;
            
            simTxt = ['VL ' titleTxt3 ' = ' num2str(PARAM3(p3_ii)) ', ' titleTxt4 ' = ' num2str(PARAM4(p4_ii)) ', ' titleTxt5 ' = ' num2str(PARAM5(p5_ii))];
            
            spkBin = BasalAct.WT.All.spktrain;
             baseline_fr_WT = get_avg_baseline(spkBin, CUTTIME, PhotoInjT);
%             [freq_all_WT,spkTime_all_WT, fg_handle_WT ] = raster_from_spkbin_BurstRange( spkBin,PhotoStop, Tstop,PhotoInjT, PhotoStop, DelayT, BurstRange, [ simTxt ' : WT']);
%             set(fg_handle_WT, 'position',[  449   450   791   528])
            
            spkBin = BasalAct.KO.All.spktrain;
             baseline_fr_KO = get_avg_baseline(spkBin, CUTTIME, PhotoInjT);
%             [freq_all_KO,spkTime_all_KO, fg_handle_KO ] = raster_from_spkbin_BurstRange( spkBin,PhotoStop, Tstop,PhotoInjT, PhotoStop, DelayT, BurstRange, [ simTxt ' : KO']);
%             set(fg_handle_KO, 'position',[  449   450   791   528])
            
            if (SAVE_RASTER_FIG)
                tmpTxt = get_Parameters_saveText(PARAMETERS, [3:5], [ p3_ii p4_ii p5_ii]);
                ffig = [ dirLoc dirFig 'RasterPlot_' tmpTxt];
                ttt = tic(); saveas( fg_handle_WT, [ffig '_WT.jpg'], 'jpg'); toc(ttt);
                ttt = tic(); saveas( fg_handle_WT, [ffig '_WT.fig'], 'fig'); toc(ttt);
                saveas( fg_handle_KO, [ffig '_KO.jpg'], 'jpg')
                saveas( fg_handle_KO, [ffig '_KO.fig'], 'fig')
                if (Close_Fig_aftr_save)
                   close(fg_handle_WT);                        close(fg_handle_KO);
                end
                
            end
            
            WT_expBaselineSpk = baseline_fr_WT/1000*BurstRange;
            WT_burstSpk = BasalAct.WT.All.BurstSpk; % BurstSpk  = Number of spike during BurstRange
            CntBurstSpkWT(p1_ii,p2_ii,p3_ii,p1_ii,p2_ii,p3_ii) = mean(WT_burstSpk - WT_expBaselineSpk); % Number of spike during BurstRange - expected number of spike from baseline activity
            
            KO_expBaselineSpk = baseline_fr_KO/1000*BurstRange;
            KO_burstSpk = BasalAct.KO.All.BurstSpk; % BurstSpk  = Number of spike during BurstRange
            CntBurstSpkKO(p1_ii,p2_ii,p3_ii,p4_ii,p5_ii) = mean(KO_burstSpk - KO_expBaselineSpk); % Number of spike during BurstRange - expected number of spike from baseline activity
            
            ChanceOfBurstingWT(p1_ii,p2_ii,p3_ii,p4_ii,p5_ii) =  sum( round(WT_burstSpk -WT_expBaselineSpk) > 0)/(nE+nI);
            ChanceOfBurstingKO(p1_ii,p2_ii,p3_ii,p4_ii,p5_ii) =  sum( round(KO_burstSpk -KO_expBaselineSpk) > 0)/(nE+nI);
            
            cnt = cnt+1;
            figure(fLA1);
            subplot( nR,nC,cnt)
            [nb,xb]=hist(BasalAct.WT.All.fr_AfterLightOff,10); bh=bar(xb,nb); hold on;
            set(bh,'facecolor',[0 0 0]);
            xlim([0 max(xb)+1])
            title( get_Parameters_titleText(PARAMETERS, [4,5], [ p4_ii, p5_ii ]))
            figure(fLA2);
            subplot( nR,nC,cnt)
            [nb,xb]=hist(BasalAct.KO.All.fr_AfterLightOff,10); bh=bar(xb,nb); hold on;
            set(bh,'facecolor',[1 0 0]);
            xlim([0 max(xb)+1])
            title( get_Parameters_titleText(PARAMETERS, [4,5], [ p4_ii, p5_ii ]))
            
            figure(fB1);
            subplot( nR,nC,cnt)
            [nb,xb]=hist(BasalAct.WT.All.BurstSpk); bh=bar(xb,nb); hold on;
            set(bh,'facecolor',[0 0 0]);
            xlim([0 max(xb)+1])
            title( get_Parameters_titleText(PARAMETERS, [4,5], [ p4_ii, p5_ii ]))
            figure(fB2);
            subplot( nR,nC,cnt)
            [nb,xb]=hist(BasalAct.KO.All.BurstSpk,10); bh=bar(xb,nb); hold on;
            set(bh,'facecolor',[1 0 0]);
            xlim([0 max(xb)+1])
            title( get_Parameters_titleText(PARAMETERS, [4,5], [ p4_ii, p5_ii ]))
                     
            clear BasalAct spkBin
        end
    end
    figure(fLA1); suptitle(['Avg Fr after light off WT: ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]);
    figure(fB1); suptitle(['Number of Burst spike WT: ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]);
    figure(fLA2); suptitle(['Avg Fr after light off KO: ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]);
    figure(fB2);  suptitle(['Number of Burst spike KO: ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]);

    if (SAVE_FIG)
        tmpTxt = get_Parameters_saveText(PARAMETERS, [3], [ p3_ii ]);
        ffig = [ dirLoc dirFig 'AvgFrLightOff' tmpTxt];
        saveas( fLA1, [ffig '_WT.jpg'], 'jpg')
        saveas( fLA1, [ffig '_WT.fig'], 'fig')
        saveas( fLA2, [ffig '_KO.jpg'], 'jpg')
        saveas( fLA2, [ffig '_KO.fig'], 'fig')        
        ffig = [ dirLoc dirFig 'NumBurstSpk'  tmpTxt];
        saveas( fB1, [ffig '_WT.jpg'], 'jpg')
        saveas( fB1, [ffig '_WT.fig'], 'fig')
        saveas( fB2, [ffig '_KO.jpg'], 'jpg')
        saveas( fB2, [ffig '_KO.fig'], 'fig')
   
        if (Close_Fig_aftr_save)
            close(fLA1); close(fLA2); close(fB1); close(fB2); 
        end
        
    end
end

cnt =0;
p1_ii = 1; p2_ii=1;
for p3_ii = 1 : length(PARAM3)
    Bfg = figure;  set(Bfg,'position',[680   283   712   695]); set(gcf,'PaperPositionMode','auto');
    Bfg2 = figure;  set(Bfg2,'position',[680   283   712   695]);  set(gcf,'PaperPositionMode','auto');
     nR = length(PARAM3); nC = length(PARAM4);
    for p4_ii = 1 : length(PARAM4)
        cnt = cnt+1;
        tmpTxt = get_Parameters_titleText(PARAMETERS, [4], [ p4_ii ]);
        figure(Bfg)
        subplot(nR,nC, cnt); hold on;
        plot( PARAM5, squeeze(CntBurstSpkWT(p1_ii,p2_ii,p3_ii,p4_ii,:)),'*-k');
        plot( PARAM5, squeeze(CntBurstSpkKO(p1_ii,p2_ii,p3_ii,p4_ii,:)),'*-r');
        legend('WT','KO','location','Best')
        title(tmpTxt);
        xlabel(titleTxt5); ylabel('Avg Burst spike')
        
        figure(Bfg2)
        subplot(nR,nC, cnt); hold on;
        plot( PARAM5, squeeze(ChanceOfBurstingWT(p1_ii,p2_ii,p3_ii,p4_ii,:)),'*-k');
        plot( PARAM5, squeeze(ChanceOfBurstingKO(p1_ii,p2_ii,p3_ii,p4_ii,:)),'*-r');
         ylim([ 0 0.7])
        legend('WT','KO','location','Best')
        title(tmpTxt);
        xlabel(titleTxt5); ylabel('Chance of Burst spike')
        
    end
     tmpTxt = get_Parameters_titleText(PARAMETERS, [3], [ p3_ii ]);
    figure(Bfg)
    suptitle(tmpTxt);
    figure(Bfg2)
    suptitle([ 'Chance of Bursting, '  tmpTxt]);
    
    if (SAVE_FIG)
        tmpTxt = get_Parameters_saveText(PARAMETERS, [3], [ p3_ii ]);
        ffig = [ dirLoc dirFig 'NumBurstSpk' tmpTxt ];
        saveas(  Bfg, [ffig '.jpg'], 'jpg');        saveas(  Bfg, [ffig '.fig'], 'fig');
        ffig = [ dirLoc dirFig 'ChanceOfBurst' tmpTxt];
        saveas(  Bfg2, [ffig '.jpg'], 'jpg');        saveas(  Bfg2, [ffig '.fig'], 'fig');
        if (Close_Fig_aftr_save)
            close(Bfg); close(Bfg2);
        end
    end
end

%%
%Make only the neuron activity during burst of VL layer
Close_Fig_aftr_save = 1;

p1_ii = 1; p2_ii=1;
for p3_ii = 1 : length(PARAM3)
    fRB = figure; set(fRB,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
    fBspk = figure; set(fBspk,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
    fcum = figure; set(fcum,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
    fcumspk = figure; set(fcumspk,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
    
    cnt = 0;
    nR = length(PARAM4); nC = length(PARAM5);
    if(nR == 1)&&(nC ==1)
         set(fRB,'position',[ 1228         336         592         398]);
         set(fBspk,'position',[ 1228         336         592         398]);
         set(fcum,'position',[ 1228         336         592         398]);
         set(fcumspk,'position',[ 1228         336         592         398]);
    end
    for p4_ii = 1 : length(PARAM4)
        for p5_ii = 1 : length(PARAM5)
            la_ii = p4_ii; m_ii = p3_ii; s_ii = 1; ld_ii = p5_ii;
            
            cnt = cnt +1;
  BasalAct = ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii}.VL;
            hbin = 2 ;
            histBin = 1:hbin:size(BasalAct.WT.All.BurstSpkTrain,2); 
          


            % Distribution of spikes in bursting range 
           
             tmpRebound_WT = sum(BasalAct.WT.All.BurstSpkTrain);
             tmpRebound_KO = sum(BasalAct.KO.All.BurstSpkTrain);
             
            ll =length(tmpRebound_WT);
            if mod(ll,hbin) == 0
                % Can use reshape fn.
                tmpRebound_WT = sum(reshape(tmpRebound_WT,[hbin ll/hbin]));
                tmpRebound_KO = sum(reshape(tmpRebound_KO,[hbin ll/hbin]));
                xx = histBin;
            else %If not, just set bin to 1
                xx = 1:1:size(BasalAct.KO.All.BurstSpkTrain,2);
            end
            
            figure(fBspk); 
            subplot(nR,nC,cnt);
            pH=bar(xx,tmpRebound_WT/ncells*100,'k'); hold on; 
            cH = get(pH,'Children');
            set(cH,'FaceAlpha',0.5); % 0 = transparent, 1 = opaque.
            pH=bar(xx,tmpRebound_KO/ncells*100,'FaceColor','r','EdgeColor','r');
           cH = get(pH,'Children');
            set(cH,'FaceAlpha',0.5); % 0 = transparent, 1 = opaque.
            xlim([0 BurstRange+10]);
            ylabel('% Neurons'); xlabel('Time(ms)')
            title([ titleTxt4 ' = ' num2str(PARAM4(p4_ii)) ', ' titleTxt5 ' = ' num2str(PARAM5(p5_ii))])
            legend('WT','KO','location','northwest')
             
            NWT = tmpRebound_WT; NKO = tmpRebound_KO;
            figure(fcumspk); 
            subplot(nR,nC,cnt);
            plot(xx,cumsum(NWT)./sum(NWT),'k'); hold on;
            plot(xx,cumsum(NKO)./sum(NKO),'r'); hold on;
            xlim([0 BurstRange+10]);
            xlabel('Latency of Spikes in bursting range'); ylabel('Cumulative Probability');
           legend('WT','KO','location','northwest')
           
              % Distribution of first rebound spike timing 
            tmpRebound_WT = zeros(ncells,1); tmpRebound_KO = zeros(ncells,1);
            for id = 1 : ncells
                tmp = find(BasalAct.WT.All.BurstSpkTrain(id,:),1,'first');
                if ~isempty(tmp) 
                    tmpRebound_WT(id) = tmp;                  
                end
                
                tmp = find(BasalAct.KO.All.BurstSpkTrain(id,:),1,'first');
                if ~isempty(tmp) 
                    tmpRebound_KO(id) = tmp;                  
                end                
            end      
             tmpRebound_WT = tmpRebound_WT(tmpRebound_WT ~= 0);
             tmpRebound_KO = tmpRebound_KO(tmpRebound_KO ~= 0);
            figure(fRB); 
            
             histBin = xx;
            subplot(nR,nC,cnt);
            [NWT, X] = hist(tmpRebound_WT,histBin);
            pH= bar(X,NWT/ncells*100,'FaceColor','k','EdgeColor','k'); hold on; 
            cH = get(pH,'Children');
            set(cH,'FaceAlpha',0.5); % 0 = transparent, 1 = opaque.

            [NKO,X] = hist(tmpRebound_KO,histBin);
            pH=bar(X,NKO/ncells*100,'FaceColor','r','EdgeColor','r'); hold on; 
            cH = get(pH,'Children');
            set(cH,'FaceAlpha',0.5); % 0 = transparent, 1 = opaque.
            xlim([0 BurstRange+10]);
            ylabel('% Neurons'); xlabel('Time(ms)')
            title([ titleTxt4 ' = ' num2str(PARAM4(p4_ii)) ', ' titleTxt5 ' = ' num2str(PARAM5(p5_ii))])
             legend('WT','KO','location','northwest')

             
            figure(fcum); 
            subplot(nR,nC,cnt);            
            plot(X, cumsum(NWT)./sum(NWT),'k'); hold on;
             plot(X, cumsum(NKO)./sum(NKO),'r'); hold on; 
             xlim([0 BurstRange+10]);
             xlabel('First rebound spike timing'); ylabel('Cumulative Probability');
         legend('WT','KO','location','northwest')
          
            clear BasalAct spkBin
        end
    end
    figure(fRB); suptitle(['Latency of first Rebound spike : ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]);  legend('WT','KO','location','northwest') % Note equal to the latency of peak firing rate --> latency of peak firing rate calculate from instantaneous firing rate
    figure(fcum); suptitle(['Latency of first Rebound spike : ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]); legend('WT','KO','location','northwest')
    figure(fBspk); suptitle(['Percent of neuron that fire during bursting period : ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]); legend('WT','KO','location','northwest') % Note equal to the latency of peak firing rate --> latency of peak firing rate calculate from instantaneous firing rate
    figure(fcumspk); suptitle(['Percent of neuron that fire during bursting period : ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]); legend('WT','KO','location','northwest') %

end

if (SAVE_FIG)
        tmpTxt = get_Parameters_saveText(PARAMETERS, [3], [ p3_ii ]);
        ffig = [ dirLoc dirFig 'FirstReboundSpike'  tmpTxt];
        saveas(fRB, [ffig '.fig'] , 'fig');
        saveas(fRB, [ffig '.jpg'] , 'jpg');
        ffig = [ dirLoc dirFig 'AllBurstSpike'  tmpTxt]; % Better minus by expected baseline activity
        saveas(fBspk, [ffig '.fig'] , 'fig');
        saveas(fBspk, [ffig '.jpg'] , 'jpg');
        ffig = [ dirLoc dirFig 'FirstReboundSpikeCumSum'  tmpTxt];
        saveas(fcum, [ffig '.fig'] , 'fig');
        saveas(fcum, [ffig '.jpg'] , 'jpg');
        ffig = [ dirLoc dirFig 'AllBurstSpikeCumSum'  tmpTxt]; % Better minus by expected baseline activity
        saveas(fcumspk, [ffig '.fig'] , 'fig');
        saveas(fcumspk, [ffig '.jpg'] , 'jpg');
        
        
        if (Close_Fig_aftr_save)
            close(fRB); close(fBspk); close(fcum); close(fcumspk);
        end 
        
end

%%
% for Poster , only one case    

%Make only the neuron activity during burst of VL layer
Close_Fig_aftr_save = 0;
TimeofFirstReboundSpk = zeros(ACT_Rec_size);
p1_ii = 1; p2_ii=1; p3_ii = 1; p4_ii=1; p5_ii = 1;   
         set(fRB,'position',[ 1228         336         592         398]);
         set(fBspk,'position',[ 1228         336         592         398]);
         set(fcum,'position',[ 1228         336         592         398]);
         set(fcumspk,'position',[ 1228         336         592         398]);

            la_ii = p3_ii; m_ii = p4_ii; s_ii = p5_ii;
            

            BasalAct = ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii}.VL;
            hbin = 2 ;
            histBin = 1:hbin:size(BasalAct.WT.All.BurstSpkTrain,2); 
          


            % Distribution of spikes in bursting range 
           
             tmpRebound_WT = sum(BasalAct.WT.All.BurstSpkTrain);
             tmpRebound_KO = sum(BasalAct.KO.All.BurstSpkTrain);
             
            ll =length(tmpRebound_WT);
            if mod(ll,hbin) == 0
                % Can use reshape fn.
                tmpRebound_WT = sum(reshape(tmpRebound_WT,[hbin ll/hbin]));
                tmpRebound_KO = sum(reshape(tmpRebound_KO,[hbin ll/hbin]));
                xx = histBin;
            else %If not, just set bin to 1
                xx = 1:1:size(BasalAct.KO.All.BurstSpkTrain,2);
            end
            
            figure(fBspk); 
           
            pH=bar(xx,tmpRebound_WT/ncells*100,'k'); hold on; 
            cH = get(pH,'Children');
            set(cH,'FaceAlpha',0.5); % 0 = transparent, 1 = opaque.
            pH=bar(xx,tmpRebound_KO/ncells*100,'FaceColor','r','EdgeColor','r');
           cH = get(pH,'Children');
            set(cH,'FaceAlpha',0.5); % 0 = transparent, 1 = opaque.
            xlim([0 BurstRange+10]);
            ylabel('% Neurons'); xlabel('Time(ms)')
            title([ titleTxt4 ' = ' num2str(PARAM4(p4_ii)) ', ' titleTxt5 ' = ' num2str(PARAM5(p5_ii))])
            legend('WT','KO','location','northwest')
             
            NWT = tmpRebound_WT; NKO = tmpRebound_KO;
            figure(fcumspk); 
  
            plot(xx,cumsum(NWT)./sum(NWT),'k'); hold on;
            plot(xx,cumsum(NKO)./sum(NKO),'r'); hold on;
            xlim([0 BurstRange+10]);
            xlabel('Latency of Spikes in bursting range'); ylabel('Cumulative Probability');
           legend('WT','KO','location','northwest')
           
              % Distribution of first rebound spike timing 
            tmpRebound_WT = zeros(ncells,1); tmpRebound_KO = zeros(ncells,1);
            for id = 1 : ncells
                tmp = find(BasalAct.WT.All.BurstSpkTrain(id,:),1,'first');
                if ~isempty(tmp) 
                    tmpRebound_WT(id) = tmp;                  
                end
                
                tmp = find(BasalAct.KO.All.BurstSpkTrain(id,:),1,'first');
                if ~isempty(tmp) 
                    tmpRebound_KO(id) = tmp;                  
                end                
            end      
             tmpRebound_WT = tmpRebound_WT(tmpRebound_WT ~= 0);
             tmpRebound_KO = tmpRebound_KO(tmpRebound_KO ~= 0);
            figure(fRB);             
             histBin = xx;
          
            [NWT, X] = hist(tmpRebound_WT,histBin);
            pH= bar(X,NWT/ncells*100,'FaceColor','k','EdgeColor','k'); hold on; 
            cH = get(pH,'Children');
            set(cH,'FaceAlpha',0.5); % 0 = transparent, 1 = opaque.

            [NKO,X] = hist(tmpRebound_KO,histBin);
            pH=bar(X,NKO/ncells*100,'FaceColor','r','EdgeColor','r'); hold on; 
            cH = get(pH,'Children');
            set(cH,'FaceAlpha',0.5); % 0 = transparent, 1 = opaque.
            xlim([0 BurstRange+10]);
            ylabel('% Neurons'); xlabel('Time(ms)')
            title([ titleTxt4 ' = ' num2str(PARAM4(p4_ii)) ', ' titleTxt5 ' = ' num2str(PARAM5(p5_ii))])
             legend('WT','KO','location','northwest')

             
            figure(fcum); 
    
            plot(X, cumsum(NWT)./sum(NWT),'k'); hold on;
             plot(X, cumsum(NKO)./sum(NKO),'r'); hold on; 
             xlim([0 BurstRange+10]);
             xlabel('First rebound spike timing'); ylabel('Cumulative Probability');
         legend('WT','KO','location','northwest')

