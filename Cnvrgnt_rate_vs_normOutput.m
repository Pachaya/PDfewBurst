
%% M1 baseline activity
dataYLeg = {CnvrgntTypes{1}.leg, CnvrgntTypes{2}.leg, CnvrgntTypes{3}.leg};
dataYLeg2 = {CnvrgntTypes{1}.leg, CnvrgntTypes{2}.leg, CnvrgntTypes{3}.leg};
dataY1 = cell(NumCnvrgntTypes ,1);  dataY2 = cell(NumCnvrgntTypes ,1);

    
    
    PhotoInjT  = TSTOP;
    if ~exist('sumW_LST','var')  || ~exist('VLperM1','var')
        Get_info_TC_convergence_numVLperM1
    end
    
    for p3_ii = 1 : length(PARAM3) 

        figSize =  [  1         401        1920         963];
        
        fscatterNumVL1 =   figure; set(gcf, 'position', figSize); set(gcf,'PaperPositionMode','auto');
        fscatterNumVL2 =   figure; set(gcf, 'position', figSize); set(gcf,'PaperPositionMode','auto');
        
        cntcnt = 0;
        nR = length(PARAM4); nC = length(PARAM5);
        
        
        for p4_ii = 1 : length(PARAM4)
%             newFgforAmp = figure; 
            for p5_ii = 1 : length(PARAM5)
                
            for ct_ii = 1 : NumCnvrgntTypes
                ACT_Record = CnvrgntTypes{ct_ii}.ACT_Record;
                 if ct_ii == 1 
                    dirLoc = [PATH 'OscInput_Sim/'];
                    dirFig = ['../OscInput_varyTCtype_Sim/'    tmpdirFig ];    
                 else     
                    dirLoc = [PATH 'OscInput_varyTCtype_Sim/']; % [PATH 'OscInput_Sim/']; % /OscInput_varyTCtype_Sim/
                    dirFig = [   tmpdirFig ];
                 end
                
                tmpWTtrain_VL =  ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii,p6_ii}.VL.WT.All.spktrain;
                tmpKOtrain_VL =  ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii,p6_ii}.VL.KO.All.spktrain;
                tmpT = CUTTIME +1 : PhotoInjT;
                VLavgFrWT = mean(sum(tmpWTtrain_VL(:,tmpT),2) /length(tmpT)*1000);
                VLavgFrKO = mean(sum(tmpKOtrain_VL(:,tmpT),2) /length(tmpT)*1000);             
                
                M1_baselineWT = zeros(length(PARAM1),length(PARAM2));
                M1_baselineKO = zeros(length(PARAM1),length(PARAM2));
                tmpT = CUTTIME +1 : PhotoInjT;
                normalActrange =  length(tmpT);
                for p1_ii = 1 : length(PARAM1)
                    for p2_ii = 1 : length(PARAM2)
                        BasalAct = ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii}.M1; %% for M1
                        M1_baselineWT(p1_ii,p2_ii) = mean(mean(BasalAct.WT.All.spktrain(:,tmpT))).*1000 -VLavgFrWT;
                        M1_baselineKO(p1_ii,p2_ii) = mean(mean(BasalAct.KO.All.spktrain(:,tmpT))).*1000- VLavgFrKO;
                    end
                end
                dataY1{ct_ii} = M1_baselineWT;   dataY2{ct_ii} = M1_baselineKO;
            end
                
                
      cntcnt = cntcnt + 1;
      figure; cntcnt = 1;
             xt = PARAM2; yt = round(VLperM1);
            xl = 'Average sum of weight per cell'; yl =  'Normalized Neuron Activity at target layer(Hz)';
            suptitleTxt = {get_Parameters_titleText(PARAMETERS, [4:5], [ p4_ii,p5_ii])};
          
    
    figure(fscatterNumVL1); subplot(nR,nC,cntcnt);  plot_scatter( sumW_LST, dataY1, dataYLeg, suptitleTxt ,xl,yl);
    figure(fscatterNumVL2); subplot(nR,nC,cntcnt);  plot_scatter( sumW_LST, dataY2, dataYLeg, suptitleTxt ,xl,yl);
    
    
%      figure(newFgforAmp); plot_scatter_dotStyl( sumW_LST, dataY2, dataYLeg, suptitleTxt ,xl,yl,[ getLineStyle(p5_ii) getDotStyle(p5_ii)); hold on;
    
            end
        end
    end
    
    
    htt  = 'Rate of convergent vs. Target Activity'; 
    figure(fscatterNumVL1); suptitle({htt, [ tt1 ', ' CtypeTxt]}); set(gcf,'PaperPositionMode','auto');
    figure(fscatterNumVL2); suptitle({htt,[ tt2  ', ' CtypeTxt]}); set(gcf,'PaperPositionMode','auto');
    
    
    if(SAVE_FIG)
        tmpTxt = ['normVL' get_Parameters_saveText(PARAMETERS, [3], [ p3_ii]) '_allTypes'  ];
        
        
        fg = fscatterNumVL1;
        figname = ['Scatter_avgSumW_normM1actWT' tmpTxt];
        saveas(fg,[dirLoc dirFig figname '.fig'],'fig'); saveas(fg,[dirLoc dirFig figname '.jpg'],'jpg');
        fg = fscatterNumVL2;
        figname = ['Scatter_avgSumW_normM1actKO' tmpTxt];
        saveas(fg,[dirLoc dirFig figname '.fig'],'fig'); saveas(fg,[dirLoc dirFig figname '.jpg'],'jpg');
        
    end
    

