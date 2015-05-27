
%% M1 baseline activity

% For param 3, 4 , 5
PhotoInjT  = TSTOP;

for p3_ii = 1 : length(PARAM3)
    figSize =  [ 16         138        2551        1158];
    fCombine1 =  figure; set(gcf, 'position',figSize); set(gcf,'PaperPositionMode','auto');
    fCombine2 =  figure; set(gcf, 'position', figSize); set(gcf,'PaperPositionMode','auto');
    fContour1 =  figure; set(gcf, 'position', figSize); set(gcf,'PaperPositionMode','auto');
    fContour2 =  figure; set(gcf, 'position', figSize); set(gcf,'PaperPositionMode','auto');
    
    fCombine1numVL =  figure; set(gcf, 'position',figSize); set(gcf,'PaperPositionMode','auto');
    fCombine2numVL =  figure; set(gcf, 'position', figSize); set(gcf,'PaperPositionMode','auto');
    fContour1numVL =  figure; set(gcf, 'position', figSize); set(gcf,'PaperPositionMode','auto');
    fContour2numVL =  figure; set(gcf, 'position', figSize); set(gcf,'PaperPositionMode','auto');
    cntcnt = 0;
    nR = length(PARAM4); nC = length(PARAM5);
    for p4_ii = 1 : length(PARAM4)
        for p5_ii = 1 : length(PARAM5)
            M1_baselineWT = zeros(length(PARAM1),length(PARAM2));
            M1_baselineKO = zeros(length(PARAM1),length(PARAM2));
            tmpT = CUTTIME +1 : PhotoInjT;
            normalActrange =  length(tmpT);
            for p1_ii = 1 : length(PARAM1)
                for p2_ii = 1 : length(PARAM2)
                    BasalAct = ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii}.M1; %% for M1
                    M1_baselineWT(p1_ii,p2_ii) = mean(mean(BasalAct.WT.All.spktrain(:,tmpT))).*1000;
                    M1_baselineKO(p1_ii,p2_ii) = mean(mean(BasalAct.KO.All.spktrain(:,tmpT))).*1000;
                end
            end
            
            SAVE_FIG = 1;
            figLoc =[  194         333        1760         719];
            tt1 = 'WT'; tt2 = 'KO';
            xt = PARAM2; yt = PARAM1;
            xl = titleTxt2; yl =titleTxt1; %yl =  'Average number of VL per M1 (cell)';
            suptitleTxt = {'Average firing rate of M1 baseline activity', get_Parameters_titleText(PARAMETERS, [3:5], [ p3_ii,p4_ii,p5_ii])};
            %suptitleTxt = get_Parameters_titleText(PARAMETERS, [4:5], [p4_ii,p5_ii]);
            [ fg_M1FR ] = plotFigure_paramMat( M1_baselineWT,M1_baselineKO, tt1,tt2, xl,yl, xt,yt, suptitleTxt,figLoc );
            cntcnt = cntcnt+1;
            Ncontour =4;
            data1=M1_baselineWT; data2 = M1_baselineKO;
            figure(fCombine1); subplot(nR,nC,cntcnt);  plot_paramMat( data1, suptitleTxt ,xl,yl, xt,yt); 
            figure(fCombine2); subplot(nR,nC,cntcnt);  plot_paramMat( data2, suptitleTxt ,xl,yl, xt,yt); 
            figure(fContour1); subplot(nR,nC,cntcnt);  plot_contourparamMat( data1, suptitleTxt ,xl,yl, xt,yt, Ncontour)
            figure(fContour2); subplot(nR,nC,cntcnt);  plot_contourparamMat( data2, suptitleTxt ,xl,yl, xt,yt, Ncontour)
            
                        
            if(SAVE_FIG)
                fg = fg_M1FR;
                tmpTxt = get_Parameters_saveText(PARAMETERS, [1:5], [ p1_ii,p2_ii, p3_ii,p4_ii, p5_ii]);
                figname = ['M1FR_atBaseline' tmpTxt];
                saveas(fg,[dirLoc dirFig figname '.fig'],'fig');
                saveas(fg,[dirLoc dirFig figname '.jpg'],'jpg');
            end
            
            xt = PARAM2; yt = round(VLperM1);
            xl = titleTxt2; yl =  'Average number of VL per M1 (cell)';
            suptitleTxt = {'Average firing rate of M1 baseline activity', get_Parameters_titleText(PARAMETERS, [3:5], [ p3_ii,p4_ii,p5_ii])};
            [ fg_M1FRnumVL ] = plotFigure_paramMat( M1_baselineWT,M1_baselineKO, tt1,tt2, xl,yl, xt,yt, suptitleTxt,figLoc );
%                         Ncontour =4;
            data1=M1_baselineWT; data2 = M1_baselineKO;
            figure(fCombine1numVL); subplot(nR,nC,cntcnt);  plot_paramMat( data1, suptitleTxt ,xl,yl, xt,yt); 
            figure(fCombine2numVL); subplot(nR,nC,cntcnt);  plot_paramMat( data2, suptitleTxt ,xl,yl, xt,yt); 
            figure(fContour1numVL); subplot(nR,nC,cntcnt);  plot_contourparamMat( data1, suptitleTxt ,xl,yl, xt,yt, Ncontour)
            figure(fContour2numVL); subplot(nR,nC,cntcnt);  plot_contourparamMat( data2, suptitleTxt ,xl,yl, xt,yt, Ncontour)
            
            if(SAVE_FIG)
                fg = fg_M1FRnumVL;
                tmpTxt = get_Parameters_saveText(PARAMETERS, [1:5], [ p1_ii,p2_ii, p3_ii,p4_ii, p5_ii]);
                figname = ['M1FR_atBaselineNumVL' tmpTxt];
                saveas(fg,[dirLoc dirFig figname '.fig'],'fig');
                saveas(fg,[dirLoc dirFig figname '.jpg'],'jpg');
            end
        end
    end
    tt1 = 'WT'; tt2 = 'KO';  txtp3 = get_Parameters_titleText(PARAMETERS, [3], [p3_ii]);
    stt = 'Average firing rate of M1 baseline activity';
    
       figure(fCombine1);  suptitle({stt, [ tt1 ' : ' txtp3] })
            figure(fCombine2); suptitle({stt, [ tt2 ' : ' txtp3] })
            figure(fContour1);  suptitle({stt, [ tt1 ' : ' txtp3] })
            figure(fContour2);  suptitle({stt, [ tt2 ' : ' txtp3] })
         figure(fCombine1numVL);  suptitle({stt, [ tt1 ' : ' txtp3] })
            figure(fCombine2numVL); suptitle({stt, [ tt2 ' : ' txtp3] }); set(gcf,'PaperPositionMode','auto');
            figure(fContour1numVL);  suptitle({stt, [ tt1 ' : ' txtp3] }) ; 
            figure(fContour2numVL);  suptitle({stt, [ tt2 ' : ' txtp3] }); set(gcf,'PaperPositionMode','auto');
                  
            
            if(SAVE_FIG)
                tmpTxt = get_Parameters_saveText(PARAMETERS, [3], [ p3_ii]);  

                fg = fCombine1;                
                figname = ['CombineM1FR_atBaselineWT' tmpTxt];
                saveas(fg,[dirLoc dirFig figname '.fig'],'fig'); saveas(fg,[dirLoc dirFig figname '.jpg'],'jpg');
                fg = fCombine2;                
                figname = ['CombineM1FR_atBaselineKO' tmpTxt];
                saveas(fg,[dirLoc dirFig figname '.fig'],'fig'); saveas(fg,[dirLoc dirFig figname '.jpg'],'jpg');
                fg = fContour1;                
                figname = ['ContourM1FR_atBaselineWT' tmpTxt];
                saveas(fg,[dirLoc dirFig figname '.fig'],'fig'); saveas(fg,[dirLoc dirFig figname '.jpg'],'jpg');
                fg = fContour2;                
                figname = ['ContourM1FR_atBaselineKO' tmpTxt];
                saveas(fg,[dirLoc dirFig figname '.fig'],'fig'); saveas(fg,[dirLoc dirFig figname '.jpg'],'jpg');
                fg = fCombine1numVL;                
                figname = ['CombineM1FRnumVL_atBaselineWT' tmpTxt];
                saveas(fg,[dirLoc dirFig figname '.fig'],'fig'); saveas(fg,[dirLoc dirFig figname '.jpg'],'jpg');
                fg = fCombine2numVL;                
                figname = ['CombineM1FRnumVL_atBaselineKO' tmpTxt];
                saveas(fg,[dirLoc dirFig figname '.fig'],'fig'); saveas(fg,[dirLoc dirFig figname '.jpg'],'jpg');
                fg = fContour1numVL;                
                figname = ['ContourM1FRnumVL_atBaselineWT' tmpTxt];
                saveas(fg,[dirLoc dirFig figname '.fig'],'fig'); saveas(fg,[dirLoc dirFig figname '.jpg'],'jpg');
                fg = fContour2numVL;                
                figname = ['ContourM1FRnumVL_atBaselineKO' tmpTxt];
                saveas(fg,[dirLoc dirFig figname '.fig'],'fig'); saveas(fg,[dirLoc dirFig figname '.jpg'],'jpg');
            end
     
end
