% Different  in M1 Activity in various condition

% For param 3, 4 , 5
PhotoInjT  = TSTOP;
p3_ii = 1; p6_ii = 1;

% for p3_ii = 1 : length(PARAM3)  % Trial
dataY1 =[]; data2 = [];
if ~ exist('CellsOfMatrixWT','var')
    CellsOfMatrixWT = cell(length(PARAM4) , length(PARAM5));
    CellsOfMatrixKO = cell(length(PARAM4) , length(PARAM5));
    for p4_ii = 1 : length(PARAM4) % OSC F
        for p5_ii = 1 : length(PARAM5) % OSC Amp
            p1_ii = 1; p2_ii = 1;
            
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
                    M1_baselineKO(p1_ii,p2_ii) = mean(mean(BasalAct.KO.All.spktrain(:,tmpT))).*1000 - VLavgFrKO;
                end
            end
            CellsOfMatrixWT{p4_ii, p5_ii} =  M1_baselineWT;
            CellsOfMatrixKO{p4_ii,p5_ii} =  M1_baselineKO;
        end
    end
end
%% Diff in Osc input 40Hz - 20Hz

figSize = [145         577        1614         655];
figSize2 =  [  57         634        1624         590];
fCombine1numVL =  figure; set(gcf, 'position',figSize); set(gcf,'PaperPositionMode','auto');
fCombine2numVL =  figure; set(gcf, 'position', figSize); set(gcf,'PaperPositionMode','auto');
fContour1numVL =  figure; set(gcf, 'position', figSize2); set(gcf,'PaperPositionMode','auto');
fContour2numVL =  figure; set(gcf, 'position', figSize2); set(gcf,'PaperPositionMode','auto');
cntcnt = 0;
nR = 2; nC = length(PARAM5)-1;

for p5_ii = 2 : length(PARAM5) % OSC Amp
    cntcnt = cntcnt + 1;
    
    xt = PARAM2; yt = round(VLperM1);
    xl = titleTxt2; yl =  'Average number of VL per M1 (cell)';
    suptitleTxt = get_Parameters_titleText(PARAMETERS, [5], [ p5_ii]);
    
    tmpMatAmp =  CellsOfMatrixWT{2, p5_ii}; tmpMat0 =  CellsOfMatrixWT{1,1};
    diff_M1_WT =  tmpMatAmp - tmpMat0;
    data1= diff_M1_WT;
    figure(fCombine1numVL); subplot(nR,nC,cntcnt);  plot_paramMat( data1, suptitleTxt ,xl,yl, xt,yt);
    figure(fContour1numVL); subplot(nR,nC,cntcnt);  plot_contourparamMat( data1, suptitleTxt ,xl,yl, xt,yt, Ncontour);
    
    tmpMatAmp =  CellsOfMatrixKO{2, p5_ii}; tmpMat0 =  CellsOfMatrixKO{1, p5_ii};
    diff_M1_KO =  tmpMatAmp - tmpMat0;
    data2= diff_M1_KO;
    figure(fCombine2numVL); subplot(nR,nC,cntcnt);  plot_paramMat( data2, suptitleTxt ,xl,yl, xt,yt);
    figure(fContour2numVL); subplot(nR,nC,cntcnt);  plot_contourparamMat( data2, suptitleTxt ,xl,yl, xt,yt, Ncontour);
    
    
    
end
tt1 = '[WT] '; tt2 = '[KO] ';
stt = 'Different in M1 activity when oscillation frequency increase (40Hz - 20Hz)';
figure(fCombine1numVL);  suptitle( {stt, [tt1 CtypeTxt]}); set(gcf,'PaperPositionMode','auto');
figure(fCombine2numVL);  suptitle( {stt, [tt2 CtypeTxt]}); set(gcf,'PaperPositionMode','auto');
figure(fContour1numVL);  suptitle( {stt, [tt1 CtypeTxt]}); set(gcf,'PaperPositionMode','auto');
figure(fContour2numVL);  suptitle( {stt, [tt2 CtypeTxt]}); set(gcf,'PaperPositionMode','auto');

if(SAVE_FIG)
    tmpTxt = ['normVL' get_Parameters_saveText(PARAMETERS, [3], [ p3_ii]) '_' codeTxt ];
    fg = fCombine1numVL;
    figname = ['DiffM1FR_cmprOscWT' tmpTxt];
    saveas(fg,[dirLoc dirFig figname '.fig'],'fig'); saveas(fg,[dirLoc dirFig figname '.jpg'],'jpg');
    fg = fCombine2numVL;
    figname = ['DiffM1FR_cmprOscKO' tmpTxt];
    saveas(fg,[dirLoc dirFig figname '.fig'],'fig'); saveas(fg,[dirLoc dirFig figname '.jpg'],'jpg');
    fg = fContour1numVL;
    figname = ['DiffM1FRcontour_cmprOscWT' tmpTxt];
    saveas(fg,[dirLoc dirFig figname '.fig'],'fig'); saveas(fg,[dirLoc dirFig figname '.jpg'],'jpg');
    fg = fContour2numVL;
    figname = ['DiffM1FRcontour_cmprOscKO' tmpTxt];
    saveas(fg,[dirLoc dirFig figname '.fig'],'fig'); saveas(fg,[dirLoc dirFig figname '.jpg'],'jpg');
end

%%  Diff in Osc input amp

figSize = [246   577   975   739];
figSize2 =  [ 246   503   800   787];
fCombine1numVL =  figure; set(gcf, 'position',figSize); set(gcf,'PaperPositionMode','auto');
fCombine2numVL =  figure; set(gcf, 'position', figSize); set(gcf,'PaperPositionMode','auto');
fContour1numVL =  figure; set(gcf, 'position', figSize2); set(gcf,'PaperPositionMode','auto');
fContour2numVL =  figure; set(gcf, 'position', figSize2); set(gcf,'PaperPositionMode','auto');
cntcnt = 0;
nR = length(PARAM4); nC = 2;

for p4_ii = 1 : length(PARAM4) % OSC F
    for amp = 2 : 3
        cntcnt = cntcnt + 1;
        
        xt = PARAM2; yt = round(VLperM1);
        xl = titleTxt2; yl =  'Average number of VL per M1 (cell)';
        suptitleTxt = get_Parameters_titleText(PARAMETERS, [4,5], [ p4_ii,amp]);
        
        tmpMat1 =  CellsOfMatrixWT{p4_ii, amp}; tmpMat0 =  CellsOfMatrixWT{p4_ii, 1};
        diff_M1_WT =  tmpMat1 - tmpMat0;
        data1= diff_M1_WT;
        figure(fCombine1numVL); subplot(nR,nC,cntcnt);  plot_paramMat( data1, suptitleTxt ,xl,yl, xt,yt);
        figure(fContour1numVL); subplot(nR,nC,cntcnt);  plot_contourparamMat( data1, suptitleTxt ,xl,yl, xt,yt, Ncontour)
        
        tmpMat1 =  CellsOfMatrixKO{p4_ii, amp}; tmpMat0 =  CellsOfMatrixWT{p4_ii, 1};
        diff_M1_KO =  tmpMat1 - tmpMat0;
        data2= diff_M1_KO;
        figure(fCombine2numVL); subplot(nR,nC,cntcnt);  plot_paramMat( data2, suptitleTxt ,xl,yl, xt,yt);
        figure(fContour2numVL); subplot(nR,nC,cntcnt);  plot_contourparamMat( data2, suptitleTxt ,xl,yl, xt,yt, Ncontour)
    end
end
tt1 = '[WT] '; tt2 = '[KO] ';
stt = 'Different in M1 activity when oscillation amplitude increase (compare with non-oscillating case)';
figure(fCombine1numVL);  suptitle( {stt, [tt1 CtypeTxt]}); set(gcf,'PaperPositionMode','auto');
figure(fCombine2numVL);  suptitle( {stt, [tt2 CtypeTxt]}); set(gcf,'PaperPositionMode','auto');
figure(fContour1numVL);  suptitle( {stt, [tt1 CtypeTxt]}); set(gcf,'PaperPositionMode','auto');
figure(fContour2numVL);  suptitle( {stt, [tt2 CtypeTxt]}); set(gcf,'PaperPositionMode','auto');

if(SAVE_FIG)
    tmpTxt = ['normVL' get_Parameters_saveText(PARAMETERS, [3], [ p3_ii]) '_' codeTxt];
    fg = fCombine1numVL;
    figname = ['DiffM1FR_cmprOscAmpWT' tmpTxt];
    saveas(fg,[dirLoc dirFig figname '.fig'],'fig'); saveas(fg,[dirLoc dirFig figname '.jpg'],'jpg');
    fg = fCombine2numVL;
    figname = ['DiffM1FR_cmprOscAmpKO' tmpTxt];
    saveas(fg,[dirLoc dirFig figname '.fig'],'fig'); saveas(fg,[dirLoc dirFig figname '.jpg'],'jpg');
    fg = fContour1numVL;
    figname = ['DiffM1FRcontour_cmprOsAmpcWT' tmpTxt];
    saveas(fg,[dirLoc dirFig figname '.fig'],'fig'); saveas(fg,[dirLoc dirFig figname '.jpg'],'jpg');
    fg = fContour2numVL;
    figname = ['DiffM1FRcontour_cmprOscAmpKO' tmpTxt];
    saveas(fg,[dirLoc dirFig figname '.fig'],'fig'); saveas(fg,[dirLoc dirFig figname '.jpg'],'jpg');
end
%%
% end
