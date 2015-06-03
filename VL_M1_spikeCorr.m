%%   The VL-M1 spike correlation
       RUN_KO = 0;
    
%% 
cuttime = 500;
p3_ii =1 ; %trial no
SAVE_FIG = 0;
% Collect_sumAC_WT = cell(length(PARAM4),length(PARAM5));
% Collect_sumAC_KO = cell(length(PARAM4),length(PARAM5));
Collect_sumCC_WT = cell(length(PARAM1),length(PARAM2),length(PARAM4),length(PARAM5));
Collect_sumCC_KO = cell(length(PARAM1),length(PARAM2),length(PARAM4), length(PARAM5));
 figPos = [ 1          41        1920         1000];
fg_autoCor_VL_WT = figure;  set(gcf, 'position', figPos); set(gcf,'PaperPositionMode','auto');
fg_autoCor_VL_KO = figure;  set(fg_autoCor_VL_KO, 'position', figPos); set(fg_autoCor_VL_KO,'PaperPositionMode','auto');

RR = length(PARAM4); CC = length(PARAM5); cnt_aut = 0;
Trange_LST = [25 50];
% for t_ii = 1 : length(Trange_LST)
%     Trange = Trange_LST(t_ii);

tt_run = tic();
for p4_ii = 1 : length(PARAM4) %osc f
Trange = 1000/PARAM4(p4_ii); 
    for p5_ii = 1 : length(PARAM5)
	
        tosc = tic();
        figPos = [ 1          41        1920         1000];
        nR = length(PARAM1); nC = length(PARAM2);
        fp1_WT= figure;  set(fp1_WT, 'position', figPos); set(fp1_WT,'PaperPositionMode','auto');
        subplot(nR,nC,1); suptitle({CtypeTxt,['WT : '  get_Parameters_titleText(PARAMETERS,[4,5],[p4_ii, p5_ii])  ', T_r_a_n_g_e = ' num2str(Trange)]});
       if(RUN_KO)
        fp1_KO= figure;  set(fp1_KO, 'position',figPos); set(fp1_KO,'PaperPositionMode','auto');
        subplot(nR,nC,1); suptitle({CtypeTxt , ['KO : ' get_Parameters_titleText(PARAMETERS,[4,5],[p4_ii, p5_ii])   ', T_r_a_n_g_e = ' num2str(Trange)]});
       end
        cnt = 0;
%         % Auto Corr at VL
%         cnt_aut = cnt_aut + 1;
%         
%         p1_ii = 1; p2_ii = 1; 
%          tmptt = get_Parameters_titleText(PARAMETERS, [4 5], [p4_ii p5_ii]);
%                     tt = tic();
%                     spkBin = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii,p5_ii}.VL.WT.All.spktrain;
%                     [box,sumAC, cntSample] = AutocorrFromSpktrain(spkBin, Trange, cuttime );
%                     figure(fg_autoCor_VL_WT); subplot(RR,CC,cnt_aut);
%                     bar(box,sumAC); title(tmptt);
%                     Collect_sumAC_WT{p4_ii,p5_ii} = sumAC;
%                     toc(tt);
%                     tt = tic();
%                     spkBin = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii,p5_ii}.VL.KO.All.spktrain;
%                     [box,sumAC, cntSample] = AutocorrFromSpktrain(spkBin, Trange, cuttime );
%                  figure(fg_autoCor_VL_KO); subplot(RR,CC,cnt_aut);
%                     bar(box,sumAC,'r'); title(tmptt);
%                     Collect_sumAC_KO{p4_ii,p5_ii} = sumAC;
%                     toc(tt);
                    
        for p1_ii = 1 : length(PARAM1)
            for p2_ii = 1 : length(PARAM2)
                
                    % Cross Corr of VL and M1
                    
                    cnt = cnt+1;
                        tmptt = get_Parameters_titleText(PARAMETERS, [1,2], [p1_ii p2_ii]);
                        tt = tic();
                    spkBin_VL = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii,p5_ii}.VL.WT.All.spktrain;
                    spkBin_M1 = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii,p5_ii}.M1.WT.All.spktrain;
                    [box,sumAC, cntSample] = CrossCorrFromSpktrain(spkBin_M1,spkBin_VL, Trange, cuttime );
                    avgsumAC = sum(sumAC(:)) / length(box); 
					sumAC = sumAC - avgsumAC;
					figure(fp1_WT); subplot(nR,nC,cnt);
                    bar(box,sumAC); title(tmptt);
                    Collect_sumCC_WT{p1_ii,p2_ii,p4_ii,p5_ii} = sumAC;
                    toc(tt);
                    if (RUN_KO)
                    tt = tic();
                    spkBin_VL = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii,p5_ii}.VL.KO.All.spktrain;
                    spkBin_M1 = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii,p5_ii}.M1.KO.All.spktrain;
                    [box,sumAC, cntSample] = CrossCorrFromSpktrain(spkBin_M1,spkBin_VL, Trange, cuttime );
                    avgsumAC = sum(sumAC(:)) / length(box); 
					sumAC = sumAC - avgsumAC;
					figure(fp1_KO); subplot(nR,nC,cnt);
                    bar(box,sumAC,'r'); title(tmptt);
                    Collect_sumCC_KO{p1_ii,p2_ii,p4_ii,p5_ii} = sumAC;
                    toc(tt);    
                    end

            end
        end

        %% 

        if (SAVE_FIG)
            ffig = [ dirLoc dirFig 'CorrrelationM1VL_Trange' num2str(Trange) '_' get_Parameters_saveText(PARAMETERS,[4,5],[p4_ii, p5_ii]) codeTxt];

            saveas( fp1_WT, [ffig '_WT.jpg'], 'jpg')
            saveas( fp1_WT, [ffig '_WT.fig'], 'fig')
            if(RUN_KO)
            saveas( fp1_KO, [ffig '_KO.jpg'], 'jpg')
            saveas( fp1_KO, [ffig '_KO.fig'], 'fig')
            end

        end
        disp('Run time for one Osc case ');
        toc(tosc);
    end
end
disp('RunTime for all')
toc(tt_run);
  
  
% if (SAVE_FIG)
%             ffig = [ dirLoc dirFig 'AutoCorrrVL_Trange' num2str(Trange) '_' get_Parameters_saveText(PARAMETERS,[4,5],[p4_ii, p5_ii]) codeTxt];
% 
%             saveas( fg_autoCor_VL_WT, [ffig '_WT.jpg'], 'jpg');             saveas( fg_autoCor_VL_WT, [ffig '_WT.fig'], 'fig');
%             saveas(fg_autoCor_VL_KO, [ffig '_KO.jpg'], 'jpg');            saveas( fg_autoCor_VL_KO, [ffig '_KO.fig'], 'fig');
% end

%
% %%  Auto Corr with different level of synchronization
%
%
% % PARAM1 = InputFR_LST;
% % legTxt1 = 'Input Frequency';
% % saveTxt1 = 'InputFR';
% % titleTxt1 = 'Input Frequency';
% % PARAM2 = SynchLvl_LST;
% % legTxt2 = 'Synchronization Level[Range = 0-1]';
% % saveTxt2 = 'SynchLvl';
% % titleTxt2 = 'Synchronization Level';
% % PARAM3 = W_VL_M1_LST;
% % legTxt3 = 'VL-M1 connection weight';
% % titleTxt3 = 'VL-M1 weight';
% % saveTxt3 = 'wVLM1';
%
% Trange_LST = [5 10 25 30];
% for f_ii = 1 : length(InputFR_LST)
%     p1_ii = f_ii; p3_ii = 1;
%     for t_ii = 1 : length(Trange_LST)
%         Trange = Trange_LST(t_ii);
%
%
%         LEG = cell(length(SynchLvl_LST),1);
%         tmpf = figure; set(tmpf,'position',[279         447        1377         420]);  set(tmpf,'PaperPositionMode','auto');
%         ccc = [0:1/length(SynchLvl_LST):1];
%
%         for s_ii = 1 : length(SynchLvl_LST)
%             p2_ii = s_ii;
%             %WT
%             spkBin = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii}.VL.WT.All.spktrain;
%             [box,sumAC, cntSample] = AutocorrFromSpktrain(spkBin, Trange, cuttime );
%             subplot(121); hold on;
%             plot(box,sumAC,'Color',[0 0 ccc(s_ii)]); LEG{s_ii} = [saveTxt2 ' ' num2str(PARAM2(p2_ii)) ];
%             %KO
%             spkBin = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii}.VL.KO.All.spktrain;
%             [box,sumAC, cntSample] = AutocorrFromSpktrain(spkBin, Trange, cuttime );
%             subplot(122); hold on;
%             plot(box,sumAC,'Color',[ccc(s_ii) 0 0]); LEG{s_ii} = [saveTxt2 ' ' num2str(PARAM2(p2_ii)) ];
%         end
%         subplot(121); legend(LEG); title('WT'); subplot(122); legend(LEG); title('KO')
%         suptitle([titleTxt1 ' = ' num2str(InputFR_LST(f_ii)) ' , '  'T_r_a_n_g_e = ' num2str(Trange)]);
%          if (SAVE_FIG)
%             ffig = [ dirLoc dirFig 'CombinedSynchLevel_Trange' num2str(Trange) '_' saveTxt1 num2str(PARAM1(p1_ii))];
%             saveas(  tmpf, [ffig '.jpg'], 'jpg')
%             saveas(  tmpf, [ffig '.fig'], 'fig')
%
%         end
%     end
% end