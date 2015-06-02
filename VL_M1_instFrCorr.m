%%   The VL-M1 spike correlation
       RUN_KO = 0;
    
%% 
p3_ii =1 ; %trial no
SAVE_FIG = 0;
% Collect_sumAC_WT = cell(length(PARAM4),length(PARAM5));
% Collect_sumAC_KO = cell(length(PARAM4),length(PARAM5));
Collect_corrcoef_WT = cell(length(PARAM1),length(PARAM2),length(PARAM4),length(PARAM5));
Collect_corrcoef_KO  = cell(length(PARAM1),length(PARAM2),length(PARAM4), length(PARAM5));
 figPos = [ 1          41        1920         1000];
fg_autoCor_VL_WT = figure;  set(gcf, 'position', figPos); set(gcf,'PaperPositionMode','auto');
fg_autoCor_VL_KO = figure;  set(fg_autoCor_VL_KO, 'position', figPos); set(fg_autoCor_VL_KO,'PaperPositionMode','auto');
RR = length(PARAM4); CC = length(PARAM5); cnt_aut = 0;

tt_run = tic();
for p4_ii = 1 : length(PARAM4)
    for p5_ii = 1 : length(PARAM5)
        tosc = tic();
        figPos = [ 1          41        1920         1000];
        nR = length(PARAM1); nC = length(PARAM2); 
        fp1_WT= figure;  set(fp1_WT, 'position', figPos); set(fp1_WT,'PaperPositionMode','auto');
        subplot(nR,nC,1); suptitle({CtypeTxt,['WT : '  get_Parameters_titleText(PARAMETERS,[4,5],[p4_ii, p5_ii]) ]});
       if(RUN_KO)
        fp1_KO= figure;  set(fp1_KO, 'position',figPos); set(fp1_KO,'PaperPositionMode','auto');
        subplot(nR,nC,1); suptitle({CtypeTxt , ['KO : ' get_Parameters_titleText(PARAMETERS,[4,5],[p4_ii, p5_ii]) });
       end
        cnt = 0;

        for p1_ii = 1 : length(PARAM1)
            for p2_ii = 1 : length(PARAM2)
                
                    % Cross Corr of VL and M1
                    
                    cnt = cnt+1;
                    tmptt = get_Parameters_titleText(PARAMETERS, [1,2], [p1_ii p2_ii]);
                    tt = tic();
					
                    spkBin_VL = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii,p5_ii}.VL.WT.All.spktrain;
                    spkBin_M1 = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii,p5_ii}.M1.WT.All.spktrain;
					
					instAvgFR_VL = mean(spkBin_VL, 1)*1000;
					instAvgFR_M1 = mean(spkBin_M1, 1)*1000;
					
					gsig = 3;   % filter [ms]
					gfil = exp( -[-round(4*gsig):round(4*gsig)].^2 /2/gsig^2);
					gfil = gfil/sum(gfil(:));
        
					smoothAvgFR_VL = conv(instAvgFR_VL , gfil, 'same');
					smoothAvgFR_M1 = conv(instAvgFR_M1 , gfil, 'same');
					t_axis = 0:length(instAvgFR_VL)-1;
					

                    figure(fp1_WT); subplot(nR,nC,cnt);
					plot(t_Axis, smoothAvgFR_VL); hold on; plot(t_Axis, smoothAvgFR_M1); hold on; 
                    Collect_corrcoef_WT{p1_ii,p2_ii,p4_ii,p5_ii} = corrcoef(smoothAvgFR_VL , smoothAvgFR_M1);
                    toc(tt);
                    
					if (RUN_KO)
                    tt = tic();
                    spkBin_VL = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii,p5_ii}.VL.KO.All.spktrain;
                    spkBin_M1 = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii,p5_ii}.M1.KO.All.spktrain;	
					instAvgFR_VL = mean(spkBin_VL, 1)*1000;
					instAvgFR_M1 = mean(spkBin_M1, 1)*1000;
        
					smoothAvgFR_VL = conv(instAvgFR_VL , gfil, 'same');
					smoothAvgFR_M1 = conv(instAvgFR_M1 , gfil, 'same');
					t_axis = 0:length(instAvgFR_VL)-1;
					
                    figure(fp1_KO); subplot(nR,nC,cnt);
					plot(t_Axis, smoothAvgFR_VL); hold on; plot(t_Axis, smoothAvgFR_M1); hold on; 
                    Collect_corrcoef_KO{p1_ii,p2_ii,p4_ii,p5_ii} = corrcoef(smoothAvgFR_VL , smoothAvgFR_M1);
                    toc(tt);    
                    end

            end
        end

        %% 

        if (SAVE_FIG)
            ffig = [ dirLoc dirFig 'InstVL_M1_' get_Parameters_saveText(PARAMETERS,[4,5],[p4_ii, p5_ii]) codeTxt];

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