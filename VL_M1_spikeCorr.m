%%   The VL-M1 spike correlation

%
% Collect_sumAC_WT = cell(length(PARAM1),length(PARAM2));
% Collect_sumAC_KO = cell(length(PARAM1),length(PARAM2));
% Collect_sumCC_WT = cell(length(PARAM1),length(PARAM2),length(PARAM3));
% Collect_sumCC_KO = cell(length(PARAM1),length(PARAM2),length(PARAM3));
% Trange_LST = [5 10 25 30];
% for t_ii = 1 : length(Trange_LST)
%     Trange = Trange_LST(t_ii);
%     for p1_ii = 1 : length(PARAM1)
%         fp1_WT= figure;  set(fp1_WT, 'position',[ 1          41        1920         1000]); set(fp1_WT,'PaperPositionMode','auto');
%         subplot(length(PARAM3)+1,length(PARAM2),1); suptitle(['WT : ' saveTxt1 '=' num2str(PARAM1(p1_ii)) ', T_r_a_n_g_e = ' num2str(Trange)]);
%         fp1_KO= figure;  set(fp1_KO, 'position',[  1          41        1920        1000]); set(fp1_KO,'PaperPositionMode','auto');
%         subplot(length(PARAM3)+1,length(PARAM2),1); suptitle(['KO : ' saveTxt1 '=' num2str(PARAM1(p1_ii)) ', T_r_a_n_g_e = ' num2str(Trange)]);
%         cnt = 0;
%         for p3_i = 0 : length(PARAM3)
%             for p2_ii = 1 : length(PARAM2)
%                 if (p3_i == 0)
%                     % Auto Corr of VL
%                     p3_ii = 1;
%                     cnt = cnt+1;
%                     spkBin = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii}.WT.VL.All.spktrain;
%                     [box,sumAC, cntSample] = AutocorrFromSpktrain(spkBin, Trange, cuttime );
%                     figure(fp1_WT); subplot(length(PARAM3)+1,length(PARAM2),cnt);
%                     bar(box,sumAC); title([saveTxt2 ' ' num2str(PARAM2(p2_ii)) ]);
%                     Collect_sumAC_WT{p1_ii,p2_ii} = sumAC;
%
%                     spkBin = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii}.KO.VL.All.spktrain;
%                     [box,sumAC, cntSample] = AutocorrFromSpktrain(spkBin, Trange, cuttime );
%                     figure(fp1_KO); subplot(length(PARAM3)+1,length(PARAM2),cnt);
%                     bar(box,sumAC,'r'); title([saveTxt2 ' ' num2str(PARAM2(p2_ii)) ]);
%                     Collect_sumAC_KO{p1_ii,p2_ii} = sumAC;
%                     %                 k = waitforbuttonpress;
%                 else
%                     % Cross Corr of VL and M1
%                     p3_ii = p3_i;
%                     cnt = cnt+1;
%
%                     spkBin_VL = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii}.WT.VL.All.spktrain;
%                     spkBin_M1 = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii}.WT.M1.All.spktrain;
%                     [box,sumAC, cntSample] = CrossCorrFromSpktrain(spkBin_M1,spkBin_VL, Trange, cuttime )
%                     figure(fp1_WT); subplot(length(PARAM3)+1,length(PARAM2),cnt);
%                     bar(box,sumAC); title([saveTxt3 ' ' num2str(PARAM3(p3_ii)) ]);
%                     Collect_sumCC_WT{p1_ii,p2_ii,p3_ii} = sumAC;
%
%                     spkBin_VL = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii}.KO.VL.All.spktrain;
%                     spkBin_M1 = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii}.KO.M1.All.spktrain;
%                     [box,sumAC, cntSample] = CrossCorrFromSpktrain(spkBin_M1,spkBin_VL, Trange, cuttime )
%                     figure(fp1_KO); subplot(length(PARAM3)+1,length(PARAM2),cnt);
%                     bar(box,sumAC,'r'); title([saveTxt3 ' ' num2str(PARAM3(p3_ii)) ]);
%                     Collect_sumCC_KO{p1_ii,p2_ii,p3_ii} = sumAC;
%                     %                 k = waitforbuttonpress;
%                 end
%
%             end
%         end
%
%         if (SAVE_FIG)
%             ffig = [ dirLoc dirFig 'CorrrelationM1VL_Trange' num2str(Trange) '_' saveTxt1 num2str(PARAM1(p1_ii))];
%
%             saveas( fp1_WT, [ffig '_WT.jpg'], 'jpg')
%             saveas( fp1_WT, [ffig '_WT.fig'], 'fig')
%             saveas( fp1_KO, [ffig '_KO.jpg'], 'jpg')
%             saveas( fp1_KO, [ffig '_KO.fig'], 'fig')
%
%         end
%     end
% end
%
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
%             spkBin = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii}.WT.VL.All.spktrain;
%             [box,sumAC, cntSample] = AutocorrFromSpktrain(spkBin, Trange, cuttime );
%             subplot(121); hold on;
%             plot(box,sumAC,'Color',[0 0 ccc(s_ii)]); LEG{s_ii} = [saveTxt2 ' ' num2str(PARAM2(p2_ii)) ];
%             %KO
%             spkBin = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii}.KO.VL.All.spktrain;
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