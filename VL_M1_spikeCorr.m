%%   The VL-M1 spike correlation
RUN_KO = 0;
SAVE_FIG = 0;
%%
cuttime = 500;
p3_ii =1 ; %trial no
Collect_sumCC_WT = cell(length(PARAM1),length(PARAM2),length(PARAM4),length(PARAM5));
Collect_peakDel_WT = zeros(length(PARAM1),length(PARAM2),length(PARAM4), length(PARAM5));
if(RUN_KO)
    Collect_sumCC_KO = cell(length(PARAM1),length(PARAM2),length(PARAM4), length(PARAM5));
    Collect_peakDel_KO = zeros(length(PARAM1),length(PARAM2),length(PARAM4), length(PARAM5));
end
RR = length(PARAM4); CC = length(PARAM5); cnt_aut = 0;

Trange_LST = [25 50];
% for t_ii = 1 : length(Trange_LST)
%     Trange = Trange_LST(t_ii);

tt_run = tic();
for p4_ii = 1 : length(PARAM4) %osc f
    Trange = 0.5*1000/PARAM4(p4_ii);
    for p5_ii = 1 : length(PARAM5)
        
        tosc = tic();
        figPos = [    399         105        1712        1245];
        nR = length(PARAM1); nC = length(PARAM2);
        fp1_WT= figure;  set(fp1_WT, 'position', figPos); set(fp1_WT,'PaperPositionMode','auto');
        subplot(nR,nC,1); suptitle({CtypeTxt,['WT : '  get_Parameters_titleText(PARAMETERS,[4,5],[p4_ii, p5_ii])  ', T_r_a_n_g_e = ' num2str(Trange)]});
        if(RUN_KO)
            fp1_KO= figure;  set(fp1_KO, 'position',figPos); set(fp1_KO,'PaperPositionMode','auto');
            subplot(nR,nC,1); suptitle({CtypeTxt , ['KO : ' get_Parameters_titleText(PARAMETERS,[4,5],[p4_ii, p5_ii])   ', T_r_a_n_g_e = ' num2str(Trange)]});
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
                [box,sumAC, cntSample] = CrossCorrFromSpktrain(spkBin_M1,spkBin_VL, Trange, cuttime );
                avgsumAC = sum(sumAC(:)) / length(box);
                sumAC = sumAC - avgsumAC;
                Collect_sumCC_WT{p1_ii,p2_ii,p4_ii,p5_ii} = sumAC;
                [pks,locs]=findpeaks(sumAC);
                if ~isempty(locs)
                    Collect_peakDel_WT(p1_ii,p2_ii,p4_ii,p5_ii) = min(box(locs) + 1i) -1i;
                end
                figure(fp1_WT); subplot(nR,nC,cnt);
                bar(box,sumAC,'r'); title({tmptt,['closest pks = ' num2str(Collect_peakDel_WT(p1_ii,p2_ii,p4_ii,p5_ii))]});
                
                toc(tt);
                if (RUN_KO)
                    
                    tt = tic();
                    spkBin_VL = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii,p5_ii}.VL.KO.All.spktrain;
                    spkBin_M1 = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii,p5_ii}.M1.KO.All.spktrain;
                    [box,sumAC, cntSample] = CrossCorrFromSpktrain(spkBin_M1,spkBin_VL, Trange, cuttime );
                    avgsumAC = sum(sumAC(:)) / length(box);
                    sumAC = sumAC - avgsumAC;
                    Collect_sumCC_KO{p1_ii,p2_ii,p4_ii,p5_ii} = sumAC;
                    [pks,locs]=findpeaks(sumAC);
                    if ~isempty(locs)
                        Collect_peakDel_WT(p1_ii,p2_ii,p4_ii,p5_ii) = min(box(locs) + 1i) -1i;
                    end
                    figure(fp1_KO); subplot(nR,nC,cnt);
                    bar(box,sumAC,'r'); title({tmptt,['closest pks = ' num2str(Collect_peakDel_WT(p1_ii,p2_ii,p4_ii,p5_ii))]});
                    
                    toc(tt);
                end
                
            end
        end
        
        %%
        set(fp1_WT, 'position', figPos); set(fp1_WT,'PaperPositionMode','auto');
        suptitle({CtypeTxt,['WT : '  get_Parameters_titleText(PARAMETERS,[4,5],[p4_ii, p5_ii])  ', T_r_a_n_g_e = ' num2str(Trange)]});
        if(RUN_KO)
           set(fp1_KO, 'position',figPos); set(fp1_KO,'PaperPositionMode','auto');
           suptitle({CtypeTxt , ['KO : ' get_Parameters_titleText(PARAMETERS,[4,5],[p4_ii, p5_ii])   ', T_r_a_n_g_e = ' num2str(Trange)]});
        end

        
        
        
        
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

%%
if(0)
    Backup_Collect_sumCC_WT = Collect_sumCC_WT;
    if(RUN_KO)
        Backup_Collect_sumCC_KO = Collect_sumCC_KO;
    end
end
