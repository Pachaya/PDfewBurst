%%   The VL spike correlation
       
    
%% 
RUN_KO = 0;
p3_ii =1 ; %trial no
SAVE_FIG = 0;
Collect_sumAC_WT = cell(length(PARAM4),length(PARAM5));
Collect_sumAC_KO = cell(length(PARAM4),length(PARAM5));
% Collect_sumCC_WT = cell(length(PARAM1),length(PARAM2),length(PARAM4),length(PARAM5));
% Collect_sumCC_KO = cell(length(PARAM1),length(PARAM2),length(PARAM4), length(PARAM5));
 figPos = [ 1          41        1920         1000];
fg_autoCor_VL_WT = figure;  set(gcf, 'position', figPos); set(gcf,'PaperPositionMode','auto');
fg_autoCor_VL_KO = figure;  set(fg_autoCor_VL_KO, 'position', figPos); set(fg_autoCor_VL_KO,'PaperPositionMode','auto');

RR = length(PARAM4); CC = length(PARAM5); cnt_aut = 0;
Trange_LST = [25 50];
% for t_ii = 1 : length(Trange_LST)
%     Trange = Trange_LST(t_ii);
Trange = 50;
tt_run = tic();
for p4_ii = 1 : length(PARAM4)
    for p5_ii = 1 : length(PARAM5)
%         figPos = [ 1          41        1920         1000];
%         nR = length(PARAM1); nC = length(PARAM2);
%         fp1_WT= figure;  set(fp1_WT, 'position', figPos); set(fp1_WT,'PaperPositionMode','auto');
%         subplot(nR,nC,1); suptitle(['WT : '  get_Parameters_titleText(PARAMETERS,[4,5],[p4_ii, p5_ii])  ', Connection Type :' CtypeTxt   ', T_r_a_n_g_e = ' num2str(Trange)]);
%         fp1_KO= figure;  set(fp1_KO, 'position',figPos); set(fp1_KO,'PaperPositionMode','auto');
%         subplot(nR,nC,1); suptitle(['KO : ' get_Parameters_titleText(PARAMETERS,[4,5],[p4_ii, p5_ii])  ', Connection Type :' CtypeTxt   ', T_r_a_n_g_e = ' num2str(Trange)]);
%         cnt = 0;

%         % Auto Corr at VL
        cnt_aut = cnt_aut + 1;
        
        p1_ii = 1; p2_ii = 1; 
         tmptt = get_Parameters_titleText(PARAMETERS, [4 5], [p4_ii p5_ii]);
                    tt = tic();
                    spkBin = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii,p5_ii}.VL.WT.All.spktrain;
                    [box,sumAC, cntSample] = AutocorrFromSpktrain(spkBin, Trange, cuttime );
                    figure(fg_autoCor_VL_WT); subplot(RR,CC,cnt_aut);
                    bar(box,sumAC); title(tmptt);
                    Collect_sumAC_WT{p4_ii,p5_ii} = sumAC;
                    toc(tt);
                    if(RUN_KO)
                    tt = tic();
                    spkBin = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii,p5_ii}.VL.KO.All.spktrain;
                    [box,sumAC, cntSample] = AutocorrFromSpktrain(spkBin, Trange, cuttime );
                 figure(fg_autoCor_VL_KO); subplot(RR,CC,cnt_aut);
                    bar(box,sumAC,'r'); title(tmptt);
                    Collect_sumAC_KO{p4_ii,p5_ii} = sumAC;
                    toc(tt);
                    end
                    
       
    end
end
figure(fg_autoCor_VL_WT); suptitle({CtypeTxt,['WT : ' 'T_r_a_n_g_e = ' num2str(Trange)]});
if(RUN_KO)
figure(fg_autoCor_VL_KO); suptitle({CtypeTxt,['KO : ' 'T_r_a_n_g_e = ' num2str(Trange)]});
end
 toc(tt_run);
if (SAVE_FIG)
            ffig = [ dirLoc dirFig 'AutoCorrrVL_Trange' num2str(Trange) '_' codeTxt];

            saveas( fg_autoCor_VL_WT, [ffig '_WT.jpg'], 'jpg');             saveas( fg_autoCor_VL_WT, [ffig '_WT.fig'], 'fig');
            if(RUN_KO)
            saveas(fg_autoCor_VL_KO, [ffig '_KO.jpg'], 'jpg');            saveas( fg_autoCor_VL_KO, [ffig '_KO.fig'], 'fig');
            end
end
