% Quantitatively Test Synchronization Level
close all
clear all

dirLoc = 'TestSynchLevel/';
dirFig = 'Fig150319/';
simCode  = 'test';
CoreName = '';

SAVE_FIG =1;

Nsample = 20;
SynchLvl = 0.5;
Tcut = 200;
Tmeasure = 1000;
TSTOP = Tcut + Tmeasure;
tstop = TSTOP;
InputFR = 40;
W_SPK = 0.005;
W_VL_M1 = 0.002/Nsample ;


InputFR_LST = [40; 100; 250;];
SynchLvl_LST = [0:0.2:1];%[0:1/Nsample:1];
W_VL_M1_LST = [0.0001 : 0.0002 : 0.0001*Nsample]; %[0.0001 : 0.0001 : 0.0001*Nsample];

PARAM1 = InputFR_LST;
legTxt1 = 'Input Frequency';
saveTxt1 = 'InputFR';
titleTxt1 = 'Input Frequency';
PARAM2 = SynchLvl_LST;
legTxt2 = 'Synchronization Level[Range = 0-1]';
saveTxt2 = 'SynchLvl';
titleTxt2 = 'Synchronization Level';
PARAM3 = W_VL_M1_LST;
legTxt3 = 'VL-M1 connection weight';
titleTxt3 = 'VL-M1 weight';
saveTxt3 = 'wVLM1';

ACT_Record = cell(length(PARAM1),length(PARAM2),length(PARAM3));
coreName = '';


for p1_ii = 1 : length(PARAM1)
    for p2_ii = 1 : length(PARAM2)
        for p3_ii = 1 : length(PARAM3)
            inF_ii = p1_ii; sl_ii = p2_ii; wVM_ii = p3_ii;
            
            W_VL_M1= W_VL_M1_LST(wVM_ii);
            SynchLvl = SynchLvl_LST(sl_ii);
            InputFR = InputFR_LST(inF_ii);
            CUTTIME = 200;
            
            simCode =[ coreName 'Nsample' num2str(Nsample) '_SynchLvl' num2str(SynchLvl) '_TSTOP' num2str(TSTOP) '_InputFR' num2str(InputFR) ...
                '_wSPK' num2str(W_SPK) '_wVLM1_' num2str(W_VL_M1)];
            
            %WT_VL
            Name_postfix  = ['WT_VL_' simCode ];
            VL_BaselineAct
            tmp.E = E;      tmp.I = I;  tmp.All = All; tmp.All.somaVall = somaVall;
            WT.VL = tmp;
            
            %KO_VL
            Name_postfix  = ['KO_VL_' simCode ];
            VL_BaselineAct
            tmp.E = E;      tmp.I = I;  tmp.All = All; tmp.All.somaVall = somaVall;
            KO.VL = tmp;
            
            %WT_M1
            Name_postfix  = ['WT_M1_' simCode ];
            VL_BaselineAct
            tmp.E = E;      tmp.I = I;  tmp.All = All; tmp.All.somaVall = somaVall;
            WT.M1 = tmp;
            
            %KO_M1
            Name_postfix  = ['KO_M1_' simCode ];
            VL_BaselineAct
            tmp.E = E;      tmp.I = I;  tmp.All = All; tmp.All.somaVall = somaVall;
            KO.M1 = tmp;
            
            ACT.WT =  WT; ACT.KO = KO;
            
            ACT_Record{p1_ii,p2_ii,p3_ii} = ACT;
        end
    end
end


%%

for p1_ii = 1 : length(PARAM1)
    for p2_ii = 1 : length(PARAM2)
        for p3_ii = 1 : length(PARAM3)
            
            simTxt = [titleTxt1 ' ' num2str(PARAM1(p1_ii)) ', ' titleTxt2 ' ' num2str(PARAM2(p2_ii)) ', ' titleTxt3 ' ' num2str(PARAM3(p3_ii))];
            %raster
            spkBin = ACT_Record{p1_ii,p2_ii,p3_ii}.WT.VL.All.spktrain;
            [freq_all_WT,spkTime_all_WT, fg_handle_WT ] = raster_from_spkbin( spkBin,cuttime, Tstop, [ simTxt ' : WT']);
            set(fg_handle_WT, 'position',[  449   450   791   528])
            spkBin = ACT_Record{p1_ii,p2_ii,p3_ii}.KO.VL.All.spktrain;
            [freq_all_KO,spkTime_all_KO, fg_handle_KO ] = raster_from_spkbin( spkBin,cuttime, Tstop, [ simTxt ' : KO']);
            set(fg_handle_KO, 'position',[  449   450   791   528])
            disp(['########### Trial = ' simTxt ])
            disp(['### WT avg freq = ' num2str(mean(freq_all_WT))])
            disp(['### KO avg freq = ' num2str(mean(freq_all_KO))])
            if (SAVE_FIG)
                ffig = [ dirLoc dirFig 'RasterPlot_' saveTxt1 num2str(PARAM1(p1_ii)) '_' saveTxt2 num2str(PARAM2(p2_ii)) '_' saveTxt3 num2str(PARAM3(p3_ii))];
                %         ffig = [ dirLoc dirFig 'RasterPlot_Trial' num2str(TRIAL_NO)];
                saveas( fg_handle_WT, [ffig '_WT.jpg'], 'jpg')
                saveas( fg_handle_WT, [ffig '_WT.fig'], 'fig')
                saveas( fg_handle_KO, [ffig '_KO.jpg'], 'jpg')
                saveas( fg_handle_KO, [ffig '_KO.fig'], 'fig')
                
            end
            
        end
    end
end
%%
SAVE_NEURON_ACT =0;
coreFileName ='';
if (SAVE_NEURON_ACT)
    save([dirLoc dirFig 'Result_ACT_Record_' coreFileName '.mat'], 'ACT_Record','-v7.3');
end
%%

Collect_sumAC_WT = cell(length(PARAM1),length(PARAM2));
Collect_sumAC_KO = cell(length(PARAM1),length(PARAM2));
Collect_sumCC_WT = cell(length(PARAM1),length(PARAM2),length(PARAM3));
Collect_sumCC_KO = cell(length(PARAM1),length(PARAM2),length(PARAM3));
Trange_LST = [5 10 25 30];
for t_ii = 1 : length(Trange_LST)
    Trange = Trange_LST(t_ii);
    for p1_ii = 1 : length(PARAM1)
        fp1_WT= figure;  set(fp1_WT, 'position',[ 1          41        1920         1000]); set(fp1_WT,'PaperPositionMode','auto');
        subplot(length(PARAM3)+1,length(PARAM2),1); suptitle(['WT : ' saveTxt1 '=' num2str(PARAM1(p1_ii)) ', T_r_a_n_g_e = ' num2str(Trange)]);
        fp1_KO= figure;  set(fp1_KO, 'position',[  1          41        1920        1000]); set(fp1_KO,'PaperPositionMode','auto');
        subplot(length(PARAM3)+1,length(PARAM2),1); suptitle(['KO : ' saveTxt1 '=' num2str(PARAM1(p1_ii)) ', T_r_a_n_g_e = ' num2str(Trange)]);
        cnt = 0;
        for p3_i = 0 : length(PARAM3)
            for p2_ii = 1 : length(PARAM2)
                if (p3_i == 0)
                    % Auto Corr of VL
                    p3_ii = 1;
                    cnt = cnt+1;
                    spkBin = ACT_Record{p1_ii,p2_ii,p3_ii}.WT.VL.All.spktrain;
                    [box,sumAC, cntSample] = AutocorrFromSpktrain(spkBin, Trange, cuttime );
                    figure(fp1_WT); subplot(length(PARAM3)+1,length(PARAM2),cnt);
                    bar(box,sumAC); title([saveTxt2 ' ' num2str(PARAM2(p2_ii)) ]);
                    Collect_sumAC_WT{p1_ii,p2_ii} = sumAC;
                    
                    spkBin = ACT_Record{p1_ii,p2_ii,p3_ii}.KO.VL.All.spktrain;
                    [box,sumAC, cntSample] = AutocorrFromSpktrain(spkBin, Trange, cuttime );
                    figure(fp1_KO); subplot(length(PARAM3)+1,length(PARAM2),cnt);
                    bar(box,sumAC,'r'); title([saveTxt2 ' ' num2str(PARAM2(p2_ii)) ]);
                    Collect_sumAC_KO{p1_ii,p2_ii} = sumAC;
                    %                 k = waitforbuttonpress;
                else
                    % Cross Corr of VL and M1
                    p3_ii = p3_i;
                    cnt = cnt+1;
                    
                    spkBin_VL = ACT_Record{p1_ii,p2_ii,p3_ii}.WT.VL.All.spktrain;
                    spkBin_M1 = ACT_Record{p1_ii,p2_ii,p3_ii}.WT.M1.All.spktrain;
                    [box,sumAC, cntSample] = CrossCorrFromSpktrain(spkBin_M1,spkBin_VL, Trange, cuttime )
                    figure(fp1_WT); subplot(length(PARAM3)+1,length(PARAM2),cnt);
                    bar(box,sumAC); title([saveTxt3 ' ' num2str(PARAM3(p3_ii)) ]);
                    Collect_sumCC_WT{p1_ii,p2_ii,p3_ii} = sumAC;
                    
                    spkBin_VL = ACT_Record{p1_ii,p2_ii,p3_ii}.KO.VL.All.spktrain;
                    spkBin_M1 = ACT_Record{p1_ii,p2_ii,p3_ii}.KO.M1.All.spktrain;
                    [box,sumAC, cntSample] = CrossCorrFromSpktrain(spkBin_M1,spkBin_VL, Trange, cuttime )
                    figure(fp1_KO); subplot(length(PARAM3)+1,length(PARAM2),cnt);
                    bar(box,sumAC,'r'); title([saveTxt3 ' ' num2str(PARAM3(p3_ii)) ]);
                    Collect_sumCC_KO{p1_ii,p2_ii,p3_ii} = sumAC;
                    %                 k = waitforbuttonpress;
                end
                
            end
        end
        
        if (SAVE_FIG)
            ffig = [ dirLoc dirFig 'CorrrelationM1VL_Trange' num2str(Trange) '_' saveTxt1 num2str(PARAM1(p1_ii))];
            
            saveas( fp1_WT, [ffig '_WT.jpg'], 'jpg')
            saveas( fp1_WT, [ffig '_WT.fig'], 'fig')
            saveas( fp1_KO, [ffig '_KO.jpg'], 'jpg')
            saveas( fp1_KO, [ffig '_KO.fig'], 'fig')
            
        end
    end
end


%%  Auto Corr with different level of synchronization


% PARAM1 = InputFR_LST;
% legTxt1 = 'Input Frequency';
% saveTxt1 = 'InputFR';
% titleTxt1 = 'Input Frequency';
% PARAM2 = SynchLvl_LST;
% legTxt2 = 'Synchronization Level[Range = 0-1]';
% saveTxt2 = 'SynchLvl';
% titleTxt2 = 'Synchronization Level';
% PARAM3 = W_VL_M1_LST;
% legTxt3 = 'VL-M1 connection weight';
% titleTxt3 = 'VL-M1 weight';
% saveTxt3 = 'wVLM1';

Trange_LST = [5 10 25 30];
for f_ii = 1 : length(InputFR_LST)
    p1_ii = f_ii; p3_ii = 1;
    for t_ii = 1 : length(Trange_LST)
        Trange = Trange_LST(t_ii);
        
        
        LEG = cell(length(SynchLvl_LST),1);
        tmpf = figure; set(tmpf,'position',[279         447        1377         420]);  set(tmpf,'PaperPositionMode','auto');
        ccc = [0:1/length(SynchLvl_LST):1];
        
        for s_ii = 1 : length(SynchLvl_LST)
            p2_ii = s_ii;
            %WT
            spkBin = ACT_Record{p1_ii,p2_ii,p3_ii}.WT.VL.All.spktrain;
            [box,sumAC, cntSample] = AutocorrFromSpktrain(spkBin, Trange, cuttime );
            subplot(121); hold on;
            plot(box,sumAC,'Color',[0 0 ccc(s_ii)]); LEG{s_ii} = [saveTxt2 ' ' num2str(PARAM2(p2_ii)) ];
            %KO
            spkBin = ACT_Record{p1_ii,p2_ii,p3_ii}.KO.VL.All.spktrain;
            [box,sumAC, cntSample] = AutocorrFromSpktrain(spkBin, Trange, cuttime );
            subplot(122); hold on;
            plot(box,sumAC,'Color',[ccc(s_ii) 0 0]); LEG{s_ii} = [saveTxt2 ' ' num2str(PARAM2(p2_ii)) ];
        end
        subplot(121); legend(LEG); title('WT'); subplot(122); legend(LEG); title('KO')
        suptitle([titleTxt1 ' = ' num2str(InputFR_LST(f_ii)) ' , '  'T_r_a_n_g_e = ' num2str(Trange)]);
         if (SAVE_FIG)
            ffig = [ dirLoc dirFig 'CombinedSynchLevel_Trange' num2str(Trange) '_' saveTxt1 num2str(PARAM1(p1_ii))];
            saveas(  tmpf, [ffig '.jpg'], 'jpg')
            saveas(  tmpf, [ffig '.fig'], 'fig')
            
        end
    end
end

