% Quantitatively Test Synchronization Level
clc
close all
clear all


simCode  = 'ResFunc';
CoreName = '';

SAVE_FIG = 1 ; Close_Fig_aftr_save = 1;

Nsample = 100;
SynchLvl = 0;
Tcut = 500;
Tmeasure = 5000;
TSTOP = Tcut + Tmeasure;
tstop = TSTOP;
InputFR = 40;
W_SPK = 0.029;
W_VL_M1 = 0.002/Nsample ;


InputFR_LST = [ 10 25:25:450];
SynchLvl_LST = [0];%[0:1/Nsample:1];
W_VL_M1_LST = [0.0001 : 0.0002 : 0.0001*Nsample]; %[0.0001 : 0.0001 : 0.0001*Nsample];

WspkMult_LST = [1]; %[1 1.5 2 2.5 3  5];
IGmean_LST = 0; %[0 -0.1 -0.5 -1 -1.5]; % 0
IGsig_LST = [0];
OSC_F_LST =  [20, 40, 10];
OSC_AMP_LST = [0 0.5 1];

PARAM1 = InputFR_LST;
lblTxt1 = 'Input Frequency';
saveTxt1 = 'InputFR';
titleTxt1 = 'Input Frequency';
PARAM2 = OSC_F_LST;
lblTxt2 = 'Oscilation Frequencyt';
saveTxt2 = 'OscF';
titleTxt2 = 'OSC Freq';
PARAM3 = OSC_AMP_LST;
lblTxt3 = 'Oscilation Amplitude relative to mean input FR';
saveTxt3 = 'OscAmp';
titleTxt3 = 'OSC Amp';
PARAM4 = WspkMult_LST;
lblTxt4 = 'Weighting factor';
saveTxt4 = 'WspkMult';
titleTxt4 = 'W_m_u_l_t';

N_Param = 4;
ACT_Rec_size = zeros(1,N_Param);
for mm = 1 : N_Param
    PARAMETERS{mm}.PARAM = eval(sprintf('PARAM%d',mm)) ;
    PARAMETERS{mm}.lblTxt = eval(sprintf('lblTxt%d',mm));
    PARAMETERS{mm}.titleTxt = eval(sprintf('titleTxt%d',mm));
    PARAMETERS{mm}.saveTxt = eval(sprintf('saveTxt%d',mm));
    ACT_Rec_size(mm) = eval(sprintf('length(PARAM%d)',mm));
end
if(length(ACT_Rec_size) == 1)
    ACT_Rec_size = [ACT_Rec_size 1];
end
ACT_Record = cell(ACT_Rec_size);
Check_Status = zeros(ACT_Rec_size);
coreName = '';

% Directory
PATH = SetPath;
dirLoc = [PATH 'ResFunc_cellType/'];
dirFig = ['TestSmallResFunc/']; %['smallResFunc_' get_Parameters_RangeTxt( PARAMETERS,[1:4]) '/'];
mkdir([dirLoc dirFig])

for p1_ii = 1 : length(PARAM1)
    for p2_ii = 1 : length(PARAM2)
        for p3_ii = 1 : length(PARAM3)
            for p4_ii = 1 : length(PARAM4)
                plist = [p1_ii, p2_ii, p3_ii,p4_ii];
                
                inF_ii = p1_ii; sl_ii = 1; wVM_ii = 1;
                igm_ii = 1; igs_ii = 1; ww_ii = p4_ii;
                of_ii = p2_ii; oa_ii = p3_ii;
                
                W_VL_M1= W_VL_M1_LST(wVM_ii);
                SynchLvl = SynchLvl_LST(sl_ii);
                InputFR = InputFR_LST(inF_ii);
                CUTTIME = 500;
                TSTOP = 5500;
                
                IGmean = IGmean_LST(igm_ii);
                IGsig = IGsig_LST(igs_ii);
                TestW = WspkMult_LST(ww_ii);
                Wscale = 0.001;
                if(IGmean == 0)
                    Wspk = Wscale*TestW;
                else
                    Wspk =  IGmean*-10*Wscale*TestW;
                end
                
                Osc_F = OSC_F_LST(of_ii);
                Osc_amp = OSC_AMP_LST(oa_ii);
                PHASE =0;
                ampTxt = sprintf('%2.2f',Osc_amp);
                %RecordSpkPvalues_VL_Nsample100_TSTOP5500_InputFR10_Wspk0.003_IGmean0_IGsig0_oscF20Hz_amp0.10_phase0
                simCode =[ coreName 'Nsample' num2str(Nsample) '_TSTOP' num2str(TSTOP) '_InputFR' num2str(InputFR)  ...
                    '_Wspk' num2str(Wspk) '_IGmean' num2str(IGmean) '_IGsig' num2str(IGsig) ...
                    '_oscF' num2str(Osc_F) 'Hz_amp' ampTxt '_phase' num2str(PHASE) ];
                disp(simCode)
                
                %                 '_wSPK' num2str(W_SPK) '_wVLM1_' num2str(W_VL_M1)];
                %             SomaVolt_WT_VL_Nsample100_TSTOP5500_InputFR15
                
                %WT_VL
                Name_postfix  = ['WT_VL_' simCode ];
                VL_BaselineAct
                tmp.E = E;      tmp.I = I;  tmp.All = All; tmp.All.somaVall = somaVall;
                WT.VL = tmp;
                
%                 %KO_VL
%                 Name_postfix  = ['KO_VL_' simCode ];
%                 VL_BaselineAct
%                 tmp.E = E;      tmp.I = I;  tmp.All = All; tmp.All.somaVall = somaVall;
%                 KO.VL = tmp;
%                 
                %KOp_VL
                Name_postfix  = ['KOp_VL_' simCode ];
                VL_BaselineAct
                tmp.E = E;      tmp.I = I;  tmp.All = All; tmp.All.somaVall = somaVall;              
                KOp.VL = tmp;
                
                %Input Spiketrain (same input to WT and KO)
                fInput = [dirLoc 'RecordSpkTrain_VL_' simCode  '.txt'];
                [nE, nI, Tstop,InspkTrain] = get_save_vec_file_BigNet(fInput );
                WT.VL.InputSpikeTrain = InspkTrain;
                KOp.VL.InputSpikeTrain = InspkTrain;
                
                
                ACT.WT =  WT;  %ACT.KO = KO;  
                ACT.KOp = KOp;
                
                ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii} = ACT;
            end
        end
    end
end

%%
PLOT_RASTER = 0;
SAVE_RASTER_FIG =0;
CLOSE_RASTER = 1;

Ncelltype =2;
avgFR = zeros(ACT_Rec_size, Ncelltype ) ;
stdFR = zeros(ACT_Rec_size, Ncelltype ) ;
avgEff = zeros(ACT_Rec_size, Ncelltype) ;
stdEff = zeros(ACT_Rec_size, Ncelltype) ;

avgInputFR = zeros(ACT_Rec_size) ;
stdInputFR = zeros(ACT_Rec_size) ;


for p1_ii = 1 : length(PARAM1)
    for p2_ii = 1 : length(PARAM2)
        for p3_ii = 1 : length(PARAM3)
            for p4_ii = 1 : length(PARAM4)
                
                simTxt = [titleTxt1 ' ' num2str(PARAM1(p1_ii))];
                if(PLOT_RASTER)
                    %raster
                    spkBin = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii}.WT.VL.All.spktrain;
                    [freq_all_WT,spkTime_all_WT, fg_handle_WT ] = raster_from_spkbin( spkBin,cuttime, Tstop, [ simTxt ' : WT']);
                    set(fg_handle_WT, 'position',[  449   450   791   528])
                    spkBin = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii}.KOp.VL.All.spktrain;
                    [freq_all_KOp,spkTime_all_KOp, fg_handle_KOp ] = raster_from_spkbin( spkBin,cuttime, Tstop, [ simTxt ' : KOp']);
                    set(fg_handle_KOp, 'position',[  449   450   791   528])
                    
                    if (SAVE_RASTER_FIG)
                        ffig = [ dirLoc dirFig 'RasterPlot_' get_Parameters_saveText(PARAMETERS, [1:4], [p1_ii p2_ii p3_ii p4_ii]) ];
                        %         ffig = [ dirLoc dirFig 'RasterPlot_Trial' num2str(TRIAL_NO)];
                        saveas( fg_handle_WT, [ffig '_WT.jpg'], 'jpg')
                        saveas( fg_handle_WT, [ffig '_WT.fig'], 'fig')
                        saveas( fg_handle_KOp, [ffig '_KO.jpg'], 'jpg')
                        saveas( fg_handle_KOp, [ffig '_KO.fig'], 'fig')
                    end
                    if(CLOSE_RASTER)
                        close(fg_handle_KOp); close(fg_handle_KO); close(fg_handle_WT);
                    end
                    
                end
                       tmpT = cuttime+1:Tstop;
                       %Input Spike
                spkBin = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii}.KOp.VL.InputSpikeTrain;
                freq_InFR = sum(spkBin(:,tmpT),2)/length(tmpT)*1000;   
                avgInputFR(p1_ii, p2_ii,p3_ii,p4_ii) = mean(freq_InFR);
                stdInputFR(p1_ii, p2_ii,p3_ii,p4_ii) = std(freq_InFR);
                       
                spkBin = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii}.WT.VL.All.spktrain;         
                freq_all_WT = sum(spkBin(:,tmpT),2)/length(tmpT)*1000;
                efficacy_all_WT = freq_all_WT./freq_InFR;
                
                spkBin = ACT_Record{p1_ii, p2_ii, p3_ii,p4_ii}.KOp.VL.All.spktrain;
                freq_all_KOp = sum(spkBin(:,tmpT),2)/length(tmpT)*1000;
                efficacy_all_KOp = freq_all_KOp./freq_InFR;
                
                avgFR(p1_ii, p2_ii,p3_ii,p4_ii, 1) = mean(freq_all_WT); stdFR(p1_ii, p2_ii,p3_ii,p4_ii, 1) = std(freq_all_WT);
                avgFR(p1_ii, p2_ii,p3_ii,p4_ii, 2) = mean(freq_all_KOp);  stdFR(p1_ii, p2_ii,p3_ii,p4_ii, 2) = std(freq_all_KOp);
                
                avgEff(p1_ii, p2_ii,p3_ii,p4_ii, 1) = mean(efficacy_all_WT); stdEff(p1_ii, p2_ii,p3_ii,p4_ii, 1) = std(efficacy_all_WT);
                avgEff(p1_ii, p2_ii,p3_ii,p4_ii, 2) = mean(efficacy_all_KOp);  stdEff(p1_ii, p2_ii,p3_ii,p4_ii, 2) = std(efficacy_all_KOp);
                
            

                disp(['########### Trial = ' simTxt ])                
                disp(['### Input avg freq = ' num2str(mean(freq_InFR))])
                disp(['### WT avg freq = ' num2str(mean(freq_all_WT))])
                disp(['### KOp avg freq = ' num2str(mean(freq_all_KOp))])
                
                
                
            end
        end
    end
end


%% Response Function

% SAVE_FIG =1;
% 
% % Weight and Response function when vary sigma p2_ii vs p3 , p2 vs p4, p3
% % vs p4 
% % p3 vs p4 
% p2_ii = 1; 
% 
% mainPARAM = PARAM4;
% subPARAM = PARAM3;
% for pp = 1 :  length(subPARAM)
%     rf = figure;  set(gcf, 'position',[     372         -94        1160         753])
%     ef = figure;  set(gcf, 'position',[     372         -94        1160         753])  
%     LEG = cell(length(mainPARAM)*2,1);
%     for ii = 1 : length(mainPARAM)        
%         p4_ii = ii; p3_ii = pp;
%         p2_ii = 1;
%         
%         figure(rf);
%         errorbar(PARAM1, avgFR(:,p2_ii,p3_ii,p4_ii,1),stdFR(:,p2_ii,p3_ii,p4_ii,1),[getDotStyle(ii) getLineStyle(1) getColorCode(ii)], 'LineWidth',2); hold on;
%         errorbar(PARAM1, avgFR(:,p2_ii,p3_ii,p4_ii,2),stdFR(:,p2_ii,p3_ii,p4_ii,2),[getDotStyle(ii) getLineStyle(2) getColorCode(ii)], 'LineWidth',2);
%         figure(ef);
%         errorbar(PARAM1, avgEff(:,p2_ii,p3_ii,p4_ii,1),stdEff(:,p2_ii,p3_ii,p4_ii,1),[getDotStyle(ii) getLineStyle(1) getColorCode(ii)], 'LineWidth',2); hold on;
%         errorbar(PARAM1, avgEff(:,p2_ii,p3_ii,p4_ii,2),stdEff(:,p2_ii,p3_ii,p4_ii,2),[getDotStyle(ii) getLineStyle(2) getColorCode(ii)], 'LineWidth',2);
%         
%         if(mainPARAM == PARAM4)
%                 
%         TestW = PARAM4(p4_ii);  IGmean = PARAM2(p2_ii);
%         Wscale = 0.001;
%         if(IGmean == 0)
%             Wspk = Wscale*TestW;
%         else
%             Wspk =  IGmean*-10*Wscale*TestW;
%         end
%         
%         LEG{2*ii-1} =  ['[WT]   W_s_p_k  = ' num2str(Wspk)];
%         LEG{2*ii} = ['[KO2]  W_s_p_k  = ' num2str(Wspk)];
%         end
%     end
%     figure(rf);
%     xlim([0 PARAM1(end)+10]); ylim([0 max(avgFR(:))+10])
%     ylabel('Average Output fring rate (Hz)' ); xlabel('Input Firing rate(Hz)');
%     title({'Response Function of each cell type', get_Parameters_titleText(PARAMETERS, [2,3], [ p2_ii, p3_ii]) })
%     legend(LEG,'location','best');
%     
%     figure(ef);
%     xlim([0 PARAM1(end)+5]); 
%     ylabel('Efficacy (output FR /input FR)' ); xlabel('Input Firing rate(Hz)');
%     title({'Response Function of each cell type', get_Parameters_titleText(PARAMETERS, [2,3], [ p2_ii, p3_ii]) })
%     legend(LEG,'location','best');
%     Line = 1/7 * ones(length(PARAM1),1);
%     plot(PARAM1,Line,':','color',[0.5 0.5 0.5])
%     if(SAVE_FIG)
%         tmpTxt = get_Parameters_saveText(PARAMETERS,[2,3], [ p2_ii, p3_ii]);
%         fg = rf;  figName = ['ResFunc' tmpTxt];
%         saveas(fg,[dirLoc dirFig figName '.fig'],'fig'); saveas(fg, [dirLoc dirFig figName '.jpg'],'jpg');
%         fg = ef;  figName = ['Efficacy' tmpTxt];
%         saveas(fg,[dirLoc dirFig figName '.fig'],'fig'); saveas(fg, [dirLoc dirFig figName '.jpg'],'jpg')
%     end
% end

%% 
if(0)
    
% p4 vs p2 
% close all;
mainPARAM = PARAM2; mainID = 2; 
subPARAM = PARAM4; subID = 4;
dummy = 3;  dd =1;
for pp = 1 :  length(subPARAM)
    rf = figure;  set(gcf, 'position',[     372         -94        1160         753])
    ef = figure;  set(gcf, 'position',[     372         -94        1160         753])  
    LEG = cell(length(mainPARAM)*2,1);
    for mm = 1 : length(mainPARAM)        
%         p4_ii = pp; p2_ii = ii;
%         p3_ii = 1;
        eval(sprintf('p%d_ii = mm; p%d_ii = pp; p%d_ii = dd; ',mainID, subID,dummy ));
        
        figure(rf);
        errorbar(PARAM1, avgFR(:,p2_ii,p3_ii,p4_ii,1),stdFR(:,p2_ii,p3_ii,p4_ii,1),[getDotStyle(mm) getLineStyle(1) getColorCode(mm)], 'LineWidth',2); hold on;
        errorbar(PARAM1, avgFR(:,p2_ii,p3_ii,p4_ii,2),stdFR(:,p2_ii,p3_ii,p4_ii,2),[getDotStyle(mm) getLineStyle(2) getColorCode(mm)], 'LineWidth',2);
        figure(ef);
        errorbar(PARAM1, avgEff(:,p2_ii,p3_ii,p4_ii,1),stdEff(:,p2_ii,p3_ii,p4_ii,1),[getDotStyle(mm) getLineStyle(1) getColorCode(mm)], 'LineWidth',2); hold on;
        errorbar(PARAM1, avgEff(:,p2_ii,p3_ii,p4_ii,2),stdEff(:,p2_ii,p3_ii,p4_ii,2),[getDotStyle(mm) getLineStyle(2) getColorCode(mm)], 'LineWidth',2);
        
        if(mainID == 4)
                
        TestW = PARAM4(p4_ii);  IGmean = IGmean_LST(1);
        Wscale = 0.001;
        if(IGmean == 0)
            Wspk = Wscale*TestW;
        else
            Wspk =  IGmean*-10*Wscale*TestW;
        end
        
        LEG{2*mm-1} =  ['[WT]   W_s_p_k  = ' num2str(Wspk)];
        LEG{2*mm} = ['[KO2]  W_s_p_k  = ' num2str(Wspk)];
        else
             tmpTxt =  get_Parameters_titleText(PARAMETERS,[subID, dummy],[pp, dd]);
             LEG{2*mm-1} =  ['[WT] ' tmpTxt];
             LEG{2*mm} = ['[KO2] ' tmpTxt];
        end
    end
    figure(rf);
    xlim([0 PARAM1(end)+10]); ylim([0 max(avgFR(:))+10])
    ylabel('Average Output fring rate (Hz)' ); xlabel('Input Firing rate(Hz)');
    title({'Response Function of each cell type', get_Parameters_titleText(PARAMETERS, [subID, dummy],[pp, dd]) })
    legend(LEG,'location','best');
    
    figure(ef);
    xlim([0 PARAM1(end)+5]); 
    ylabel('Efficacy (output FR /input FR)' ); xlabel('Input Firing rate(Hz)');
    title({'Response Function of each cell type', get_Parameters_titleText(PARAMETERS, [subID, dummy],[pp, dd]) })
    legend(LEG,'location','best');
    
    if(SAVE_FIG)
        tmpTxt = get_Parameters_saveText(PARAMETERS,[subID, dummy],[pp, dd]);
        fg = rf;  figName = ['ResFunc' tmpTxt];
        saveas(fg,[dirLoc dirFig figName '.fig'],'fig'); saveas(fg, [dirLoc dirFig figName '.jpg'],'jpg');
        fg = ef;  figName = ['Efficacy' tmpTxt];
        saveas(fg,[dirLoc dirFig figName '.fig'],'fig'); saveas(fg, [dirLoc dirFig figName '.jpg'],'jpg')
    end
end
end
%% 
% main 2 sub 4 dummy 2 

testlist = perms(2:4);

% for idid  = 1 : size(testlist,1)
%     tmpid = testlist(idid,:);
tmpid = [2 3 4];
    mainPARAM = eval(sprintf('PARAM%d',tmpid(1)));  mainID = tmpid(1); 
    subPARAM = eval(sprintf('PARAM%d',tmpid(2))); subID = tmpid(2);
    dummy =  tmpid(3);  dd =1;
    disp('============================================================')
    disp(['mainID = ' num2str(mainID) '    subID = ' num2str(subID) '   dummy = ' num2str(dummy)])
    disp('============================================================')
    if (length(mainPARAM)  > 1)
        PLOT_Response_Function_getInSpk;
    end
% k = waitforbuttonpress;
% end


%%
%  tmpRange = 1:7;
%  figure;
% errorbar(PARAM1(tmpRange), avgFR(tmpRange,1),stdFR(tmpRange,1),'.-k','LineWidth',2); hold on;
% errorbar(PARAM1(tmpRange), avgFR(tmpRange,2),stdFR(tmpRange,2),'.-r','LineWidth',2);
% errorbar(PARAM1(tmpRange), avgFR(tmpRange,3),stdFR(tmpRange,3),'.-m','LineWidth',2);
% plot(PARAM1(tmpRange), PARAM1(tmpRange), ':b');
% xlim([0 PARAM1(tmpRange(end))+1]); ylim([0 max(max(avgFR((tmpRange),:)))+10])
% legend('WT','KO1','KO2','Linear relationship');
% ylabel('Average Output fring rate (Hz)' ); xlabel('Input Firing rate(Hz)');
%  title('Response Function of each cell type')
%%
SAVE_NEURON_ACT =0;
coreFileName ='';
if (SAVE_NEURON_ACT)
    save([dirLoc dirFig 'Result_ACT_Record_' coreFileName '.mat'], 'ACT_Record','-v7.3');
end
%%          The VL-M1 correlation

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

