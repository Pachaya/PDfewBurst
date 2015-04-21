
% clc
% close all
% clear all
WT = [];
KO = [];
tmp = [];
%NoiseSTDEV_List = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
RES = 1; %1 bin = 1 ms
VariesINOISE = 0;
VariesInSPK = 0;

SAVE_WORKSPACE = 0;

DoTTest = 1;
SAVE_BASAL_ACT = 0;
FNAME_SIM = 'GPmVL' ; %'PreLim_INOISE0.2_InSPK50_EE1_EI2_IE20_II12_mult125_w0.0006';
plotSampleV = 0;
SPECIFIED_RATIO = 0; %TestSim_WT_InGauss0.2_W0.0005_50.00Hz_rEE0_rEI1_rIE2_rII1_Wmult10
Basal_Act_set  = [];

NoiseMEAN = 0;

SAVE_FIG = 1;
SAVE_RASTER_FIG = 0;

dirLoc = 'Model1_MinParam150324_fullsim/' ; %'TestSim_BasalAct/';  %'VL_LocalConn/';
rEE_rEI_rIE_rII = [1 1 1 1];
w_MULT = 125;
avgFR_CUTTIME = 500;

pulseHz = 1; N_pulse = 1;
Input_I_amp_list = [0];
Input_I_dur_list = [0];
Input_I_amp = -0.3; %[-1; -2; -3; -4; -5;];% [-1; -5]; %
Input_I_dur = 1000; %[10; 20; 30; 40; 50;];
ReboundPeakAmp_WTlist = zeros(length(Input_I_amp_list),length(Input_I_dur_list));
ReboundPeakDel_WTlist = zeros(length(Input_I_amp_list),length(Input_I_dur_list));
ReboundPeakAmp_KOlist = zeros(length(Input_I_amp_list),length(Input_I_dur_list));
ReboundPeakDel_KOlist = zeros(length(Input_I_amp_list),length(Input_I_dur_list));
Save_BasalAct =  cell(length(Input_I_amp_list),length(Input_I_dur_list));
InjStartT = 1500;
InjStopT = InjStartT+Input_I_dur;
reboundDelay = 15;
burstRange = 100;
cutTime = avgFR_CUTTIME;

ADD_I_to_M1 = 0;
ADD_I_to_VL = 0;

FIG_ALL = 0;
dirFig = 'Fig_CheckWeightOnly_delT15_BurstRange100/';
mkdir([dirLoc dirFig ])

TSTOP = 4000;

rTC_LST = [120];
wmTC_LST = [10;];

LightAmp_LST = [0.3; 0.4; 0.5];
%GPmVLw_mean_LST = [0.02:0.02:0.1];
% GPmVLw_sig_LST =[0.005; 0.01; 0.02; 0.03;];
% LightAmp_LST = [0.3; ];
GPmVLw_mean_LST = [0.04:0.02:0.08];
GPmVLw_sig_LST =[ 0.01; 0.02; ];

% Candidates Model 1 :
% I amp = 0.3, mean 0.06, sigma = 0.01, 0.005 – 0.03
% I amp = 0.5, mean 0.04, sigma = 0.03, 0.01 – 0.03
%PDfewBurst_GPmVLmd1_rTC120_wmTC10_WT_GPmInput_Amp0.3_Dur1000_GPmVLw_m0.06_sig0.01_InGauss0.2_IGmean-0.15_IGmeanSig0_W0.0015_SpecifiedPoisSpk_sig0.00Hz_T4000_trial3
NUM_TRIAL = 5;
ncells = 1150;
TRIAL_LST = 1 : NUM_TRIAL;

PARAM1 = LightAmp_LST; % TRIAL_LST; 
lblTxt1 = 'Light Stimulus amplitude (nA)'; %'Trial'; %
saveTxt1 = 'GPmAmp'; %'trial'; 
titleTxt1 = 'Light Amp';
PARAM2 = GPmVLw_mean_LST;
lblTxt2 = 'GPm-VL weight mean';
saveTxt2 = 'GPmVLw_m';
titleTxt2 = 'W_m_e_a_n';
PARAM3 = GPmVLw_sig_LST;
lblTxt3 = 'GPm-VL weight sigma';
saveTxt3 = 'GPmVLw_sig';
titleTxt3 = 'W_s_i_g';

ACT_Record = cell(length(PARAM1),length(PARAM2),length(PARAM3));

N_Param = 3;
for ii = 1 : N_Param
    PARAMETERS{ii}.PARAM = eval(sprintf('PARAM%d',ii)) ;
    PARAMETERS{ii}.lblTxt = eval(sprintf('lblTxt%d',ii));
    PARAMETERS{ii}.titleTxt = eval(sprintf('titleTxt%d',ii));
    PARAMETERS{ii}.saveTxt = eval(sprintf('saveTxt%d',ii));
end

CUTTIME = 500;
PhotoInjT = 1500;
PhotoStop = 2500;
DelayT = 15;
BurstRange = 100;



for p1_ii = 1 : length(PARAM1)
    for p2_ii = 1 : length(PARAM2)
        for p3_ii = 1 : length(PARAM3)
             la_ii = p1_ii; 
            m_ii = p2_ii; s_ii = p3_ii;
            TRIAL_NO = 1;
            
            for cell_type = 1 : 2
                
                if (cell_type == 1)
                    cTxt = 'WT';
                elseif (cell_type == 2)
                    cTxt = 'KO';
                end
                
                % PDfewBurst_GPmVLmd1_rTC120_wmTC10_WT_GPmInput_Amp0.3_Dur1000_GPmVLw_m0.06_sig0.01_InGauss0.2_IGmean-0.15_IGmeanSig0_W0.0015_SpecifiedPoisSpk_sig0.00Hz_T4000_trial3
                
                coreFileName = 'PDfewBurst_GPmVLmd1' ;
                
                InGauss_STDEV = 0.2; %0.2;, 0.3
                NoiseMEAN = -0.15;
                IGmeanSig = 0;
                W_Weight = 0.0012;
                PoisInputFr = 0;
                TSTOP = 4000;
                GPmLightDur = 1000;
                
%                 rTC = 120;
%                 wmTC = 10;
%                 GPmLight = 0.3;
              
                
                                GPmLight = LightAmp_LST(la_ii);
                                GPm_w_mn = GPmVLw_mean_LST(m_ii);
                                GPm_w_sg = GPmVLw_sig_LST(s_ii);
                
                
                Name_postfix = [coreFileName '_rTC' num2str(rTC) '_wmTC' num2str(wmTC) '_' cTxt '_' 'GPmInput_Amp' num2str(GPmLight) '_Dur' num2str(GPmLightDur) '_GPmVLw_m' num2str(GPm_w_mn) '_sig' num2str(GPm_w_sg) ...
                    '_InGauss' num2str(InGauss_STDEV) '_IGmean' num2str(NoiseMEAN) '_IGmeanSig' num2str(IGmeanSig) '_W' num2str(W_Weight) '_SpecifiedPoisSpk_sig0.00Hz_T' num2str(TSTOP) '_trial' num2str(TRIAL_NO)];
                
                %PDfewBurst_GPmVLmd1_rTC120_wmTC10_WT_GPmInput_Amp0.3_Dur1000_GPmVLw_m0.06_sig0.01_InGauss0.2_IGmean-0.15_IGmeanSig0_W0.0015_SpecifiedPoisSpk_sig0.00Hz_T4000_trial3
                
                
                
                disp('==================================================================================================')
                disp(Name_postfix)
                if (cell_type == 1)
                    cTxt = 'WT';
                else
                    cTxt = 'KO';
                end
                figNameCode = [ saveTxt1 num2str(PARAM1(p1_ii))  '_' saveTxt2 num2str(PARAM2(p2_ii)) '_' cTxt];
                
                %                 CheckFileExist( dirLoc, Name_postfix  )
                
                VL_PhotoactivationAll
                
                %                 if(0)
                
                tmp.E = E;      tmp.I = I;  tmp.All = All; %tmp.VL = VL; tmp.M1 = M1;
                if (cell_type == 1)
                    WT = tmp;
                    cTxt = 'WT';
                    Name_postfix_WT = Name_postfix;
                elseif (cell_type == 2)
                    KO = tmp;
                    cTxt = 'KO';
                    Name_postfix_KO = Name_postfix;
                end
                
                
            end
            if (plotSampleV)
                    %soma's volt
                    [NEcell,NIcell, Tstop,somaVall] = get_save_vec_file_BigNet([dirLoc 'SomaVolt_' Name_postfix_WT '.txt'] );
                    samE = 385;
                    samI = 129;
                    samTstr = PhotoInjT-200; samTstp = PhotoStop + DelayT + BurstRange +200; %500+ Input_I_dur +200;
                    fgSmpl = figure; plot(samTstr:samTstp, somaVall(samE,samTstr:samTstp), 'r'); hold on;
                    plot(samTstr:samTstp, somaVall(samE+1,samTstr:samTstp), 'm');
                    %title('WT')
                    
                    if(ADD_I_to_VL)
                        plot(samTstr:samTstp,somaVall(N_E+samI,samTstr:samTstp),'b');
                        plot(samTstr:samTstp,somaVall(N_E+samI+1,samTstr:samTstp),'k');
                    end
                    title([cTxt ': sample membrane potential, Inject I : amp = ' num2str(Input_I_amp) ', dur = ' num2str(Input_I_dur) ])
                    if(SAVE_FIG)
                        saveas(fgSmpl, [dirLoc dirFig 'SamV' figNameCode '.fig'], 'fig')
                        saveas(fgSmpl, [dirLoc dirFig 'SamV' figNameCode '.jpg'], 'jpg')
                    end
            end
                
            Basal_Act.INOISE.MEAN = NoiseMEAN;
            Basal_Act.INOISE.STDEV = InGauss_STDEV;
            
            Basal_Act.WT = WT;
            Basal_Act.KO = KO;
            DoTTest = 0;
            if (DoTTest)
                % two-tailed t-test % during normal
                
                % disp('PoisSPk')
                % [hP,pP] = ttest(WT.PoisSpk.fr_data, KO.PoisSpk.fr_data)
                disp('--------------------------------------------------------------------------------------------------')
                disp('E cells : two-tailed t-test')
                [hE,pE] = ttest2(WT.VL.normal.fr_data, KO.VL.normal.fr_data, 'Tail','both');
                disp([ ' p-value = ' num2str(pE)])
                if(hE)
                    disp('H0 rejected : average normal firing rate of E cells in WT significantly different from KO (p < 0.05)')
                else
                    disp('Do not rejected H0 : average normal firing rate of E cells in WT does not significantly different from KO (p > 0.05)')
                end
                Basal_Act.ttest2.hE = hE;
                Basal_Act.ttest2.pE = pE;
                if(ADD_I_to_VL)
                    disp('--------------------------------------------------------------------------------------------------')
                    disp('I cells : two-tailed t-test')
                    [hI,pI] = ttest2(WT.I.fr_data, KO.I.fr_data);
                    disp([ ' p-value = ' num2str(pI)])
                    if(hI)
                        disp('H0 rejected : average normal firing rate of I cells in WT significantly different from KO (p < 0.05)')
                    else
                        disp('Do not rejected H0 : average normal firing rate of I cells in WT does not significantly different from KO (p > 0.05)')
                    end
                    disp('--------------------------------------------------------------------------------------------------')
                    disp('All cells : two-tailed t-test')
                    [hA,pA]  = ttest2([WT.VL.normal.fr_data; WT.I.fr_data], [KO.VL.normal.fr_data; KO.I.fr_data]);
                    disp([ ' p-value = ' num2str(pA)])
                    if(hA)
                        disp('H0 rejected : average normal firing rate of all cells in WT significantly different from KO (p < 0.05)')
                    else
                        disp('Do not rejected  H0 : average normal firing rate of all cells in WT does not significantly different from KO (p > 0.05)')
                    end
                    Basal_Act.ttest2.hI = hI;
                    Basal_Act.ttest2.pI = pI;
                    Basal_Act.ttest2.hA = hA;
                    Basal_Act.ttest2.pA = pA;
                end
                
            end
            %% Sample T-Test
            DoSampleTtest = 0;
            if(DoSampleTtest)
                
                Nsample = 15;
                Nrepeat = 100;
                disp(['===== Sample Test : #sample = ' num2str(Nsample) ' for ' num2str(Nrepeat) ' trials ====='])
                tresult = cell(Nrepeat,1);
                plist = zeros(Nrepeat,1);
                hlist = zeros(Nrepeat,1);
                tresult_t = cell(Nrepeat,1);
                plist_t = zeros(Nrepeat,1);
                hlist_t = zeros(Nrepeat,1);
                
                for ii = 1: Nrepeat
                    sample = randperm(ncells, Nsample);
                    [ht,pt] = ttest2(WT.VL.normal.fr_data(sample), KO.VL.normal.fr_data(sample));
                    [h,p] = ranksum(WT.VL.normal.fr_data(sample), KO.VL.normal.fr_data(sample));
                    
                    tresult{ii}.sampleID = sample;
                    tresult{ii}.h = h;
                    tresult{ii}.p = p;
                    %     disp(['h:' num2str(h) ', p:' num2str(p)])
                    plist(ii) = p;
                    hlist(ii) = h;
                    
                    tresult_t{ii}.sampleID = sample;
                    tresult_t{ii}.h = ht;
                    tresult_t{ii}.p = pt;
                    %     disp(['ht:' num2str(ht) ', p:' num2str(pt)])
                    plist_t(ii) = pt;
                    hlist_t(ii) = ht;
                end
                
                disp('##Ranksum##')
                %  sum((plist>0.05))
                disp([num2str(sum((plist>0.05))) ' cases do not rejected  H0 : average normal firing rate of all cells in WT does not significantly different from KO (p > 0.05)'])
                disp('##T-test#')
                %  sum((plist_t>0.05))
                disp([num2str(sum((plist_t>0.05))) ' cases do not rejected  H0 : average normal firing rate of all cells in WT does not significantly different from KO (p > 0.05)'])
                % Common cases
                disp('The other cases are rejecting H0 : average normal firing rate of all cells in WT significantly different from KO (p < 0.05)')
                disp('')
                disp('common cases in Ranksum and T-test')
                find((plist>0.05) == (plist_t>0.05))
                Basal_Act.sampleTest.ranksum = tresult;
                Basal_Act.sampleTest.ttest = tresult_t;
                
            end
            
            ACT_Record{p1_ii,p2_ii,p3_ii } = Basal_Act;
            
        end
    end
end




%%


CntBurstSpkWT = zeros(length(PARAM1),length(PARAM2), length(PARAM3));
CntBurstSpkKO = zeros(length(PARAM1),length(PARAM2), length(PARAM3));
ChanceOfBurstingWT = zeros(length(PARAM1),length(PARAM2), length(PARAM3));
ChanceOfBurstingKO = zeros(length(PARAM1),length(PARAM2), length(PARAM3));
for p1_ii = 1 : length(PARAM1)
    fLA1 = figure; set(fLA1,'position',[302          49        1290         948]) ;  set(gcf,'PaperPositionMode','auto')
    fLA2 = figure; set(fLA2,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
    fB1 = figure; set(fB1,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
    fB2 = figure; set(fB2,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
    cnt = 0;
    for p2_ii = 1 : length(PARAM2)
        for p3_ii = 1 : length(PARAM3)
            la_ii = p1_ii; m_ii = p2_ii; s_ii = p3_ii;
            GPmLight = LightAmp_LST(la_ii);
            GPm_w_mn = GPmVLw_mean_LST(m_ii);
            GPm_w_sg = GPmVLw_sig_LST(s_ii);
%               TRIAL_NO = p1_ii;
%                 GPmLight = 0.3;
%                 GPm_w_mn = 0.06; GPm_w_sg = 0.01;
                
            simTxt = [titleTxt1 ' = ' num2str(PARAM1(p1_ii)) ', ' titleTxt2 ' = ' num2str(PARAM2(p2_ii)) ', ' titleTxt3 ' = ' num2str(PARAM3(p3_ii))];  
            
            spkBin = ACT_Record{p1_ii,p2_ii,p3_ii}.WT.All.spktrain;
            baseline_fr_WT = get_avg_baseline(spkBin, cutTime, InjStartT);
            [freq_all_WT,spkTime_all_WT, fg_handle_WT ] = raster_from_spkbin_BurstRange( spkBin,PhotoStop, Tstop,PhotoInjT, PhotoStop, DelayT, BurstRange, [ simTxt ' : WT']);
            set(fg_handle_WT, 'position',[  449   450   791   528])
                      
            spkBin = ACT_Record{p1_ii,p2_ii,p3_ii}.KO.All.spktrain;
            baseline_fr_KO = get_avg_baseline(spkBin, cutTime, InjStartT);
            [freq_all_KO,spkTime_all_KO, fg_handle_KO ] = raster_from_spkbin_BurstRange( spkBin,PhotoStop, Tstop,PhotoInjT, PhotoStop, DelayT, BurstRange, [ simTxt ' : KO']);
            set(fg_handle_KO, 'position',[  449   450   791   528])
                      
            if (SAVE_RASTER_FIG)
                ffig = [ dirLoc dirFig 'RasterPlot_' saveTxt1 num2str(PARAM1(p1_ii)) '_' saveTxt2 num2str(PARAM2(p2_ii)) '_' saveTxt3 num2str(PARAM3(p3_ii))];
                saveas( fg_handle_WT, [ffig '_WT.jpg'], 'jpg')
                saveas( fg_handle_WT, [ffig '_WT.fig'], 'fig')
                saveas( fg_handle_KO, [ffig '_KO.jpg'], 'jpg')
                saveas( fg_handle_KO, [ffig '_KO.fig'], 'fig')
            end
            
            WT_expBaselineSpk = baseline_fr_WT/1000*burstRange;
            WT_burstSpk = ACT_Record{p1_ii,p2_ii,p3_ii}.WT.All.BurstSpk; % BurstSpk  = Number of spike during burstRange
            CntBurstSpkWT(p1_ii,p2_ii,p3_ii) = mean(WT_burstSpk -WT_expBaselineSpk); % Number of spike during burstRange - expected number of spike from baseline activity
            
            KO_expBaselineSpk = baseline_fr_WT/1000*burstRange;
            KO_burstSpk = ACT_Record{p1_ii,p2_ii,p3_ii}.KO.All.BurstSpk; % BurstSpk  = Number of spike during burstRange
            CntBurstSpkKO(p1_ii,p2_ii,p3_ii) = mean(KO_burstSpk - KO_expBaselineSpk); % Number of spike during burstRange - expected number of spike from baseline activity
            
            ChanceOfBurstingWT(p1_ii,p2_ii,p3_ii) =  sum( round(WT_burstSpk -WT_expBaselineSpk) > 0)/(nE+nI);
            ChanceOfBurstingKO(p1_ii,p2_ii,p3_ii) =  sum( round(KO_burstSpk -KO_expBaselineSpk) > 0)/(nE+nI); 
            
            cnt = cnt+1;
            figure(fLA1);
            subplot(length(PARAM2), length(PARAM3),cnt)
            [nb,xb]=hist(ACT_Record{p1_ii,p2_ii,p3_ii}.WT.All.fr_AfterLightOff,10); bh=bar(xb,nb); hold on;
            set(bh,'facecolor',[0 0 0]);
            xlim([0 max(xb)+1])
            title([ 'mean = ' num2str(GPm_w_mn) ', sig = ' num2str(GPm_w_sg)])
            figure(fLA2);
            subplot(length(PARAM2), length(PARAM3),cnt)
            [nb,xb]=hist(ACT_Record{p1_ii,p2_ii,p3_ii}.KO.All.fr_AfterLightOff,10); bh=bar(xb,nb); hold on;
            set(bh,'facecolor',[1 0 0]);
            xlim([0 max(xb)+1])
            title([ 'mean = ' num2str(GPm_w_mn) ', sig = ' num2str(GPm_w_sg)])
            
            figure(fB1);
            subplot(length(PARAM2), length(PARAM3),cnt)
            [nb,xb]=hist(ACT_Record{p1_ii,p2_ii,p3_ii}.WT.All.BurstSpk); bh=bar(xb,nb); hold on;
            set(bh,'facecolor',[0 0 0]);
            xlim([0 max(xb)+1])
            title([ 'mean = ' num2str(GPm_w_mn) ', sig = ' num2str(GPm_w_sg)])
            figure(fB2);
            subplot(length(PARAM2), length(PARAM3),cnt)
            [nb,xb]=hist(ACT_Record{p1_ii,p2_ii,p3_ii}.KO.All.BurstSpk,10); bh=bar(xb,nb); hold on;
            set(bh,'facecolor',[1 0 0]);
            xlim([0 max(xb)+1])
            title([ 'mean = ' num2str(GPm_w_mn) ', sig = ' num2str(GPm_w_sg)])
        end
    end
    figure(fLA1); suptitle(['Avg Fr after light off WT: ' titleTxt1 '=' num2str(PARAM1(p1_ii)) ]);
    figure(fB1); suptitle(['Number of Burst spike WT: ' titleTxt1 '=' num2str(PARAM1(p1_ii)) ]);
    figure(fLA2); suptitle(['Avg Fr after light off KO: ' titleTxt1 '=' num2str(PARAM1(p1_ii)) ]);
    figure(fB2);  suptitle(['Number of Burst spike KO: ' titleTxt1 '=' num2str(PARAM1(p1_ii)) ]);
    if (SAVE_FIG)
        ffig = [ dirLoc dirFig 'AvgFrLightOff_' saveTxt1 num2str(PARAM1(p1_ii))];
        saveas( fLA1, [ffig '_WT.jpg'], 'jpg')
        saveas( fLA1, [ffig '_WT.fig'], 'fig')
        saveas( fLA2, [ffig '_KO.jpg'], 'jpg')
        saveas( fLA2, [ffig '_KO.fig'], 'fig')
        ffig = [ dirLoc dirFig 'NumBurstSpk_' saveTxt1 num2str(PARAM1(p1_ii))];
        saveas( fB1, [ffig '_WT.jpg'], 'jpg')
        saveas( fB1, [ffig '_WT.fig'], 'fig')
        saveas( fB2, [ffig '_KO.jpg'], 'jpg')
        saveas( fB2, [ffig '_KO.fig'], 'fig')
    end
end

%
for p1_ii = 1 : length(PARAM1)
    Bfg = figure;  set(Bfg,'position',[680   283   712   695]); set(gcf,'PaperPositionMode','auto');
    Bfg2 = figure;  set(Bfg2,'position',[680   283   712   695]);  set(gcf,'PaperPositionMode','auto');
    
    for p2_ii = 1 : length(PARAM2)
        figure(Bfg)
        subplot(3,2, p2_ii); hold on;
        plot( PARAM3, squeeze(CntBurstSpkWT(p1_ii,p2_ii,:)),'*-k');
        plot( PARAM3, squeeze(CntBurstSpkKO(p1_ii,p2_ii,:)),'*-r');
        legend('WT','KO','location','Best')
        title([titleTxt2 ' = ' num2str(PARAM2(p2_ii))]);
        xlabel(titleTxt3); ylabel('Avg Burst spike')
        
        figure(Bfg2)
        subplot(3,2, p2_ii); hold on;
        plot( PARAM3, squeeze(ChanceOfBurstingWT(p1_ii,p2_ii,:)),'*-k');
        plot( PARAM3, squeeze(ChanceOfBurstingKO(p1_ii,p2_ii,:)),'*-r');
        legend('WT','KO','location','Best')
        title([titleTxt3 ' = ' num2str(PARAM2(p2_ii))]);
        xlabel(titleTxt3); ylabel('Chance of Burst spike')
        
    end
    figure(Bfg)
    suptitle([ lblTxt1 ' = ' num2str(PARAM1(p1_ii))]);
    figure(Bfg2)
    suptitle([ 'Chance of Bursting, ' lblTxt1 ' = ' num2str(PARAM1(p1_ii))]);
    
    if (SAVE_FIG)
        ffig = [ dirLoc dirFig 'NumBurstSpk_' saveTxt1 num2str(PARAM1(p1_ii)) ];
        saveas(  Bfg, [ffig '.jpg'], 'jpg');        saveas(  Bfg, [ffig '.fig'], 'fig');
        ffig = [ dirLoc dirFig 'ChanceOfBurst_' saveTxt1 num2str(PARAM1(p1_ii)) ];
        saveas(  Bfg2, [ffig '.jpg'], 'jpg');        saveas(  Bfg2, [ffig '.fig'], 'fig');
    end
end

