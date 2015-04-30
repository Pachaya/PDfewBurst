clc
close all
clear all
WT = [];
KO = [];
tmp = [];
%NoiseSTDEV_List = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
RES = 1; %1 bin = 1 ms
VariesINOISE = 0;
VariesInSPK = 0;

SAVE_WORKSPACE = 0; SPECIFIED_BASAL_ACT_CODENAME = '';

DoTTest = 1; DoSampleTtest = 1;
SAVE_BASAL_ACT = 1;
plotSampleV = 0;
SPECIFIED_RATIO = 0; %TestSim_WT_InGauss0.2_W0.0005_50.00Hz_rEE0_rEI1_rIE2_rII1_Wmult10


NoiseMEAN = 0;

SAVE_FIG = 1; Close_Fig_aftr_save = 1;

dirLoc = 'Min150403_SumWTC/';
% dirLoc = 'Model1_MinParam150403_testTC/' ;
rEE_rEI_rIE_rII = [1 1 1 1];
w_MULT = 125;
avgFR_CUTTIME = 500;

pulseHz = 1; N_pulse = 1;
Input_I_amp_list = [0];
Input_I_dur_list = [0];
Input_I_amp = -0.3; %[-1; -2; -3; -4; -5;];% [-1; -5]; %
Input_I_dur = 500; %[10; 20; 30; 40; 50;];
ReboundPeakAmp_WTlist = zeros(length(Input_I_amp_list),length(Input_I_dur_list));
ReboundPeakDel_WTlist = zeros(length(Input_I_amp_list),length(Input_I_dur_list));
ReboundPeakAmp_KOlist = zeros(length(Input_I_amp_list),length(Input_I_dur_list));
ReboundPeakDel_KOlist = zeros(length(Input_I_amp_list),length(Input_I_dur_list));
Save_BasalAct =  cell(length(Input_I_amp_list),length(Input_I_dur_list));
InjStartT = 1500;
InjStopT = InjStartT+Input_I_dur;
reboundDelay = 0;
burstRange = 100;
cutTime = avgFR_CUTTIME;

ADD_I_to_M1 = 0;
ADD_I_to_VL = 0;

FIG_ALL = 0;
dirFig = 'TestRunDelete/' ; %'specified_wmTC0.002div_forVL/';
mkdir([dirLoc dirFig])
TSTOP = 4000;

rTC_LST = [10 50:50:300];
wmTC_LST = [1];

sumW = 0.002;
specified_wmTC_LST = sumW./[1 0.6 0.2];

LightAmp_LST = [0.5; ];
GPmVLw_mean_LST = [0.5;];
GPmVLw_sig_LST =[ 0 ];

% Candidates Model 1 :
% I amp = 0.3, mean 0.06, sigma = 0.01, 0.005 ?0.03
% I amp = 0.5, mean 0.04, sigma = 0.03, 0.01 ?0.03
%PDfewBurst_GPmVLmd1_rTC120_wmTC10_WT_GPmInput_Amp0.3_Dur1000_GPmVLw_m0.06_sig0.01_InGauss0.2_IGmean-0.15_IGmeanSig0_W0.0015_SpecifiedPoisSpk_sig0.00Hz_T4000_trial3
NUM_TRIAL = 5;
ncells = 1150;
TRIAL_LST = 1 : NUM_TRIAL;

PARAM1 = rTC_LST;
lblTxt1 = 'Range of thalamocortical connection';
saveTxt1 = 'rTC';
titleTxt1 = 'Range_T_C';
PARAM2 = specified_wmTC_LST;
lblTxt2 = 'Specified summation of weight of thalamorcortical connection'; %'Weight of thalamocortical connection';
saveTxt2 = 'specified_wmTC';
titleTxt2 = 'sum W_T_C';
PARAM3 = LightAmp_LST; % TRIAL_LST;
lblTxt3 = 'Light Stimulus amplitude (nA)'; %'Trial'; %
saveTxt3 = 'GPmAmp'; %'trial';
titleTxt3 = 'Light Amp';
PARAM4 = GPmVLw_mean_LST;
lblTxt4 = 'GPm-VL weight mean';
saveTxt4 = 'GPmVLw_m';
titleTxt4 = 'W_m_e_a_n';
PARAM5 = GPmVLw_sig_LST;
lblTxt5 = 'GPm-VL weight sigma';
saveTxt5 = 'GPmVLw_sig';
titleTxt5 = 'W_s_i_g';

N_Param = 5;
ACT_Rec_size = zeros(1,N_Param);
for ii = 1 : N_Param
    PARAMETERS{ii}.PARAM = eval(sprintf('PARAM%d',ii)) ;
    PARAMETERS{ii}.lblTxt = eval(sprintf('lblTxt%d',ii));
    PARAMETERS{ii}.titleTxt = eval(sprintf('titleTxt%d',ii));
    PARAMETERS{ii}.saveTxt = eval(sprintf('saveTxt%d',ii));
    ACT_Rec_size(ii) = eval(sprintf('length(PARAM%d)',ii));
end
ACT_Record = cell(ACT_Rec_size);
% ACT_Rec_size = [ACT_Rec_size 2];
% Check_Status = zeros(ACT_Rec_size);

CUTTIME = 500;
PhotoInjT = 1500;
PhotoStop = 2500;
DelayT = 0;
BurstRange = 100;


%CheckSimulationByName_Combine_GpmVL_VLM1

for p1_ii = 1 : length(PARAM1)
    for p2_ii = 1 : length(PARAM2)
        for p3_ii = 1 : length(PARAM3)
            for p4_ii = 1 : length(PARAM4)
                for p5_ii = 1 : length(PARAM5)
                    
                    r_ii = p1_ii; swm_ii = p2_ii;
                    la_ii = p3_ii;  m_ii = p4_ii; s_ii = p5_ii;
                    TRIAL_NO = 1;
                    
                    for cell_type = 1 : 2
                        
                        if (cell_type == 1)
                            cTxt = 'WT';
                        elseif (cell_type == 2)
                            cTxt = 'KO';
                        end
                        
                        
                     coreFileName = 'PDfewBurst_GPmVLmd1_0del' ;
                        
                        InGauss_STDEV = 0.2; %0.2;, 0.3
                        NoiseMEAN = 0;
                        IGmeanSig = 0;
                        W_Weight = 0.029;
                        PoisInputFr = 0;
                        TSTOP = 4000;
                        GPmLightDur = 1000;
                        
                        rTC =rTC_LST(r_ii);
%                         wmTC = wmTC_LST(wm_ii);
                        specified_wmTC = specified_wmTC_LST(swm_ii);
                        
                        
                        GPmLight = LightAmp_LST(la_ii);
                        GPm_w_mn = GPmVLw_mean_LST(m_ii);
                        GPm_w_sg = GPmVLw_sig_LST(s_ii);
                        
                        
                        tmpstrW = sprintf('%g',specified_wmTC);                         
                        Simulation_Code = [coreFileName '_rTC' num2str(rTC) '_specified_wmTC' tmpstrW '_' cTxt '_' 'GPmInput_Amp' num2str(GPmLight) '_Dur' num2str(GPmLightDur) '_GPmVLw_m' num2str(GPm_w_mn) '_sig' num2str(GPm_w_sg) ...
                           '_InGauss' num2str(InGauss_STDEV) '_IGmean' num2str(NoiseMEAN) '_IGmeanSig' num2str(IGmeanSig) '_W' num2str(W_Weight) '_SpecifiedPoisSpk_sig0.00Hz_T' num2str(TSTOP) '_trial' num2str(TRIAL_NO)];                        
                        
                       % PDfewBurst_GPmVLmd1_0del_rTC300_specified_wmTC0.001454_KO_GPmInput_Amp0.5_Dur1000_GPmVLw_m0_sig0.5_InGauss0.2_IGmean0_IGmeanSig0_W0.029_SpecifiedPoisSpk_sig0.00Hz_T40_trial1
                        
                
                        disp('==================================================================================================')
                        disp(Simulation_Code)
                        disp('######  Download M1 ')
                        %M1
                        Name_postfix = [ 'M1_' Simulation_Code];
                        getCenterPart =1;
                        PhotoactivationAll% BaselineAct % Get all cell activity
                        
                        tmp.E = E;      tmp.I = I;  tmp.All = All; %tmp.VL = VL; tmp.M1 = M1;
                        if (cell_type == 1)
                            M1.WT = tmp;
                            cTxt = 'WT';
                        elseif (cell_type == 2)
                            M1.KO = tmp;
                            cTxt = 'KO';
                        end
                        
                        
                        figNameCode = [ saveTxt1 num2str(PARAM1(p1_ii))  '_' saveTxt2 num2str(PARAM2(p2_ii)) '_' cTxt];
                        
                        %                 CheckFileExist( dirLoc, Name_postfix  )
                        disp('######  Download VL ')
                        Name_postfix = [ Simulation_Code];
                        getCenterPart =0;
                        PhotoactivationAll %VL_PhotoactivationAll
                        
                        %                 if(0)
                        
                        tmp.E = E;      tmp.I = I;  tmp.All = All; %tmp.VL = VL; tmp.M1 = M1;
                        if (cell_type == 1)
                            VL.WT = tmp;
                            cTxt = 'WT'; Name_postfix_WT = Name_postfix;
                            clear tmp
                        elseif (cell_type == 2)
                            VL.KO = tmp;
                            cTxt = 'KO'; Name_postfix_KO = Name_postfix;
                            clear tmp
                        end
                        
                    end
                    
                    if (plotSampleV)
                        %soma's volt
                        tt= tic();
                        
                        samVL = 385;
                        samM1 = 87;
                        samTstr = PhotoInjT-100; samTstp = PhotoStop + DelayT+BurstRange+100;
                        fgSmpl = figure; set(fgSmpl,'position',[  377          94        1310         894]);  set(gcf,'PaperPositionMode','auto')
                        subplot(221);
                        [NEcell,NIcell, Tstop,somaVall] = get_save_vec_file_BigNet([dirLoc 'SomaVolt_M1_' Name_postfix_WT '.txt'] );
                        plot(samTstr:samTstp, somaVall(samM1,samTstr:samTstp), 'r'); hold on;
                        plot(samTstr:samTstp, somaVall(samM1+1,samTstr:samTstp), 'm');
                        ylim([-90 50]); xlim([samTstr samTstp]);
                        title(['WT M1']);
                        subplot(222);
                        [NEcell,NIcell, Tstop,somaVall] = get_save_vec_file_BigNet([dirLoc 'SomaVolt_M1_' Name_postfix_KO '.txt'] );
                        plot(samTstr:samTstp, somaVall(samM1,samTstr:samTstp), 'r'); hold on;
                        plot(samTstr:samTstp, somaVall(samM1+1,samTstr:samTstp), 'm');
                        ylim([-90 50]); xlim([samTstr samTstp]);
                        title(['KO M1']);
                        subplot(223);
                        [NEcell,NIcell, Tstop,somaVall] = get_save_vec_file_BigNet([dirLoc 'SomaVolt_' Name_postfix_WT '.txt'] );
                        plot(samTstr:samTstp, somaVall(samVL,samTstr:samTstp), 'r'); hold on;
                        plot(samTstr:samTstp, somaVall(samVL+1,samTstr:samTstp), 'm');
                        ylim([-90 50]); xlim([samTstr samTstp]);
                        title(['WT VL']);
                        subplot(224);
                        [NEcell,NIcell, Tstop,somaVall] = get_save_vec_file_BigNet([dirLoc 'SomaVolt_' Name_postfix_KO '.txt'] );
                        plot(samTstr:samTstp, somaVall(samVL,samTstr:samTstp), 'r'); hold on;
                        plot(samTstr:samTstp, somaVall(samVL+1,samTstr:samTstp), 'm');
                        ylim([-90 50]); xlim([samTstr samTstp]);
                        title(['KO VL']);
                        simTxt = get_Parameters_titleText(PARAMETERS, [1:5], [p1_ii p2_ii p3_ii p4_ii p5_ii]);
                        suptitle(simTxt)
                        disp('Finish plot sample membrain voltage');
                        toc(tt);
                        if(SAVE_FIG)
                            simTxt =  get_Parameters_saveText(PARAMETERS, [1:5], [p1_ii p2_ii p3_ii p4_ii p5_ii]);
                            saveas(fgSmpl, [dirLoc dirFig 'SamV' simTxt '.fig'], 'fig')
                            saveas(fgSmpl, [dirLoc dirFig 'SamV' simTxt '.jpg'], 'jpg')
                            if (Close_Fig_aftr_save)
                                close(fgSmpl);
                            end
                        end
                        
                        disp('Finish save sample membrain voltage');
                        toc(tt)
                    end
                    
                    
                    
                    Basal_Act.VL = VL;
                    Basal_Act.M1 = M1;
                    clear VL
                    clear M1
                    
                    if (DoTTest)  %for Baseline activity
                        % two-tailed t-test % during normal
                        
                        % disp('PoisSPk')
                        % [hP,pP] = ttest(WT.PoisSpk.fr_data, KO.PoisSpk.fr_data)
                        disp('--------------------------------------------------------------------------------------------------')
                        disp('Two-tailed t-test')
                        tmpWTfr = Basal_Act.VL.WT.All.fr_data;
                        tmpKOfr = Basal_Act.VL.KO.All.fr_data;
                        
                        [hE,pE] = ttest2(tmpWTfr, tmpKOfr, 'Tail','both');
                        disp([ ' p-value = ' num2str(pE)])
                        if(hE)
                            disp('H0 rejected : average normal firing rate of E cells in WT significantly different from KO (p < 0.05)')
                        else
                            disp('Do not rejected H0 : average normal firing rate of E cells in WT does not significantly different from KO (p > 0.05)')
                        end
                        Basal_Act.ttest2.hE = hE;
                        Basal_Act.ttest2.pE = pE;
                        
                        
                    end
                    %% Sample T-Test for Baseline activity
                    
                    if(DoSampleTtest)
                        tmpWTfr = Basal_Act.VL.WT.All.fr_data;
                        tmpKOfr = Basal_Act.VL.KO.All.fr_data;
                        
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
                            [ht,pt] = ttest2(tmpWTfr(sample), tmpKOfr(sample));
                            [h,p] = ranksum(tmpWTfr(sample), tmpKOfr(sample));
                            
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
                    
                    ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii } = Basal_Act;
                    clear Basal_Act
                end
            end
        end
    end
end

%%

if(SAVE_BASAL_ACT)
    save([dirLoc dirFig 'Saved_Activity_result_' date '.mat' ], 'ACT_Record','PARAMETERS','-v7.3');
   % load([dirLoc dirFig 'Saved_Activity_result_' date '.mat' ]) ; For
   % loading 
end

%%
%% Plot the instantaneous average firing rate

SIGMA_fr_All = 30;


for p1_ii = 1 : length(PARAM1)
    for p2_ii = 1 : length(PARAM2)
        for p3_ii = 1 : length(PARAM3)
            for p4_ii = 1 : length(PARAM4)
                for p5_ii = 1 : length(PARAM5)
                    
                    %%% for VL
                    tmpWTtrain_VL =  ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii}.VL.WT.All.spktrain;
                    tmpKOtrain_VL =  ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii}.VL.KO.All.spktrain;
                    %Compare WT and KO
                    avgFR_WT_VL = getInstantaneousFiringRate(tmpWTtrain_VL, SIGMA_fr_All, RES);
                    avgFR_KO_VL = getInstantaneousFiringRate(tmpKOtrain_VL, SIGMA_fr_All, RES);
                    
                    %%% for M1
                    tmpWTtrain_M1 =  ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii}.M1.WT.All.spktrain;
                    tmpKOtrain_M1 =  ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii}.M1.KO.All.spktrain;
                    %Compare WT and KO
                    avgFR_WT_M1 = getInstantaneousFiringRate(tmpWTtrain_M1, SIGMA_fr_All, RES);
                    avgFR_KO_M1 = getInstantaneousFiringRate(tmpKOtrain_M1, SIGMA_fr_All, RES);
                    
                    fFR = figure; set(gcf, 'position', [ 624   256   733   722]); set(gcf,'PaperPositionMode','auto');
                    subplot(211); plot(avgFR_WT_M1,'k');  hold on;   plot(avgFR_KO_M1,'m');  hold on;
                    title('M1');   ylabel('Average firing rate (Hz)');            xlabel('Time(ms)'); legend('WT','KO','location','best');
                    
                    subplot(212); plot(avgFR_WT_VL,'k');  hold on;   plot(avgFR_KO_VL,'m');  hold on;
                    title('VL');   ylabel('Average firing rate (Hz)');            xlabel('Time(ms)'); legend('WT','KO','location','best');
                    
                    tmpTxt = get_Parameters_titleText(PARAMETERS, [1:5], [p1_ii p2_ii p3_ii p4_ii p5_ii]);
                    suptitle({'Average firing rate Fr(t)', tmpTxt})
                    
                    if (SAVE_FIG)
                        tmpTxt = get_Parameters_saveText(PARAMETERS, [1:5], [p1_ii p2_ii p3_ii p4_ii p5_ii]);
                        ffig = [ dirLoc dirFig 'Inst_Avg_FR_' tmpTxt ];
                        saveas(  fFR, [ffig '.jpg'], 'jpg');        saveas(  fFR, [ffig '.fig'], 'fig');
                        if (Close_Fig_aftr_save)
                            close(fFR)
                        end
                    end
                end
            end
        end
    end
end
%%
%%%% Raster plot of M1
SAVE_RASTER_FIG = 1;
for p1_ii = 1 : length(PARAM1)
    for p2_ii = 1 : length(PARAM2)
        for p3_ii = 1 : length(PARAM3)
            for p4_ii = 1 : length(PARAM4)
                for p5_ii = 1 : length(PARAM5)
                    
                    BasalAct = ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii}.M1;
                    spkBin = BasalAct.WT.All.spktrain;
                    simTxt =  get_Parameters_titleText(PARAMETERS, [1:5], [p1_ii p2_ii p3_ii p4_ii p5_ii]);
                    [freq_all_WT,spkTime_all_WT, fg_handle_WT ] = raster_from_spkbin_BurstRange( spkBin,CUTTIME, Tstop,PhotoInjT, PhotoStop, DelayT, BurstRange, [ simTxt ' : WT']);
                    set(fg_handle_WT, 'position',[  449   450   791   528])
                    
                    spkBin = BasalAct.KO.All.spktrain;
                    baseline_fr_KO = get_avg_baseline(spkBin, cutTime, InjStartT);
                    [freq_all_KO,spkTime_all_KO, fg_handle_KO ] = raster_from_spkbin_BurstRange( spkBin,CUTTIME, Tstop,PhotoInjT, PhotoStop, DelayT, BurstRange, [ simTxt ' : KO']);
                    set(fg_handle_KO, 'position',[  449   450   791   528])
                    
                    if (SAVE_RASTER_FIG)
                        tmpTxt = get_Parameters_saveText(PARAMETERS, [1:5], [p1_ii p2_ii p3_ii p4_ii p5_ii]);
                        ffig = [ dirLoc dirFig 'RasterPlotM1_' tmpTxt];
                        saveas( fg_handle_WT, [ffig '_WT.jpg'], 'jpg');
                        saveas( fg_handle_WT, [ffig '_WT.fig'], 'fig');
                        saveas( fg_handle_KO, [ffig '_KO.jpg'], 'jpg');
                        saveas( fg_handle_KO, [ffig '_KO.fig'], 'fig');
                        if (Close_Fig_aftr_save)
                            close(fg_handle_WT)
                            close(fg_handle_KO)
                            
                        end
                        
                    end
                end
            end
        end
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    GPM - VL Connnection
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rasterplot of VL  and GPm-VL chance of burst
SAVE_RASTER_FIG =0;
% GPm VL
CntBurstSpkWT = zeros(ACT_Rec_size);
CntBurstSpkKO = zeros(ACT_Rec_size);
ChanceOfBurstingWT = zeros(ACT_Rec_size);
ChanceOfBurstingKO = zeros(ACT_Rec_size);
TimeofFirstReboundSpk = zeros(ACT_Rec_size);
p1_ii = 1; p2_ii=1;
for p3_ii = 1 : length(PARAM3)
    fLA1 = figure; set(fLA1,'position',[302          49        1290         948]) ;  set(gcf,'PaperPositionMode','auto')
    fLA2 = figure; set(fLA2,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
    fB1 = figure; set(fB1,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
    fB2 = figure; set(fB2,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
    fRB = figure; set(fRB,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
    fBspk = figure; set(fBspk,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
    cnt = 0;
    nR = length(PARAM4); nC = length(PARAM5);
    for p4_ii = 1 : length(PARAM4)
        for p5_ii = 1 : length(PARAM5)
            la_ii = p3_ii; m_ii = p4_ii; s_ii = p5_ii;
            
            
            GPmLight = LightAmp_LST(la_ii);
            GPm_w_mn = GPmVLw_mean_LST(m_ii);
            GPm_w_sg = GPmVLw_sig_LST(s_ii);
            %               TRIAL_NO = p1_ii;
            %                 GPmLight = 0.3;
            %                 GPm_w_mn = 0.06; GPm_w_sg = 0.01;
            BasalAct = ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii}.VL;
            
            simTxt = ['VL ' titleTxt3 ' = ' num2str(PARAM3(p3_ii)) ', ' titleTxt4 ' = ' num2str(PARAM4(p4_ii)) ', ' titleTxt5 ' = ' num2str(PARAM5(p5_ii))];
            
            spkBin = BasalAct.WT.All.spktrain;
             baseline_fr_WT = get_avg_baseline(spkBin, cutTime, InjStartT);
%             [freq_all_WT,spkTime_all_WT, fg_handle_WT ] = raster_from_spkbin_BurstRange( spkBin,PhotoStop, Tstop,PhotoInjT, PhotoStop, DelayT, BurstRange, [ simTxt ' : WT']);
%             set(fg_handle_WT, 'position',[  449   450   791   528])
            
            spkBin = BasalAct.KO.All.spktrain;
             baseline_fr_KO = get_avg_baseline(spkBin, cutTime, InjStartT);
%             [freq_all_KO,spkTime_all_KO, fg_handle_KO ] = raster_from_spkbin_BurstRange( spkBin,PhotoStop, Tstop,PhotoInjT, PhotoStop, DelayT, BurstRange, [ simTxt ' : KO']);
%             set(fg_handle_KO, 'position',[  449   450   791   528])
            
            if (SAVE_RASTER_FIG)
                tmpTxt = get_Parameters_saveText(PARAMETERS, [3:5], [ p3_ii p4_ii p5_ii]);
                ffig = [ dirLoc dirFig 'RasterPlot_' tmpTxt];
                ttt = tic(); saveas( fg_handle_WT, [ffig '_WT.jpg'], 'jpg'); toc(ttt);
                ttt = tic(); saveas( fg_handle_WT, [ffig '_WT.fig'], 'fig'); toc(ttt);
                saveas( fg_handle_KO, [ffig '_KO.jpg'], 'jpg')
                saveas( fg_handle_KO, [ffig '_KO.fig'], 'fig')
                if (Close_Fig_aftr_save)
                   close(fg_handle_WT);                        close(fg_handle_KO);
                end
                
            end
            
            WT_expBaselineSpk = baseline_fr_WT/1000*burstRange;
            WT_burstSpk = BasalAct.WT.All.BurstSpk; % BurstSpk  = Number of spike during burstRange
            CntBurstSpkWT(p1_ii,p2_ii,p3_ii,p1_ii,p2_ii,p3_ii) = mean(WT_burstSpk - WT_expBaselineSpk); % Number of spike during burstRange - expected number of spike from baseline activity
            
            KO_expBaselineSpk = baseline_fr_KO/1000*burstRange;
            KO_burstSpk = BasalAct.KO.All.BurstSpk; % BurstSpk  = Number of spike during burstRange
            CntBurstSpkKO(p1_ii,p2_ii,p3_ii,p4_ii,p5_ii) = mean(KO_burstSpk - KO_expBaselineSpk); % Number of spike during burstRange - expected number of spike from baseline activity
            
            ChanceOfBurstingWT(p1_ii,p2_ii,p3_ii,p4_ii,p5_ii) =  sum( round(WT_burstSpk -WT_expBaselineSpk) > 0)/(nE+nI);
            ChanceOfBurstingKO(p1_ii,p2_ii,p3_ii,p4_ii,p5_ii) =  sum( round(KO_burstSpk -KO_expBaselineSpk) > 0)/(nE+nI);
            
            cnt = cnt+1;
            figure(fLA1);
            subplot( nR,nC,cnt)
            [nb,xb]=hist(BasalAct.WT.All.fr_AfterLightOff,10); bh=bar(xb,nb); hold on;
            set(bh,'facecolor',[0 0 0]);
            xlim([0 max(xb)+1])
            title([ 'mean = ' num2str(GPm_w_mn) ', sig = ' num2str(GPm_w_sg)])
            figure(fLA2);
            subplot( nR,nC,cnt)
            [nb,xb]=hist(BasalAct.KO.All.fr_AfterLightOff,10); bh=bar(xb,nb); hold on;
            set(bh,'facecolor',[1 0 0]);
            xlim([0 max(xb)+1])
            title([ 'mean = ' num2str(GPm_w_mn) ', sig = ' num2str(GPm_w_sg)])
            
            figure(fB1);
            subplot( nR,nC,cnt)
            [nb,xb]=hist(BasalAct.WT.All.BurstSpk); bh=bar(xb,nb); hold on;
            set(bh,'facecolor',[0 0 0]);
            xlim([0 max(xb)+1])
            title([ 'mean = ' num2str(GPm_w_mn) ', sig = ' num2str(GPm_w_sg)])
            figure(fB2);
            subplot( nR,nC,cnt)
            [nb,xb]=hist(BasalAct.KO.All.BurstSpk,10); bh=bar(xb,nb); hold on;
            set(bh,'facecolor',[1 0 0]);
            xlim([0 max(xb)+1])
            title([ 'mean = ' num2str(GPm_w_mn) ', sig = ' num2str(GPm_w_sg)])
            
            
            % Distribution of first rebound spike timing 
            tmpRebound_WT = zeros(ncells,1); tmpRebound_KO = zeros(ncells,1);
            for id = 1 : ncells
                tmp = find(BasalAct.WT.All.BurstSpkTrain(id,:),1,'first');
                if ~isempty(tmp) 
                    tmpRebound_WT(id) = tmp;                  
                end
                
                tmp = find(BasalAct.KO.All.BurstSpkTrain(id,:),1,'first');
                if ~isempty(tmp) 
                    tmpRebound_KO(id) = tmp;                  
                end                
            end      
             tmpRebound_WT = tmpRebound_WT(tmpRebound_WT ~= 0);
             tmpRebound_KO = tmpRebound_KO(tmpRebound_KO ~= 0);
            figure(fRB); 
            
            subplot(nR,nC,cnt);
            [N,X] = hist(tmpRebound_WT,100);
            bar(X,N/ncells*100,'FaceColor','k','EdgeColor','w'); hold on; 
            [N,X] = hist(tmpRebound_KO,100);
            bar(X,N/ncells*100,'FaceColor','r','EdgeColor','w'); hold on; 
            ylabel('% Neurons'); xlabel('Time(ms)')
            title([ titleTxt4 ' = ' num2str(PARAM4(p4_ii)) ', ' titleTxt5 ' = ' num2str(PARAM5(p5_ii))])
            
            % Distribution of spikes in bursting range 
           
             tmpRebound_WT = sum(BasalAct.WT.All.BurstSpkTrain);
             tmpRebound_KO = sum(BasalAct.KO.All.BurstSpkTrain);
             xx = 1:1:size(BasalAct.WT.All.BurstSpkTrain,2);
            figure(fBspk); 
            subplot(nR,nC,cnt);
            bar(xx,tmpRebound_WT/ncells*100,'k'); hold on; bar(xx,tmpRebound_KO/ncells*100,'FaceColor','r','EdgeColor','r');
            ylabel('% Neurons'); xlabel('Time(ms)')
            title([ titleTxt4 ' = ' num2str(PARAM4(p4_ii)) ', ' titleTxt5 ' = ' num2str(PARAM5(p5_ii))])
            
            clear BasalAct spkBin
        end
    end
    figure(fLA1); suptitle(['Avg Fr after light off WT: ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]);
    figure(fB1); suptitle(['Number of Burst spike WT: ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]);
    figure(fLA2); suptitle(['Avg Fr after light off KO: ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]);
    figure(fB2);  suptitle(['Number of Burst spike KO: ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]);
    figure(fRB); suptitle(['Latency of first Rebound spike : ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]);  % Note equal to the latency of peak firing rate --> latency of peak firing rate calculate from instantaneous firing rate
    figure(fBspk); suptitle(['Percent of neuron that fire during bursting period : ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]);  % Note equal to the latency of peak firing rate --> latency of peak firing rate calculate from instantaneous firing rate
    
    if (SAVE_FIG)
        tmpTxt = get_Parameters_saveText(PARAMETERS, [3], [ p3_ii ]);
        ffig = [ dirLoc dirFig 'AvgFrLightOff' tmpTxt];
        saveas( fLA1, [ffig '_WT.jpg'], 'jpg')
        saveas( fLA1, [ffig '_WT.fig'], 'fig')
        saveas( fLA2, [ffig '_KO.jpg'], 'jpg')
        saveas( fLA2, [ffig '_KO.fig'], 'fig')        
        ffig = [ dirLoc dirFig 'NumBurstSpk'  tmpTxt];
        saveas( fB1, [ffig '_WT.jpg'], 'jpg')
        saveas( fB1, [ffig '_WT.fig'], 'fig')
        saveas( fB2, [ffig '_KO.jpg'], 'jpg')
        saveas( fB2, [ffig '_KO.fig'], 'fig')
        ffig = [ dirLoc dirFig 'FirstReboundSpike'  tmpTxt];
        saveas(fRB, [ffig '.fig'] , 'fig');
        saveas(fRB, [ffig '.jpg'] , 'jpg');
        ffig = [ dirLoc dirFig 'AllBurstSpike'  tmpTxt]; % Better minus by expected baseline activity
        saveas(fBspk, [ffig '.fig'] , 'fig');
        saveas(fBspk, [ffig '.jpg'] , 'jpg');
        
        if (Close_Fig_aftr_save)
            close(fLA1); close(fLA2); close(fB1); close(fB2); close(fRB); close(fBspk);
        end
        
    end
end

cnt =0;
p1_ii = 1; p2_ii=1;
for p3_ii = 1 : length(PARAM3)
    Bfg = figure;  set(Bfg,'position',[680   283   712   695]); set(gcf,'PaperPositionMode','auto');
    Bfg2 = figure;  set(Bfg2,'position',[680   283   712   695]);  set(gcf,'PaperPositionMode','auto');
     nR = length(PARAM3); nC = length(PARAM4);
    for p4_ii = 1 : length(PARAM4)
        cnt = cnt+1;
        tmpTxt = get_Parameters_titleText(PARAMETERS, [4], [ p4_ii ]);
        figure(Bfg)
        subplot(nR,nC, cnt); hold on;
        plot( PARAM5, squeeze(CntBurstSpkWT(p1_ii,p2_ii,p3_ii,p4_ii,:)),'*-k');
        plot( PARAM5, squeeze(CntBurstSpkKO(p1_ii,p2_ii,p3_ii,p4_ii,:)),'*-r');
        legend('WT','KO','location','Best')
        title(tmpTxt);
        xlabel(titleTxt5); ylabel('Avg Burst spike')
        
        figure(Bfg2)
        subplot(nR,nC, cnt); hold on;
        plot( PARAM5, squeeze(ChanceOfBurstingWT(p1_ii,p2_ii,p3_ii,p4_ii,:)),'*-k');
        plot( PARAM5, squeeze(ChanceOfBurstingKO(p1_ii,p2_ii,p3_ii,p4_ii,:)),'*-r');
         ylim([ 0 0.7])
        legend('WT','KO','location','Best')
        title(tmpTxt);
        xlabel(titleTxt5); ylabel('Chance of Burst spike')
        
    end
     tmpTxt = get_Parameters_titleText(PARAMETERS, [3], [ p3_ii ]);
    figure(Bfg)
    suptitle(tmpTxt);
    figure(Bfg2)
    suptitle([ 'Chance of Bursting, '  tmpTxt]);
    
    if (SAVE_FIG)
        tmpTxt = get_Parameters_saveText(PARAMETERS, [3], [ p3_ii ]);
        ffig = [ dirLoc dirFig 'NumBurstSpk' tmpTxt ];
        saveas(  Bfg, [ffig '.jpg'], 'jpg');        saveas(  Bfg, [ffig '.fig'], 'fig');
        ffig = [ dirLoc dirFig 'ChanceOfBurst' tmpTxt];
        saveas(  Bfg2, [ffig '.jpg'], 'jpg');        saveas(  Bfg2, [ffig '.fig'], 'fig');
        if (Close_Fig_aftr_save)
            close(Bfg); close(Bfg2);
        end
    end
end

%%
%Make only the neuron activity during burst of VL layer
Close_Fig_aftr_save = 0;
TimeofFirstReboundSpk = zeros(ACT_Rec_size);
p1_ii = 1; p2_ii=1;
for p3_ii = 1 : length(PARAM3)
    fRB = figure; set(fRB,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
    fBspk = figure; set(fBspk,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
    fcum = figure; set(fcum,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
    fcumspk = figure; set(fcumspk,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
    
    cnt = 0;
    nR = length(PARAM4); nC = length(PARAM5);
    if(nR == 1)&&(nC ==1)
         set(fRB,'position',[ 1228         336         592         398]);
         set(fBspk,'position',[ 1228         336         592         398]);
         set(fcum,'position',[ 1228         336         592         398]);
         set(fcumspk,'position',[ 1228         336         592         398]);
    end
    for p4_ii = 1 : length(PARAM4)
        for p5_ii = 1 : length(PARAM5)
            la_ii = p3_ii; m_ii = p4_ii; s_ii = p5_ii;
            
            cnt = cnt +1;
            BasalAct = ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii}.VL;
            hbin = 2 ;
            histBin = 1:hbin:size(BasalAct.WT.All.BurstSpkTrain,2); 
          


            % Distribution of spikes in bursting range 
           
             tmpRebound_WT = sum(BasalAct.WT.All.BurstSpkTrain);
             tmpRebound_KO = sum(BasalAct.KO.All.BurstSpkTrain);
             
            ll =length(tmpRebound_WT);
            if mod(ll,hbin) == 0
                % Can use reshape fn.
                tmpRebound_WT = sum(reshape(tmpRebound_WT,[hbin ll/hbin]));
                tmpRebound_KO = sum(reshape(tmpRebound_KO,[hbin ll/hbin]));
                xx = histBin;
            else %If not, just set bin to 1
                xx = 1:1:size(BasalAct.KO.All.BurstSpkTrain,2);
            end
            
            figure(fBspk); 
            subplot(nR,nC,cnt);
            pH=bar(xx,tmpRebound_WT/ncells*100,'k'); hold on; 
            cH = get(pH,'Children');
            set(cH,'FaceAlpha',0.5); % 0 = transparent, 1 = opaque.
            pH=bar(xx,tmpRebound_KO/ncells*100,'FaceColor','r','EdgeColor','r');
           cH = get(pH,'Children');
            set(cH,'FaceAlpha',0.5); % 0 = transparent, 1 = opaque.
            xlim([0 BurstRange+10]);
            ylabel('% Neurons'); xlabel('Time(ms)')
            title([ titleTxt4 ' = ' num2str(PARAM4(p4_ii)) ', ' titleTxt5 ' = ' num2str(PARAM5(p5_ii))])
            legend('WT','KO','location','best')
             
            NWT = tmpRebound_WT; NKO = tmpRebound_KO;
            figure(fcumspk); 
            subplot(nR,nC,cnt);
            plot(xx,cumsum(NWT)./sum(NWT),'k'); hold on;
            plot(xx,cumsum(NKO)./sum(NKO),'r'); hold on;
            xlim([0 BurstRange+10]);
            xlabel('Latency of Spikes in bursting range'); ylabel('Cumulative Probability');
           legend('WT','KO','location','best')
           
              % Distribution of first rebound spike timing 
            tmpRebound_WT = zeros(ncells,1); tmpRebound_KO = zeros(ncells,1);
            for id = 1 : ncells
                tmp = find(BasalAct.WT.All.BurstSpkTrain(id,:),1,'first');
                if ~isempty(tmp) 
                    tmpRebound_WT(id) = tmp;                  
                end
                
                tmp = find(BasalAct.KO.All.BurstSpkTrain(id,:),1,'first');
                if ~isempty(tmp) 
                    tmpRebound_KO(id) = tmp;                  
                end                
            end      
             tmpRebound_WT = tmpRebound_WT(tmpRebound_WT ~= 0);
             tmpRebound_KO = tmpRebound_KO(tmpRebound_KO ~= 0);
            figure(fRB); 
            
             histBin = xx;
            subplot(nR,nC,cnt);
            [NWT, X] = hist(tmpRebound_WT,histBin);
            pH= bar(X,NWT/ncells*100,'FaceColor','k','EdgeColor','k'); hold on; 
            cH = get(pH,'Children');
            set(cH,'FaceAlpha',0.5); % 0 = transparent, 1 = opaque.

            [NKO,X] = hist(tmpRebound_KO,histBin);
            pH=bar(X,NKO/ncells*100,'FaceColor','r','EdgeColor','r'); hold on; 
            cH = get(pH,'Children');
            set(cH,'FaceAlpha',0.5); % 0 = transparent, 1 = opaque.
            xlim([0 BurstRange+10]);
            ylabel('% Neurons'); xlabel('Time(ms)')
            title([ titleTxt4 ' = ' num2str(PARAM4(p4_ii)) ', ' titleTxt5 ' = ' num2str(PARAM5(p5_ii))])
             legend('WT','KO','location','best')

             
            figure(fcum); 
            subplot(nR,nC,cnt);            
            plot(X, cumsum(NWT)./sum(NWT),'k'); hold on;
             plot(X, cumsum(NKO)./sum(NKO),'r'); hold on; 
             xlim([0 BurstRange+10]);
             xlabel('First rebound spike timing'); ylabel('Cumulative Probability');
         legend('WT','KO','location','best')
          
            clear BasalAct spkBin
        end
    end
    figure(fRB); suptitle(['Latency of first Rebound spike : ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]);  legend('WT','KO','location','best') % Note equal to the latency of peak firing rate --> latency of peak firing rate calculate from instantaneous firing rate
    figure(fcum); suptitle(['Latency of first Rebound spike : ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]); legend('WT','KO','location','best')
    figure(fBspk); suptitle(['Percent of neuron that fire during bursting period : ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]); legend('WT','KO','location','best') % Note equal to the latency of peak firing rate --> latency of peak firing rate calculate from instantaneous firing rate
    figure(fcumspk); suptitle(['Percent of neuron that fire during bursting period : ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]); legend('WT','KO','location','best') %

end

if (SAVE_FIG)
        tmpTxt = get_Parameters_saveText(PARAMETERS, [3], [ p3_ii ]);
        ffig = [ dirLoc dirFig 'FirstReboundSpike'  tmpTxt];
        saveas(fRB, [ffig '.fig'] , 'fig');
        saveas(fRB, [ffig '.jpg'] , 'jpg');
        ffig = [ dirLoc dirFig 'AllBurstSpike'  tmpTxt]; % Better minus by expected baseline activity
        saveas(fBspk, [ffig '.fig'] , 'fig');
        saveas(fBspk, [ffig '.jpg'] , 'jpg');
        ffig = [ dirLoc dirFig 'FirstReboundSpikeCumSum'  tmpTxt];
        saveas(fcum, [ffig '.fig'] , 'fig');
        saveas(fcum, [ffig '.jpg'] , 'jpg');
        ffig = [ dirLoc dirFig 'AllBurstSpikeCumSum'  tmpTxt]; % Better minus by expected baseline activity
        saveas(fcumspk, [ffig '.fig'] , 'fig');
        saveas(fcumspk, [ffig '.jpg'] , 'jpg');
        
        
        if (Close_Fig_aftr_save)
            close(fRB); close(fBspk); close(fcum); close(fcumspk);
        end 
        
end

%%
% for Poster , only one case    

%Make only the neuron activity during burst of VL layer
Close_Fig_aftr_save = 0;
TimeofFirstReboundSpk = zeros(ACT_Rec_size);
p1_ii = 1; p2_ii=1; p3_ii = 1; p4_ii=1; p5_ii = 1;   
         set(fRB,'position',[ 1228         336         592         398]);
         set(fBspk,'position',[ 1228         336         592         398]);
         set(fcum,'position',[ 1228         336         592         398]);
         set(fcumspk,'position',[ 1228         336         592         398]);

            la_ii = p3_ii; m_ii = p4_ii; s_ii = p5_ii;
            

            BasalAct = ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii}.VL;
            hbin = 2 ;
            histBin = 1:hbin:size(BasalAct.WT.All.BurstSpkTrain,2); 
          


            % Distribution of spikes in bursting range 
           
             tmpRebound_WT = sum(BasalAct.WT.All.BurstSpkTrain);
             tmpRebound_KO = sum(BasalAct.KO.All.BurstSpkTrain);
             
            ll =length(tmpRebound_WT);
            if mod(ll,hbin) == 0
                % Can use reshape fn.
                tmpRebound_WT = sum(reshape(tmpRebound_WT,[hbin ll/hbin]));
                tmpRebound_KO = sum(reshape(tmpRebound_KO,[hbin ll/hbin]));
                xx = histBin;
            else %If not, just set bin to 1
                xx = 1:1:size(BasalAct.KO.All.BurstSpkTrain,2);
            end
            
            figure(fBspk); 
           
            pH=bar(xx,tmpRebound_WT/ncells*100,'k'); hold on; 
            cH = get(pH,'Children');
            set(cH,'FaceAlpha',0.5); % 0 = transparent, 1 = opaque.
            pH=bar(xx,tmpRebound_KO/ncells*100,'FaceColor','r','EdgeColor','r');
           cH = get(pH,'Children');
            set(cH,'FaceAlpha',0.5); % 0 = transparent, 1 = opaque.
            xlim([0 BurstRange+10]);
            ylabel('% Neurons'); xlabel('Time(ms)')
            title([ titleTxt4 ' = ' num2str(PARAM4(p4_ii)) ', ' titleTxt5 ' = ' num2str(PARAM5(p5_ii))])
            legend('WT','KO','location','best')
             
            NWT = tmpRebound_WT; NKO = tmpRebound_KO;
            figure(fcumspk); 
  
            plot(xx,cumsum(NWT)./sum(NWT),'k'); hold on;
            plot(xx,cumsum(NKO)./sum(NKO),'r'); hold on;
            xlim([0 BurstRange+10]);
            xlabel('Latency of Spikes in bursting range'); ylabel('Cumulative Probability');
           legend('WT','KO','location','best')
           
              % Distribution of first rebound spike timing 
            tmpRebound_WT = zeros(ncells,1); tmpRebound_KO = zeros(ncells,1);
            for id = 1 : ncells
                tmp = find(BasalAct.WT.All.BurstSpkTrain(id,:),1,'first');
                if ~isempty(tmp) 
                    tmpRebound_WT(id) = tmp;                  
                end
                
                tmp = find(BasalAct.KO.All.BurstSpkTrain(id,:),1,'first');
                if ~isempty(tmp) 
                    tmpRebound_KO(id) = tmp;                  
                end                
            end      
             tmpRebound_WT = tmpRebound_WT(tmpRebound_WT ~= 0);
             tmpRebound_KO = tmpRebound_KO(tmpRebound_KO ~= 0);
            figure(fRB);             
             histBin = xx;
          
            [NWT, X] = hist(tmpRebound_WT,histBin);
            pH= bar(X,NWT/ncells*100,'FaceColor','k','EdgeColor','k'); hold on; 
            cH = get(pH,'Children');
            set(cH,'FaceAlpha',0.5); % 0 = transparent, 1 = opaque.

            [NKO,X] = hist(tmpRebound_KO,histBin);
            pH=bar(X,NKO/ncells*100,'FaceColor','r','EdgeColor','r'); hold on; 
            cH = get(pH,'Children');
            set(cH,'FaceAlpha',0.5); % 0 = transparent, 1 = opaque.
            xlim([0 BurstRange+10]);
            ylabel('% Neurons'); xlabel('Time(ms)')
            title([ titleTxt4 ' = ' num2str(PARAM4(p4_ii)) ', ' titleTxt5 ' = ' num2str(PARAM5(p5_ii))])
             legend('WT','KO','location','best')

             
            figure(fcum); 
    
            plot(X, cumsum(NWT)./sum(NWT),'k'); hold on;
             plot(X, cumsum(NKO)./sum(NKO),'r'); hold on; 
             xlim([0 BurstRange+10]);
             xlabel('First rebound spike timing'); ylabel('Cumulative Probability');
         legend('WT','KO','location','best')


            
%%
if(0)
    % Zoom in and examine Slope
    
    fg_Allcut = figure;  set(fg_Allcut, 'position', [  221         443        1606         475]);
    avgFR_WT_All = sWT_All_VL/ncells*1000;
    avgFR_KO_All = sKO_All/ncells*1000;
    plot(avgFR_WT_All,'k','LineWidth',3);  hold on;
    plot(avgFR_KO_All,'m','LineWidth',3);  hold on;
    legend('WT','KO')
    % title(['All cell, Inject I : amp = ' num2str(Input_I_amp) ', dur = ' num2str(Input_I_dur) ])
    ylabel('Average firing rate (Hz)')
    xlabel('Time(ms)')
    zoomx1 = 100; zoomx2 = 2000;
    xlim([zoomx1 zoomx2])
    ylim( [0 max([ avgFR_WT_All avgFR_WT_All])+2])
    legend('WT','KO')
    
    % PeakAvgFr_WT_All = max(avgFR_WT_All(zoomx1:zoomx2));
    [pks, locs] = findpeaks(avgFR_WT_All); pks =  pks(locs>InjStopT); locs = locs(locs > InjStopT); locs = locs(1); pks =  pks(1); % Find the first peak
    peakvalWT_All = pks; peakTimeWT_All = locs;
    [pks, locs] = findpeaks(avgFR_KO_All);  pks =  pks(locs>InjStopT); locs = locs(locs > InjStopT); locs = locs(1); pks =  pks(1); % Find the first peak
    peakvalKO_All = pks; peakTimeKO_All = locs;
    % [peakvalWT_All, peakTimeWT_All] = max(avgFR_WT_All(InjStopT:InjStopT+burstRange)); % Consider finding after stop injection negative current only
    % [peakvalKO_All, peakTimeKO_All] = max(avgFR_KO_All(InjStopT:InjStopT+burstRange));
    % peakTimeWT_All = peakTimeWT_All + InjStopT ;
    % peakTimeKO_All = peakTimeKO_All + InjStopT ;
    % figure(17); hold on;
    scatter(peakTimeWT_All, peakvalWT_All,150,'*k','LineWidth',2); peakTimeWT_All = peakTimeWT_All-1; hold on; %plot start at time 0
    scatter(peakTimeKO_All, peakvalKO_All,150,'*m','LineWidth',2); peakTimeKO_All = peakTimeKO_All -1;
    rectangle('Position',[InjStartT , 0,Input_I_dur,  max([ avgFR_WT_All avgFR_WT_All])+2],'EdgeColor','r');
    set(gca,'box','off')
    title({['Average firing rate Fr(t) cell in WT and KO , Inject I : amp = ' num2str(Input_I_amp) ', dur = ' num2str(Input_I_dur)],['WT rebound peak = ' num2str(peakvalWT_All) ' Hz at ' num2str(peakTimeWT_All) ' ms , KO rebound peak = ' num2str(peakvalKO_All) ' Hz at ' num2str(peakTimeKO_All) ' ms']});
    %% Raster Plot and Fr(t)
    raster_from_spkbin2(  WT.E.spktrain, WT.I.spktrain,avgFR_CUTTIME, Tstop,[ 'WT, Inject I : amp = ' num2str(Input_I_amp) ', dur = ' num2str(Input_I_dur) ],RES, SIGMA_fr_E,SIGMA_fr_I,SIGMA_fr_All)
    raster_from_spkbin2(  KO.E.spktrain, KO.I.spktrain,avgFR_CUTTIME, Tstop, ['KO, Inject I : amp = ' num2str(Input_I_amp) ', dur = ' num2str(Input_I_dur) ],RES, SIGMA_fr_E,SIGMA_fr_I,SIGMA_fr_All)
    %% save
    ReboundPeakAmp_WTlist(A_ii, D_ii) = peakvalWT_All ; ReboundPeakDel_WTlist(A_ii, D_ii) = peakTimeWT_All;
    ReboundPeakAmp_KOlist(A_ii, D_ii) = peakvalKO_All; ReboundPeakDel_KOlist(A_ii, D_ii) = peakTimeKO_All;
    Basal_Act.WT.avgFR = avgFR_WT_All;
    Basal_Act.KO.avgFR = avgFR_KO_All;
    
end
%%
%%
if(SAVE_WORKSPACE)
    save([dirLoc dirFig 'Saved_Workspace_result_' SPECIFIED_BASAL_ACT_CODENAME date '.mat' ],'-v7.3');
end