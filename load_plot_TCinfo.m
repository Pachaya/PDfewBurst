
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

SAVE_WORKSPACE = 0;

DoTTest = 0; DoSampleTtest = 0;
SAVE_BASAL_ACT = 1;
plotSampleV = 0;
SPECIFIED_RATIO = 0; %TestSim_WT_InGauss0.2_W0.0005_50.00Hz_rEE0_rEI1_rIE2_rII1_Wmult10


NoiseMEAN = 0;

SAVE_FIG = 0; Close_Fig_aftr_save = 0;

dirLoc = 'Model1_MinParam150403_testTC/' ;
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
dirFig = 'Test_sumTCweight_specified_wmTC/';
mkdir([dirLoc dirFig])
TSTOP = 40;

rTC_LST = [10:10:300];
wmTC_LST = [1];

specified_wmTC_LST = [0.001454];


LightAmp_LST = [0.5; ];
GPmVLw_mean_LST = [0;];
GPmVLw_sig_LST =[ 0.5 ];

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

CUTTIME = 0;
PhotoInjT = 1;
PhotoStop = 1;
DelayT = 0;
BurstRange = 40;

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
                        
                        
                        coreFileName = 'PDfewBurst_GPmVLmd1' ;
                        
                        InGauss_STDEV = 0.2; %0.2;, 0.3
                        NoiseMEAN = 0;
                        IGmeanSig = 0;
                        W_Weight = 0.029;
                        PoisInputFr = 0;
                        TSTOP = 40;
                        GPmLightDur = 1000;
                        
                        rTC = rTC_LST(r_ii);
%                         wmTC = wmTC_LST(wm_ii);
                        specified_wmTC = specified_wmTC_LST(swm_ii);
                        
                        GPmLight = LightAmp_LST(la_ii);
                        GPm_w_mn = GPmVLw_mean_LST(m_ii);
                        GPm_w_sg = GPmVLw_sig_LST(s_ii);
                        
                        
                        Simulation_Code = [coreFileName '_rTC' num2str(rTC) '_specified_wmTC' num2str(specified_wmTC) '_' cTxt '_' 'GPmInput_Amp' num2str(GPmLight) '_Dur' num2str(GPmLightDur) '_GPmVLw_m' num2str(GPm_w_mn) '_sig' num2str(GPm_w_sg) ...
                            '_InGauss' num2str(InGauss_STDEV) '_IGmean' num2str(NoiseMEAN) '_IGmeanSig' num2str(IGmeanSig) '_W' num2str(W_Weight) '_SpecifiedPoisSpk_sig0.00Hz_T' num2str(TSTOP) '_trial' num2str(TRIAL_NO)];
                        % PDfewBurst_GPmVLmd1_0del_rTC300_specified_wmTC0.001454_KO_GPmInput_Amp0.5_Dur1000_GPmVLw_m0_sig0.5_InGauss0.2_IGmean0_IGmeanSig0_W0.029_SpecifiedPoisSpk_sig0.00Hz_T40_trial1

                        disp('==================================================================================================')
                        disp(Simulation_Code)
                        disp('######  Download M1 ')
                        %M1
                        Name_postfix = [ 'M1_' Simulation_Code];
                        getCenterPart =1;
                        BaselineAct % Get all cell activity
                        
                        tmp.E = E;      tmp.I = I;  tmp.All = All; %tmp.VL = VL; tmp.M1 = M1;
                        if (cell_type == 1)
                            M1.WT = tmp;
                            M1.WT.Simulation_Code = Simulation_Code;
                            cTxt = 'WT';
                        elseif (cell_type == 2)
                            M1.KO = tmp;
                            M1.KO.Simulation_Code = Simulation_Code;
                            cTxt = 'KO';
                        end
                        
                        
                        figNameCode = [ saveTxt1 num2str(PARAM1(p1_ii))  '_' saveTxt2 num2str(PARAM2(p2_ii)) '_' cTxt];
                        
                        %                 CheckFileExist( dirLoc, Name_postfix  )
                        disp('######  Download VL ')
                        Name_postfix = [ Simulation_Code];
                        VL_PhotoactivationAll
                        
                        %                 if(0)
                        
                        tmp.E = E;      tmp.I = I;  tmp.All = All; %tmp.VL = VL; tmp.M1 = M1;
                        clear E I All somaVall somaVallE somaVallI
                        if (cell_type == 1)
                            VL.WT = tmp;
                            cTxt = 'WT'; 
                            Name_postfix_WT = Simulation_Code;
                            VL.WT.Simulation_Code = Simulation_Code;
                            clear tmp
                        elseif (cell_type == 2)
                            VL.KO = tmp;
                            cTxt = 'KO'; 
                            Name_postfix_KO = Simulation_Code;
                            VL.KO.Simulation_Code = Simulation_Code;
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
                        
                        disp('Finish plot sample membrain voltage');
                        clear somaVall
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
    save([dirLoc dirFig 'Saved_Activity_result.mat' ], 'ACT_Record','PARAMETERS','-v7.3');
    disp('Saved Activity Result')
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   VL - M1  Connection
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Count number of VL per M1  - find mean 
%  Record Weight of Connection

% rTC_LST = [120 130 140 150 170 190 210];
% wmTC_LST = [7];

% For 2D gaussian area under the curve is N?��?SQRT(2��)


p3_ii = 1; p4_ii = 1; p5_ii =1;
ssize = 1500; dist = ssize/2; pmax =0.85; W_scale = 1E-05;
for p1_ii = 1 : length(PARAM1)
    for p2_ii = 1 : length(PARAM2)
        disp('##################################################################################################')
        disp( [lblTxt1 ' = ' num2str(PARAM1(p1_ii)) ]);
        disp( [lblTxt2 ' = ' num2str(PARAM2(p2_ii)) ]);
        rangeTC = PARAM1(p1_ii);      sigTC = rangeTC/sqrt(2); wmax = PARAM2(p2_ii) * W_scale; dd =1:1:dist;        
%         expectedWsum = sum( funcGauss_byDistance(pmax, dd, sigTC ).*funcGauss_byDistance(wmax, dd, sigTC ));
%         disp(['Expected summation of weight = ' num2str(expectedWsum) ]);
        
   Name_postfix_WT  = ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii }.VL.WT.Simulation_Code; 
   Name_postfix_KO  = ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii }.VL.KO.Simulation_Code; 
   display = 1;
  
    disp('--------------------------------------------------------------------------------------------------')
    disp('WT' );
    disp('--------------------------------------------------------------------------------------------------')
    [TC_basedOnM1_WT, TC_sumW_WT, TC_maxW_WT, TC_numVL_WT ]               = ExtractTC_info( dirLoc, Name_postfix_WT,  display);
    
    disp('--------------------------------------------------------------------------------------------------')
    disp('KO' );
    disp('--------------------------------------------------------------------------------------------------')
    [TC_basedOnM1_KO, TC_sumW_KO, TC_maxW_KO, TC_numVL_KO ]               = ExtractTC_info( dirLoc, Name_postfix_KO,  display);
    disp('##################################################################################################')
    end
end


%%    
    
    
p3_ii = 1; p4_ii = 1; p5_ii =1;
ssize = 1500; dist = ssize/2; pmax =0.85; W_scale = 1E-05;
disp('RangeTC    numVL/M1    expected sumW    sumW   ')
for p1_ii = 1 : length(PARAM1)
    for p2_ii = 1 : length(PARAM2)
%         disp('##################################################################################################')
%         disp( [lblTxt1 ' = ' num2str(PARAM1(p1_ii)) ]);
%         disp( [lblTxt2 ' = ' num2str(PARAM2(p2_ii)) ]);
        rangeTC = PARAM1(p1_ii);      sigTC = rangeTC/sqrt(2); 
        specified_wmTC = specified_wmTC_LST(p2_ii);
%         wmax = PARAM2(p2_ii) * W_scale; dd =1:1:dist;        
%         expectedWsum = sum( funcGauss_byDistance(pmax, dd, sigTC ).*funcGauss_byDistance(wmax, dd, sigTC ));
%         disp(['Expected summation of weight = ' num2str(expectedWsum) ]);
        
   Name_postfix_WT  = ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii }.VL.WT.Simulation_Code; 
%    Name_postfix_KO  = ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii }.VL.KO.Simulation_Code; 
   display = 0;
 
    [TC_basedOnM1_WT, TC_sumW_WT, TC_maxW_WT, TC_numVL_WT ]               = ExtractTC_info( dirLoc, Name_postfix_WT,  display);
    disp([ num2str(rangeTC) '        ' num2str(mean(TC_numVL_WT)) '        ' num2str(specified_wmTC) '        ' num2str(mean(TC_sumW_WT)) ])
    
    end
end

% ##################################################################################################
%{
RangeTC    numVL/M1    expected sumW    sumW   
10        1        4.9016e-05        9.9719e-06
20        1.0854        0.00010228        1.0737e-05
30        1.5217        0.00015555        1.4885e-05
40        2.3062        0.00020881        2.2171e-05
50        3.5758        0.00026208        3.3615e-05
60        5.006        0.00031535        4.5874e-05
70        6.8675        0.00036861        6.1126e-05
80        8.7831        0.00042188        7.577e-05
90        11.0301        0.00047514        9.1737e-05
100        13.7651        0.00052841        0.00010978
110        16.6687        0.00058167        0.00012755
120        19.6386        0.00063494        0.00014435
130        23.2711        0.00068821        0.000163
140        26.9578        0.00074147        0.00018047
150        31        0.00079474        0.00019764
160        35.1928        0.000848        0.00021442
170        39.4699        0.00090127        0.0002302
180        44.4157        0.00095454        0.00024614
190        49.3916        0.0010078        0.00026005
200        54.7892        0.0010611        0.0002743
210        60.3795        0.0011143        0.00028792
220        65.9157        0.0011676        0.00029992
230        71.8614        0.0012209        0.00031154
240        78.1205        0.0012741        0.00032324
250        84.7771        0.0013274        0.00033395
260        91.6867        0.0013807        0.00034391
270        99.2108        0.0014339        0.00035359
280        106.2892        0.0014872        0.00036259
290        114.3614        0.0015405        0.00037106
300        122.6084        0.0015937        0.00037974
%}