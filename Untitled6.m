
% clc
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

DoTTest = 1;
SAVE_BASAL_ACT = 0;
FNAME_SIM = 'PD4' ; %'PreLim_INOISE0.2_InSPK50_EE1_EI2_IE20_II12_mult125_w0.0006';
plotSampleV = 0;
SPECIFIED_RATIO = 0; %TestSim_WT_InGauss0.2_W0.0005_50.00Hz_rEE0_rEI1_rIE2_rII1_Wmult10
Basal_Act_set  = [];

NoiseMEAN = 0;

SAVE_FIG = 1;

dirLoc = 'GPmVL/' ; %'TestSim_BasalAct/';  %'VL_LocalConn/';
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
InjStartT = 1000;
InjStopT = InjStartT+Input_I_dur;
reboundDelay = 0;
burstRange = 100;
cutTime = avgFR_CUTTIME;

ADD_I_to_M1 = 0;
ADD_I_to_VL = 0;

FIG_ALL = 0;
dirFig = 'Fig/';
NoiseMEAN_WTKO = [0; 0;]; %[0.0291; 0.0075;]; %[0.0295; 0.016]; %%%%%%%%%%%%%%%%%%%%%%%%%%% [ WT; KO;] for 10 Hz -> [0.0295; 0.016]; for 5 Hz -> [0.0291; 0.0075;]; 5Hz 0.0318 0.3  for WT 5Hz (2) -> [0.0318; 0.0075;]
NoiseSIG_WTKO  =  [0.3; 0.3;]; %[0.8; 0.8;];  %[0.24; 0.3];  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  [WT; KO;] for 10 Hz -> [0.4; 0.5]; for 5 Hz -> [0.24; 0.3];  5Hz 0.0318 0.3  for WT 5Hz (2) -> [0.3; 0.3;]
VARIE_MEAN_InGauss = 1;
NoiseMEANsigma_WTKO = [0; 0;]; %%%%% [WT; KO;]
TSTOP = 3000;

rTC_LST = [120:10:150];
wmTC_LST = [7;];

LightAmp_LST = [0.3; 0.4; 0.5];
GPmVLw_mean_LST = [0.02:0.02:0.1];
GPmVLw_sig_LST =[0.005; 0.01; 0.02; 0.03];
% PDfewBurst_GPmVLmd1_KO_GPmInput_Amp0.4_Dur1000_GPmVLw_m0.1_sig0.03_InGauss0.2_IGmean0_IGmeanSig0_W0.001_0.00Hz_T3000
PARAM1 = LightAmp_LST;
legTxt1 = 'Light Stimulus amplitude (nA)'; %'E_L '; %'soma size';
saveTxt1 = 'GPmAmp';
PARAM2 = GPmVLw_mean_LST;
legTxt2 = 'GPm-VL weight mean';
saveTxt2 = 'GPmVLw_m';
PARAM2 = GPmVLw_sig_LST;
legTxt2 = 'GPm-VL weight sigma';
saveTxt2 = 'GPmVLw_sig';

BASAL_ACT_Record = cell(length(PARAM1),length(PARAM2),length(PARAM3));

for p1_ii = 1 : length(PARAM1)
    for p2_ii = 1 : length(PARAM2)
        for p3_ii = 1 : length(PARAM3)
            la_ii = p1_ii; m_ii = p2_ii; s_ii = p3_ii;
            
            for cell_type = 1 : 2
                
                if (cell_type == 1)
                    cTxt = 'WT';
                elseif (cell_type == 2)
                    cTxt = 'KO';
                end
                
                coreFileName = 'PDfewBurst_GPmVL' ;
                
                InGauss_STDEV = 0.2; %0.2;
                NoiseMEAN = 0;
                IGmeanSig = 0;
                W_Weight = 0.001;
                PoisInputFr = 0;
                TSTOP = 3000;
                GPmLightDur = 1000;
                
                GPmLight = LightAmp_LST(la_ii);
                GPm_w_mn = GPmVLw_mean_LST(m_ii);
                GPm_w_sg = GPmVLw_sig_LST(s_ii);
                
                
                Name_postfix = [coreFileName '_' cTxt '_' 'GPmInput_Amp' num2str(GPmLight) '_Dur' num2str(GPmLightDur) '_GPmVLw_m' num2str(GPm_w_mn) '_sig' num2str(GPm_w_sg) ...
                    '_InGauss' num2str(InGauss_STDEV) '_InGmean' num2str(NoiseMEAN) '_IGmeanSig' num2str(IGmeanSig) '_W' num2str(W_Weight) '_' num2str(PoisInputFr) '.00Hz_T' num2str(TSTOP) ];
                
                % PDfewBurst_GPmVLmd1_KO_GPmInput_Amp0.4_Dur1000_GPmVLw_m0.1_sig0.03_InGauss0.2_IGmean0_IGmeanSig0_W0.001_0.00Hz_T3000
                
                
                
                disp('==================================================================================================')
                disp(Name_postfix)
                if (cell_type == 1)
                    cTxt = 'WT';
                else
                    cTxt = 'KO';
                end
                figNameCode = [ saveTxt1 num2str(PARAM1(p1_ii))  '_' saveTxt2 num2str(PARAM2(p2_ii)) '_' cTxt];
                %             VL_Basal_ActivityCenter %call function
                %             VL_Basal_Activity
                getCenterPart_M1 = 0;
                Sim_Neuron_Activity_PD4
                
                
                
                tmp.E = E;      tmp.I = I;  tmp.All = All; tmp.VL = VL; tmp.M1 = M1;
                if (cell_type == 1)
                    WT = tmp;
                    cTxt = 'WT';
                elseif (cell_type == 2)
                    KO = tmp;
                    cTxt = 'KO';
                end
                if (plotSampleV)
                    %soma's volt
                    [NEcell,NIcell, Tstop,somaVall] = get_save_vec_file_BigNet([dirLoc 'SomaVolt_' Name_postfix '.txt'] );
                    samE = 385;
                    samI = 129;
                    samTstr = 1000; samTstp = 1800; %500+ Input_I_dur +200;
                    fgSmpl = figure; plot(samTstr:samTstp, somaVall(samE,samTstr:samTstp), 'r'); hold on;
                    plot(samTstr:samTstp, somaVall(samE+1,samTstr:samTstp), 'm');
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
            DoSampleTtest = 1;
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
            R
            BASAL_ACT_Record{p1_ii,p2_ii,p3_ii } = Basal_Act;
        end
    end
end
