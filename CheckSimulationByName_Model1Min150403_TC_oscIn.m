
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

DoTTest = 1;
SAVE_BASAL_ACT = 0;
FNAME_SIM = 'GPmVL' ; %'PreLim_INOISE0.2_InSPK50_EE1_EI2_IE20_II12_mult125_w0.0006';
plotSampleV = 0;
SPECIFIED_RATIO = 0; %TestSim_WT_InGauss0.2_W0.0005_50.00Hz_rEE0_rEI1_rIE2_rII1_Wmult10
Basal_Act_set  = [];

NoiseMEAN = 0;

SAVE_FIG = 1;

dirLoc = 'Test_Osc_Input/' ; 
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
dirFig = '';

% PDfewBurst_GPmVLmd1_0del_rTC150_wmTC20_WT_GPmInput_Amp0.6_Dur500_GPmVLw_m0.5_sig0_InGauss0.2_IGmean0_IGmeanSig0_W0.029_SpecifiedPoisSpk_sig0.00Hz_T4000_trial1
% ### Sorry ###

TSTOP = 4000;

%SomaVolt_M1_PDfewBurst_GPmVLmd1_rTC120_wmTC50_KO_GPmInput_Amp0.5_Dur1000_GPmVLw_m0.3_sig0_InGauss0.2_IGmean0_IGmeanSig0_W0.029_SpecifiedPoisSpk_sig0.00Hz_T4000_trial1
rTC_LST = [ 25 50 75 100 150 200 250 ]; %
wmTC_LST = [20 30 40 50 ]; %[10 25 50 75 100];


LightAmp_LST = [0.5];
LightDur_LST = [1000];
GPmVLw_mean_LST = [ 0.5];    
GPmVLw_sig_LST =[0];
% 

% Candidates Model 1 :
% I amp = 0.3, mean 0.06, sigma = 0.01, 0.005 ?0.03
% I amp = 0.5, mean 0.04, sigma = 0.03, 0.01 ?0.03
%PDfewBurst_GPmVLmd1_rTC120_wmTC10_WT_GPmInput_Amp0.3_Dur1000_GPmVLw_m0.06_sig0.01_InGauss0.2_IGmean-0.15_IGmeanSig0_W0.0015_SpecifiedPoisSpk_sig0.00Hz_T4000_trial3
NUM_TRIAL = 5;
ncells = 1150;
TRIAL_LST = 1 : NUM_TRIAL;


OSC_F_LST = [20 40];
OSC_Amp_LST = [1 1/2 1/3];


PARAM1 = rTC_LST;
lblTxt1 = 'Range of thalamocortical connection';
saveTxt1 = 'rTC';
titleTxt1 = 'Range_T_C';
PARAM2 = wmTC_LST;
lblTxt2 = 'Weight of thalamocortical connection';
saveTxt2 = 'wmTC';
titleTxt2 = 'W_T_C';
PARAM3 = LightAmp_LST; % TRIAL_LST;
lblTxt3 = 'Light Stimulus amplitude (nA)'; %'Trial'; %
saveTxt3 = 'GPmAmp'; %'trial';
titleTxt3 = 'Light Amp';
% PARAM4 = GPmVLw_mean_LST;
% lblTxt4 = 'GPm-VL weight mean';
% saveTxt4 = 'GPmVLw_m';
% titleTxt4 = 'W_m_e_a_n';
PARAM4 = OSC_F_LST;
lblTxt4 = 'Oscillation Frequency';
saveTxt4 = 'OSC_F';
titleTxt4 = 'Osc Freq';
PARAM5 = OSC_Amp_LST;
lblTxt5 = 'Oscillation Amplitude relative to mean FR';
saveTxt5 = 'OSC_amp';
titleTxt5 = 'Osc Amp';

% PARAM5 = GPmVLw_sig_LST;
% lblTxt5 = 'GPm-VL weight sigma';
% saveTxt5 = 'GPmVLw_sig';
% titleTxt5 = 'W_s_i_g';


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
Check_Status = zeros(ACT_Rec_size);

CUTTIME = 500;
PhotoInjT = 1500;
PhotoStop = 2500;
DelayT = 0;
BurstRange = 100;

%CheckSimulationByName_Combine_GpmVL_VLM1\
disp('==================================================================================================')
disp('============================    Check Simulation By Name    ======================================')
disp('==================================================================================================')
for p1_ii = 1 : length(PARAM1)
    for p2_ii = 1 : length(PARAM2)
        for p3_ii = 1 : length(PARAM3)
            for p4_ii = 1 : length(PARAM4)
                for p5_ii = 1 : length(PARAM5)
                    
                    r_ii = p1_ii; wm_ii = p2_ii;
                    la_ii = 1;  m_ii = p3_ii; s_ii = 1; ld_ii= 1; 
                    of_ii = p4_ii; oa_ii = p5_ii;                    
                    TRIAL_NO = 1;
                    cell_type = 1;
%                     for cell_type = 1 : 2
                        
                        if (cell_type == 1)
                            cTxt = 'WT';
                        elseif (cell_type == 2)
                            cTxt = 'KO';
                        end
                        
                            
                        coreFileName = 'GPmVLmd1_0del' ;
                        
                        InGauss_STDEV = 0.2; %0.2;, 0.3
                        NoiseMEAN = 0;
                        IGmeanSig = 0;
                        W_Weight = 0.029;
                        PoisInputFr = 0;
                        TSTOP = 4000;
                        GPmLightDur = LightDur_LST(ld_ii);
                        
                        rTC =rTC_LST(r_ii);
                        wmTC = wmTC_LST(wm_ii);
                                                
                        
                        GPmLight = LightAmp_LST(la_ii);
                        GPm_w_mn = GPmVLw_mean_LST(m_ii);
                        GPm_w_sg = GPmVLw_sig_LST(s_ii);
                        
                        osc_f = OSC_F_LST(of_ii);
                        osc_amp = OSC_Amp_LST(oa_ii);
                        phase = 0;
                        
                        tmpAmp = sprintf('%2.2f', osc_amp);
                        Simulation_Code = [coreFileName '_rTC' num2str(rTC) '_wmTC' num2str(wmTC) '_' cTxt '_' 'GPmInput_Amp' num2str(GPmLight) '_Dur' num2str(GPmLightDur) '_GPmVLw_m' num2str(GPm_w_mn) '_sig' num2str(GPm_w_sg) ...
                            '_InGauss' num2str(InGauss_STDEV) '_IGmean' num2str(NoiseMEAN) '_IGmeanSig' num2str(IGmeanSig) '_W' num2str(W_Weight) '_SpecifiedPoisSpk_sig0.00Hz_oscF'  num2str(osc_f) 'Hz_amp' tmpAmp '_phase0_T' num2str(TSTOP) '_trial' num2str(TRIAL_NO)];
                        % GPmVLmd1_0del_rTC25_wmTC30_KO_GPmInput_Amp0.5_Dur1000_GPmVLw_m0.5_sig0_InGauss0.2_IGmean0_IGmeanSig0_W0.029_SpecifiedPoisSpk_sig0.00Hz_oscF20Hz_amp1.00_phase0_T4000_trial1
                        
                        disp('==================================================================================================')
                        disp(Simulation_Code)
                        
                                                
                        figNameCode = [ saveTxt1 num2str(PARAM1(p1_ii))  '_' saveTxt2 num2str(PARAM2(p2_ii)) '_' cTxt];
                        
                        
                        
                        Name_postfix = [ Simulation_Code];
%                         Check_Status(p1_ii,p2_ii,p3_ii,p4_ii,p5_ii,cell_type) = CheckFileExist( dirLoc, Name_postfix  );
                        Check_Status(p1_ii,p2_ii,p3_ii,p4_ii,p5_ii) = CheckFileExist( dirLoc, Name_postfix  );
                        %VL_PhotoactivationAll
                        
                        
                        
%                     end
                end
            end
        end
    end
end

sum(Check_Status == 0)


