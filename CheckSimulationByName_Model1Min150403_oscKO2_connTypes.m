
clc
close all
clear all
WT = [];
KO = [];
tmp = [];
%NoiseSTDEV_List = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
RES = 1; %1 bin = 1 ms

% Controlling saveing
SAVE_WORKSPACE = 1;
SAVE_BASAL_ACT = 1; SPECIFIED_BASAL_ACT_CODENAME = '';
SAVE_FIG = 1; Close_Fig_aftr_save = 1;
% Statistical - Test
DoTTest = 1;
%raster plot
plotSampleV = 1;


% Simulation setting
avgFR_CUTTIME = 500;
% PhotoInjT = 1500;
% InjStopT = PhotoInjT+Input_I_dur;
% DelayT = 0;
% BurstRange = 100;
% CUTTIME = avgFR_CUTTIME;
TSTOP = 3000;

CUTTIME = 500;
PhotoInjT = 500;% 1500;
PhotoStop = 500;
DelayT = 0;
BurstRange = 0;

% Parameters setting
NUM_TRIAL = 5;
ncells = 1150;
M1_ncells = 166;
TRIAL_LST = 1 : NUM_TRIAL;

rTC_LST = [50 100 150 200 250]; % 
wmTC_LST = [50 100 200 300 400 500]; %[10 25 50 75 100];

LightAmp_LST = [0.5];
LightDur_LST = [1000];
GPmVLw_mean_LST = [ 0.5];
GPmVLw_sig_LST =[0];

OSC_F_LST = [20 40 ];
OSC_Amp_LST = [0 0.5 1];
OSC_phase_LST = 0;

Input_FR_LST = 50; %[50 100 250 450];

TRIAL_NO_LST = 6;

PARAM1 = rTC_LST;
lblTxt1 = 'Range of thalamocortical connection';
saveTxt1 = 'rTC';
titleTxt1 = 'Range_T_C';
PARAM2 = wmTC_LST;
lblTxt2 = 'Weight of thalamocortical connection';
saveTxt2 = 'wmTC';
titleTxt2 = 'W_T_C';
PARAM3 = Input_FR_LST ;
lblTxt3 = 'Input Freq';
saveTxt3 = 'InputFR';
titleTxt3 = 'Input Frequency';
PARAM4 = OSC_F_LST;
lblTxt4 = 'Oscillation Frequency';
saveTxt4 = 'OSC_F';
titleTxt4 = 'Osc Freq';
PARAM5 = OSC_Amp_LST;
lblTxt5 = 'Oscillation Amplitude relative to mean FR';
saveTxt5 = 'OSC_amp';
titleTxt5 = 'Osc Amp';



% PARAM5 = OSC_phase_LST; % TRIAL_LST;
% lblTxt5 = 'Phase of oscillating input (Hz)';
% saveTxt5 = 'OscPhase'; %'trial';
% titleTxt5 = 'Osc Phase';

% PARAM4 = LightAmp_LST; % TRIAL_LST;
% lblTxt4 = 'Light Stimulus amplitude (nA)'; %'Trial'; %
% saveTxt4 = 'GPmAmp'; %'trial';
% titleTxt4 = 'Light Amp';
% PARAM5 = LightDur_LST; % TRIAL_LST;
% lblTxt5 = 'Light Stimulus duration (ms)'; %'Trial'; %
% saveTxt5 = 'GPmDur'; %'trial';
% titleTxt5 = 'Light Dur';

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

% Directory


PATH = SetPath;
dirLoc = [PATH  'OscInput_varyTCtype_Sim/'];  % 'OscInput_Sim/']; 
% dirFig = ['Fig' get_Parameters_RangeTxt( PARAMETERS,[1,2,4,5]) '/'];
% mkdir([dirLoc dirFig])


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
                    la_ii = 1;  m_ii = 1; s_ii = 1; ld_ii= 1;   
                    of_ii = p4_ii; oa_ii = p5_ii; op_ii = 1;
                    TRIAL_NO = 5;
                    fr_ii = p3_ii;
                    cell_type = 1;
%                     for cell_type = 1 : 2
                        
                        if (cell_type == 1)
                            cTxt = 'WT';
                        elseif (cell_type == 2)
                            cTxt = 'KO';
                        end
                        
                        % PDfewBurst_GPmVLmd1_rTC120_wmTC10_WT_GPmInput_Amp0.3_Dur1000_GPmVLw_m0.06_sig0.01_InGauss0.2_IGmean-0.15_IGmeanSig0_W0.0015_SpecifiedPoisSpk_sig0.00Hz_T4000_trial3
                        
                        coreFileName =  'GPmVLmd1_0del_KO2_avgPnegExpWTC';% 'GPmVLmd1_0del_KO2_avgPWuniformTC';%  % GPmVLmd1_0del_KO2' ;
                        
                        InGauss_STDEV = 0; %0.2;, 0.3
                        NoiseMEAN = 0;
                        IGmeanSig = 0;
                        W_Weight = 0.001;
                        PoisInputFr = Input_FR_LST(fr_ii);
                        TSTOP = 3000;
%                         GPmLightDur = LightDur_LST(ld_ii);
                        
                        rTC =rTC_LST(r_ii);
                        wmTC = wmTC_LST(wm_ii);
                        osc_f = OSC_F_LST(of_ii);
                        osc_amp =OSC_Amp_LST(oa_ii);
                        osc_phase = OSC_phase_LST(op_ii);
                                                
                        
%                         GPmLight = LightAmp_LST(la_ii);
%                         GPm_w_mn = GPmVLw_mean_LST(m_ii);
%                         GPm_w_sg = GPmVLw_sig_LST(s_ii);
                        
                        txtFR = sprintf('%2.2f',PoisInputFr); txtAmp = sprintf('%2.2f',osc_amp);
                        if( InGauss_STDEV ==0)
                        Simulation_Code = [coreFileName '_rTC' num2str(rTC) '_wmTC' num2str(wmTC) '_' cTxt '_IGmean' num2str(NoiseMEAN) '_IGmeanSig' num2str(IGmeanSig) ...
                            '_W' num2str(W_Weight) '_' txtFR 'Hz_oscF' num2str(osc_f) 'Hz_amp' txtAmp '_phase' num2str(osc_phase) '_T' num2str(TSTOP) '_trial' num2str(TRIAL_NO)];
                        else
                                                    Simulation_Code = [coreFileName '_rTC' num2str(rTC) '_wmTC' num2str(wmTC) '_' cTxt '_InGauss' num2str(InGauss_STDEV) '_IGmean' num2str(NoiseMEAN) '_IGmeanSig' num2str(IGmeanSig) ...
                            '_W' num2str(W_Weight) '_' txtFR 'Hz_oscF' num2str(osc_f) 'Hz_amp' txtAmp '_phase' num2str(osc_phase) '_T' num2str(TSTOP) '_trial' num2str(TRIAL_NO)];
                        end
                       
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


