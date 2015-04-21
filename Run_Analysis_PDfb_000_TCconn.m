
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
FNAME_SIM = 'GPmVL' ; %'PreLim_INOISE0.2_InSPK50_EE1_EI2_IE20_II12_mult125_w0.0006';
plotSampleV = 0;
SPECIFIED_RATIO = 0; %TestSim_WT_InGauss0.2_W0.0005_50.00Hz_rEE0_rEI1_rIE2_rII1_Wmult10
Basal_Act_set  = [];

NoiseMEAN = 0;

SAVE_FIG = 1;

dirLoc = 'TestTCconn/' ;
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
dirFig = 'Fig/';
TSTOP = 4000;

rTC_LST = [50:10:300];
wmTC_LST = [10;];

LightAmp_LST = [0.3; 0.4; 0.5];
GPmVLw_mean_LST = [0.02; 0.04;];
GPmVLw_sig_LST =[0.005; 0.01; 0.02; 0.03;];

% Candidates Model 1 :
% I amp = 0.3, mean 0.06, sigma = 0.01, 0.005 – 0.03
% I amp = 0.5, mean 0.04, sigma = 0.03, 0.01 – 0.03
%PDfewBurst_GPmVLmd1_rTC120_wmTC10_WT_GPmInput_Amp0.3_Dur1000_GPmVLw_m0.06_sig0.01_InGauss0.2_IGmean-0.15_IGmeanSig0_W0.0015_SpecifiedPoisSpk_sig0.00Hz_T4000_trial3
NUM_TRIAL = 5;
ncells = 1150;
TRIAL_LST = 1 : NUM_TRIAL;

PARAM1 = rTC_LST; 
PARAM2 = wmTC_LST; 
TC_Conn_List = cell(length(PARAM1), length(PARAM2),2);


for p1_ii = 1 : length(PARAM1)
    for p2_ii = 1 : length(PARAM2)
        
            r_ii = p1_ii; wm_ii = p2_ii;
            
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
                TSTOP = 30;
                GPmLightDur = 1000;
                TRIAL_NO = 1;
                rTC = rTC_LST(r_ii);
                wmTC = wmTC_LST(wm_ii);
%                 GPmLight = 0.3;
              
                
%                                 GPmLight = LightAmp_LST(la_ii);
%                                 GPm_w_mn = GPmVLw_mean_LST(m_ii);
%                                 GPm_w_sg = GPmVLw_sig_LST(s_ii);
                
                
                Simulation_Code = [coreFileName '_rTC' num2str(rTC) '_wmTC' num2str(wmTC) '_' cTxt ...
                    '_InGauss' num2str(InGauss_STDEV) '_W' num2str(W_Weight) '_T' num2str(TSTOP) '_trial' num2str(TRIAL_NO)];
                
                %  PDfewBurst_GPmVLmd1_rTC300_wmTC10_WT_InGauss0.2_W0.0012_T30_trial1
                disp('==================================================================================================')
                disp(Simulation_Code)
                 Name_postfix = [ Simulation_Code]; 
               
                 TC_Conn = Get_TCconn_Info(  Simulation_Code, dirLoc );
                
                TC_Conn_List{p1_ii,p2_ii,cell_type} =  TC_Conn;
                %   CheckFileExist( dirLoc, Name_postfix  )                
            
        end
    end
end

%% RangeTC vs # of VL per M1

VLperM1_LST = zeros(length(PARAM1),2);
VLperM1_std = zeros(length(PARAM1),2);
p2_ii = 1;
for p1_ii = 1 : length(PARAM1)
    for cell_type = 1 : 2
        VLperM1_LST(p1_ii,cell_type) = mean(TC_Conn_List{p1_ii,p2_ii,cell_type}.NumVL_of_M1); 
        VLperM1_std(p1_ii,cell_type) = std(TC_Conn_List{p1_ii,p2_ii,cell_type}.NumVL_of_M1); 
    end
end

figure; 
plot(rTC_LST,VLperM1_LST(:,1),'*-k');hold on;
plot(rTC_LST,VLperM1_LST(:,2),'*-r');
title('Number of VL per M1 in the simulation')        
xlabel('Range of TC connection'); ylabel('#VL / M1');
hold on;

for p1_ii = 1 : length(PARAM1)
    for cell_type = 1 : 2
        tmp = TC_Conn_List{p1_ii,p2_ii,cell_type}.NumVL_of_M1; 
        scatter(rTC_LST(p1_ii).*ones(length(tmp),1), tmp,'k')
        scatter(rTC_LST(p1_ii).*ones(length(tmp),1), tmp,'r')
        VLperM1_std(p1_ii) = std(TC_Conn_List{p1_ii,p2_ii,cell_type}.NumVL_of_M1); 
    end
end

