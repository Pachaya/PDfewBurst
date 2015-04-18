
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
plotSampleV = 1;
SPECIFIED_RATIO = 0; %TestSim_WT_InGauss0.2_W0.0005_50.00Hz_rEE0_rEI1_rIE2_rII1_Wmult10
Basal_Act_set  = [];

NoiseMEAN = 0;

SAVE_FIG = 1;

dirLoc = 'Model1_MinParam150403_fullsim/' ;
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
dirFig = 'Fig_testrTC/';
mkdir([dirLoc dirFig])
TSTOP = 4000;

rTC_LST = [120; 150; 170; 190; 210;];
wmTC_LST = [10 25 50 75 ];


LightAmp_LST = [0.4; ];
GPmVLw_mean_LST = [0.08];
GPmVLw_sig_LST =[ 0.01 ];

% Candidates Model 1 :
% I amp = 0.3, mean 0.06, sigma = 0.01, 0.005 ?0.03
% I amp = 0.5, mean 0.04, sigma = 0.03, 0.01 ?0.03
%PDfewBurst_GPmVLmd1_rTC120_wmTC10_WT_GPmInput_Amp0.3_Dur1000_GPmVLw_m0.06_sig0.01_InGauss0.2_IGmean-0.15_IGmeanSig0_W0.0015_SpecifiedPoisSpk_sig0.00Hz_T4000_trial3
NUM_TRIAL = 5;
ncells = 1150;
TRIAL_LST = 1 : NUM_TRIAL;

CUTTIME = 500;
PhotoInjT = 1500;
PhotoStop = 2500;
DelayT = 0;
BurstRange = 100;


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


for p1_ii = 1 : length(PARAM1)
    for p2_ii = 1 : length(PARAM2)
        for p3_ii = 1 : length(PARAM3)
            for p4_ii = 1 : length(PARAM4)
                for p5_ii = 1 : length(PARAM5)
                    
                    r_ii = p1_ii; wm_ii = p2_ii;
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
                        TSTOP = 4000;
                        GPmLightDur = 1000;
                        
                        rTC = rTC_LST(r_ii);
                        wmTC = wmTC_LST(wm_ii);
                        
                        GPmLight = LightAmp_LST(la_ii);
                        GPm_w_mn = GPmVLw_mean_LST(m_ii);
                        GPm_w_sg = GPmVLw_sig_LST(s_ii);
                        
                        
                        Simulation_Code = [coreFileName '_rTC' num2str(rTC) '_wmTC' num2str(wmTC) '_' cTxt '_' 'GPmInput_Amp' num2str(GPmLight) '_Dur' num2str(GPmLightDur) '_GPmVLw_m' num2str(GPm_w_mn) '_sig' num2str(GPm_w_sg) ...
                            '_InGauss' num2str(InGauss_STDEV) '_IGmean' num2str(NoiseMEAN) '_IGmeanSig' num2str(IGmeanSig) '_W' num2str(W_Weight) '_SpecifiedPoisSpk_sig0.00Hz_T' num2str(TSTOP) '_trial' num2str(TRIAL_NO)];
                        % PDfewBurst_GPmVLmd1_0del_rTC120_wmTC4_KO_GPmInput_Amp0.3_Dur1000_GPmVLw_m0.04_sig0_InGauss0.2_IGmean-0.15_IGmeanSig0_W0.0012_SpecifiedPoisSpk_sig0.00Hz_T4000_trial1
                        disp('==================================================================================================')
                        disp(Simulation_Code)
                        disp('######  Download M1 ')
                        %M1
                        Name_postfix = [ 'M1_' Simulation_Code];
                        getCenterPart =1;
                        %                         BaselineAct % Get all cell activity
                        
%                         tmp.E = E;      tmp.I = I;  tmp.All = All; %tmp.VL = VL; tmp.M1 = M1;
%                         if (cell_type == 1)
%                             M1.WT = tmp;
%                             cTxt = 'WT';
%                         elseif (cell_type == 2)
%                             M1.KO = tmp;
%                             cTxt = 'KO';
%                         end
%                         
                        
                
                        
                        %                 CheckFileExist( dirLoc, Name_postfix  )
                        disp('######  Download VL ')
                        Name_postfix = [ Simulation_Code];
                        %                         VL_PhotoactivationAll
                        
                        %                 if(0)
                        
%                         tmp.E = E;      tmp.I = I;  tmp.All = All; %tmp.VL = VL; tmp.M1 = M1;
                        if (cell_type == 1)
%                             VL.WT = tmp;
                            cTxt = 'WT'; Name_postfix_WT = Name_postfix;
                        elseif (cell_type == 2)
%                             VL.KO = tmp;
                            cTxt = 'KO'; Name_postfix_KO = Name_postfix;
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
                        end
                        
                        disp('Finish plot sample membrain voltage');
                        toc(tt)
                    end
                    
                    
                    
%                     Basal_Act.VL = VL;
%                     Basal_Act.M1 = M1;
                    
                    
                end
                
                %                     ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii } = Basal_Act;
                
            end
        end
    end
end
