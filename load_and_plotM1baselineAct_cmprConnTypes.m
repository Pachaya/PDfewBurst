% Load and Get M1 base line activity only

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
plotSampleV = 0;


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

rTC_LST =  [ 50 100 150 200 250 ]; %
wmTC_LST = [ 50 100 200 300 400 500]; %[10 25 50 75 100];

LightAmp_LST = [0.5];
LightDur_LST = [1000];
GPmVLw_mean_LST = [ 0.5];
GPmVLw_sig_LST =[0];

OSC_F_LST = [20 40 ];
OSC_Amp_LST = [0 1/2 1];
OSC_phase_LST = 0;

TRIAL_NO_LST = 1;
Noise_MEAN_LST = 0; %[0 -0.5 -1 -1.5 -2 ];
InGauss_STDEV_LST = 0; %[ 0 0.2 0.5 1 ];
spkW_LST = [0.001];
PoisInputFr_LST = 50;

PARAM1 = rTC_LST;
lblTxt1 = 'Range of thalamocortical connection';
saveTxt1 = 'rTC';
titleTxt1 = 'Range_T_C';
PARAM2 = wmTC_LST;
lblTxt2 = 'Weight of thalamocortical connection';
saveTxt2 = 'wmTC';
titleTxt2 = 'W_T_C';
PARAM3 = TRIAL_NO_LST ;
lblTxt3 = 'Trial#';
saveTxt3 = 'trial';
titleTxt3 = 'Trial NO.';
PARAM4 = OSC_F_LST;
lblTxt4 = 'Oscillation Frequency';
saveTxt4 = 'OSC_F';
titleTxt4 = 'Osc Freq';
PARAM5 = OSC_Amp_LST;
lblTxt5 = 'Oscillation Amplitude relative to mean FR';
saveTxt5 = 'OSC_amp';
titleTxt5 = 'Osc Amp';
PARAM6 = spkW_LST;
lblTxt6 = 'Input spike weight';
saveTxt6 = 'Wspk';
titleTxt6 = 'W_s_p_k';


                    
N_Param = 6;
ACT_Rec_size = zeros(1,N_Param);
for ii = 1 : N_Param
    PARAMETERS{ii}.PARAM = eval(sprintf('PARAM%d',ii)) ;
    PARAMETERS{ii}.lblTxt = eval(sprintf('lblTxt%d',ii));
    PARAMETERS{ii}.titleTxt = eval(sprintf('titleTxt%d',ii));
    PARAMETERS{ii}.saveTxt = eval(sprintf('saveTxt%d',ii));
    ACT_Rec_size(ii) = eval(sprintf('length(PARAM%d)',ii));
end
ACT_Record = cell(ACT_Rec_size);

% Connection Types
NumCnvrgntTypes = 3; CnvrgntTypes = cell(NumCnvrgntTypes,1);
CnvrgntTypes{1}.FileName = 'GPmVLmd1_0del_KO2'; 
CnvrgntTypes{2}.FileName = 'GPmVLmd1_0del_KO2_avgPWuniformTC'; 
CnvrgntTypes{3}.FileName = 'GPmVLmd1_0del_KO2_avgPnegExpWTC' ; 
CnvrgntTypes{1}.TitleName = 'Connectivity: Gaussian, Strength: Gaussian'; 
CnvrgntTypes{2}.TitleName = 'Connectivity: Uniform, Strength: Uniform'; 
CnvrgntTypes{3}.TitleName = 'Connectivity: Uniform, Strength: Negative exponential'; 
CnvrgntTypes{1}.CodeName = 'gaussPgaussW'; 
CnvrgntTypes{2}.CodeName = 'avgPuniformW'; 
CnvrgntTypes{3}.CodeName = 'avgPnegecpW'; 
CnvrgntTypes{1}.leg = 'P:Gauss, W:Gauss'; 
CnvrgntTypes{2}.leg = 'P:uniform, W:uniform';
CnvrgntTypes{3}.leg = 'P:uniform, W:Neg Exp'; 

 PoisInputFr = 50;
 
% Directory

PATH = SetPath;
dirLoc = [PATH 'OscInput_varyTCtype_Sim/']; % [PATH 'OscInput_Sim/']; % /OscInput_varyTCtype_Sim/
dirFig = ['CompareConTypes_' num2str( PoisInputFr) 'Hz_' get_Parameters_RangeTxt( PARAMETERS,[1,2,4,5,6]) '/'];
mkdir([dirLoc dirFig])
NEVER_LOAD_DATA = 0;
if (NEVER_LOAD_DATA)
for ct_ii = 1 : NumCnvrgntTypes 
 ACT_Record = cell(ACT_Rec_size);
 if ct_ii == 1 
    dirLoc = [PATH 'OscInput_Sim/'];
    dirFig = ['../OscInput_varyTCtype_Sim/' 'CompareConTypes_' num2str( PoisInputFr) 'Hz_' get_Parameters_RangeTxt( PARAMETERS,[1,2,4,5,6]) '/'];    
 else
    dirLoc = [PATH 'OscInput_varyTCtype_Sim/']; % [PATH 'OscInput_Sim/']; % /OscInput_varyTCtype_Sim/
    dirFig = ['CompareConTypes_' num2str( PoisInputFr) 'Hz_' get_Parameters_RangeTxt( PARAMETERS,[1,2,4,5,6]) '/'];
 end
for p1_ii = 1 : length(PARAM1)
    for p2_ii = 1 : length(PARAM2)
        for p3_ii = 1 : length(PARAM3)
            for p4_ii = 1 : length(PARAM4)
                for p5_ii = 1 : length(PARAM5)
                    for p6_ii = 1 : length(PARAM6)
                        
                     r_ii = p1_ii; wm_ii = p2_ii;
                    la_ii = 1;  m_ii = 1; s_ii = 1; ld_ii= 1;   
                    of_ii = p4_ii; oa_ii = p5_ii; op_ii = 1;
                    TRIAL_NO = p3_ii; w_ii = p6_ii;
                    gn_ii = 1; g_ii = 1;
                        

                        for cell_type = 1 : 2
                                      if (cell_type == 1)
                            cTxt = 'WT';
                        elseif (cell_type == 2)
                            cTxt = 'KO';
                        end
                        
                        % PDfewBurst_GPmVLmd1_rTC120_wmTC10_WT_GPmInput_Amp0.3_Dur1000_GPmVLw_m0.06_sig0.01_InGauss0.2_IGmean-0.15_IGmeanSig0_W0.0015_SpecifiedPoisSpk_sig0.00Hz_T4000_trial3
                        
                        coreFileName =  CnvrgntTypes{ct_ii}.FileName; %'GPmVLmd1_0del_KO2_avgPWuniformTC' ; %% 'GPmVLmd1_0del_KO2'  for Gaussian distribution ,   'GPmVLmd1_0del_KO2_uniformP63WTC' for uniform distribution , 'GPmVLmd1_0del_KO2_negExpWTC' for constant connectivity with random w from negative exponential distribution
                        
                        InGauss_STDEV = InGauss_STDEV_LST(gn_ii); %0.2;, 0.3
                        NoiseMEAN = Noise_MEAN_LST(g_ii);
                        IGmeanSig = 0;
                        W_Weight = spkW_LST(w_ii);
                        PoisInputFr = 50;
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
                        
%                        GPmVLmd1_0del_KO2_rTC250_wmTC50_WT_IGmean-2_IGmeanSig0_W0.029_20.00Hz_oscF40Hz_amp0.50_phase0_T3000_trial1                   
%                        GPmVLmd1_0del_KO2_rTC100_wmTC40_WT_InGauss1_IGmean-2_IGmeanSig0_W0.029_10.00Hz_oscF40Hz_amp0.50_phase0_T3000_trial1 
                        Name_postfix = [ Simulation_Code];
                            
                            disp('######  Download M1 ')
                            %M1
                            Name_postfix = [ 'M1_' Simulation_Code];
                            getCenterPart = 1; % Only at center to avoide the border problem
                            PhotoactivationAll
                            
                            tmp.E = E;      tmp.I = I;  tmp.All = All; %tmp.VL = VL; tmp.M1 = M1;
                            tmp.Simulation_Code = Name_postfix; % for M1, simulation code has initial 'M1_'
                            if (cell_type == 1)
                                M1.WT = tmp;
                                cTxt = 'WT';
                            elseif (cell_type == 2)
                                M1.KO = tmp;
                                cTxt = 'KO';
                            end
                            clear tmp
                            %                 CheckFileExist( dirLoc, Name_postfix  )
                            disp('######  Download VL ')
                            Name_postfix = [ Simulation_Code];
                            getCenterPart = 0;
                            PhotoactivationAll
                            %                 if(0)
                            tmp.E = E;      tmp.I = I;  tmp.All = All; %tmp.VL = VL; tmp.M1 = M1;
                            tmp.Simulation_Code = Name_postfix;
                            if (cell_type == 1)
                                VL.WT = tmp;
                                cTxt = 'WT'; Name_postfix_WT = Name_postfix;
                            elseif (cell_type == 2)
                                VL.KO = tmp;
                                cTxt = 'KO'; Name_postfix_KO = Name_postfix;
                            end
                            clear tmp
                            
                        end
                        
                        if (plotSampleV)
                            %soma's volt
                            tt= tic();
                            
                            samVL = 385;
                            samM1 = 87;
                            samTstr = PhotoInjT-100; samTstp = PhotoStop + DelayT+BurstRange+100;
                            fgSmpl = figure; set(fgSmpl,'position',[  377          94        1310         894]);  set(gcf,'PaperPositionMode','auto')
                            subplot(221);
                            [NEcell,NIcell, Tstop,somaVallWT] = get_save_vec_file_BigNet([dirLoc 'SomaVolt_M1_' Name_postfix_WT '.txt'] );
                            plot(samTstr:samTstp, somaVallWT(samM1,samTstr:samTstp), 'r'); hold on;
                            plot(samTstr:samTstp, somaVallWT(samM1+1,samTstr:samTstp), 'm');
                            ylim([-90 50]); xlim([samTstr samTstp]);
                            title(['WT M1']);
                            subplot(222);
                            [NEcell,NIcell, Tstop,somaVallWT] = get_save_vec_file_BigNet([dirLoc 'SomaVolt_M1_' Name_postfix_KO '.txt'] );
                            plot(samTstr:samTstp, somaVallWT(samM1,samTstr:samTstp), 'r'); hold on;
                            plot(samTstr:samTstp, somaVallWT(samM1+1,samTstr:samTstp), 'm');
                            ylim([-90 50]); xlim([samTstr samTstp]);
                            title(['KO M1']);
                            subplot(223);
                            [NEcell,NIcell, Tstop,somaVallWT] = get_save_vec_file_BigNet([dirLoc 'SomaVolt_' Name_postfix_WT '.txt'] );
                            plot(samTstr:samTstp, somaVallWT(samVL,samTstr:samTstp), 'r'); hold on;
                            plot(samTstr:samTstp, somaVallWT(samVL+1,samTstr:samTstp), 'm');
                            ylim([-90 50]); xlim([samTstr samTstp]);
                            title(['WT VL']);
                            subplot(224);
                            [NEcell,NIcell, Tstop,somaVallWT] = get_save_vec_file_BigNet([dirLoc 'SomaVolt_' Name_postfix_KO '.txt'] );
                            plot(samTstr:samTstp, somaVallWT(samVL,samTstr:samTstp), 'r'); hold on;
                            plot(samTstr:samTstp, somaVallWT(samVL+1,samTstr:samTstp), 'm');
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
                                if(Close_Fig_aftr_save)
                                    close(fgSmpl);
                                end
                            end
                            
                            disp('Finish plot sample membrain voltage');
                            toc(tt)
                        end
                        
                        
                        
                        Basal_Act.VL = VL;
                        Basal_Act.M1 = M1;
                        
                        % T-Test goes here // Remove in load_and_$$$$ series
                        % for faster runtime
                        
                        ACT_Record{p1_ii, p2_ii, p3_ii, p4_ii,  p5_ii, p6_ii} = Basal_Act;
                        
                        clear Basal_Act
                    end
                    
                    
                end
            end
        end
    end
end
%% Save Data

    contype =   CnvrgntTypes{ct_ii}.CodeName;
    CnvrgntTypes{ct_ii}.ACT_Record = ACT_Record;
    save([dirLoc dirFig 'Activity_' contype '_' date '.mat' ], 'ACT_Record','PARAMETERS','-v7.3');
    CnvrgntTypes{ct_ii}.matFile = [dirLoc dirFig 'Activity_' contype '_' date '.mat' ];
    disp('============================================================================================');
    disp('============================================================================================');
    disp(CnvrgntTypes{ct_ii}.TitleName)
    disp('============================================================================================');
    disp('============================================================================================');
    clear ACT_Record
end
else
    
    
    matFile = {'Activity_gaussPgaussW_27-May-2015','Activity_avgPuniformW_27-May-2015', 'Activity_avgPnegecpW_27-May-2015'};
    for ct_ii = 1 : NumCnvrgntTypes
        CnvrgntTypes{ct_ii}.matFile = matFile{ct_ii};
        load([dirLoc dirFig matFile{ct_ii} '.mat' ]);
        CnvrgntTypes{ct_ii}.ACT_Record = ACT_Record;
        disp('============================================================================================');
        disp(CnvrgntTypes{ct_ii}.TitleName)
        disp('============================================================================================');
        
    end
end
%%
 save([dirLoc dirFig 'ActivityResult_allCnvrgntTypes_' date '.mat' ], 'CnvrgntTypes','PARAMETERS','-v7.3');


%% 
close all
for ct_ii = 1 : NumCnvrgntTypes 
 ACT_Record = CnvrgntTypes{ct_ii}.ACT_Record;
 if ct_ii == 1 
    dirLoc = [PATH 'OscInput_Sim/'];
    dirFig = ['../OscInput_varyTCtype_Sim/' 'CompareConTypes_' num2str( PoisInputFr) 'Hz_' get_Parameters_RangeTxt( PARAMETERS,[1,2,4,5,6]) '/'];    
 else
    dirLoc = [PATH 'OscInput_varyTCtype_Sim/']; % [PATH 'OscInput_Sim/']; % /OscInput_varyTCtype_Sim/
    dirFig = ['CompareConTypes_' num2str( PoisInputFr) 'Hz_' get_Parameters_RangeTxt( PARAMETERS,[1,2,4,5,6]) '/'];
 end
%% Check Osc F and Amp of one current case
% Check_OscF_and_Amp_of_one_case
%%  Get info for TC convergenc --> number of VL per M1 , average maxW , average summation of weight
Get_info_TC_convergence_numVLperM1

%%
CtypeTxt =  CnvrgntTypes{ct_ii}.TitleName;
codeTxt =   CnvrgntTypes{ct_ii}.CodeName;
% M1_BaselineActivity_Matrix_noIndvPlot  % ---> plot raw M1 activity 
M1_BaselineActivity_Matrix_noIndvPlot_normVL % ---> plot M1 activity  - VL activity (normalized with VL activity)
Diff_M1act_OSC_F_AMP_oneType_normVL % Plot different in M1 raw activity when osc F increase or when osc amp increase

end

%% Convergent rate ( = average of total weight per cell) vs  output activity
Cnvrgnt_rate_vs_normOutput

%%

% Normalize M1 activity with actual number of spikes in VL
% Gaussian vs Uniform ->  same volume of connectivity and strength
% May be explain why the brain use this rules of connection 
% Why can we use simple statistical rules to model
% Quatitatively analysis for each parameter cases ->What is the  biggest different in response when osc 40 ? osc 20 or when amp 1 ? amp 0
% -> input
% t frequency change  thn what happen
% Ex Pick one set ex r = 250 , how response different in each cases
% Pick one set of ex W, how response different in each case 

 %ACT_Record1 = CnvrgntTypes{1}.ACT_Record; 
 
% PARAM1 = rTC_LST;
% lblTxt1 = 'Range of thalamocortical connection';
% saveTxt1 = 'rTC';
% titleTxt1 = 'Range_T_C';
% PARAM2 = wmTC_LST;
% lblTxt2 = 'Weight of thalamocortical connection';
% saveTxt2 = 'wmTC';
% titleTxt2 = 'W_T_C';
% PARAM3 = TRIAL_NO_LST ;
% lblTxt3 = 'Trial#';
% saveTxt3 = 'trial';
% titleTxt3 = 'Trial NO.';
% PARAM4 = OSC_F_LST;
% lblTxt4 = 'Oscillation Frequency';
% saveTxt4 = 'OSC_F';
% titleTxt4 = 'Osc Freq';
% PARAM5 = OSC_Amp_LST;
% lblTxt5 = 'Oscillation Amplitude relative to mean FR';
% saveTxt5 = 'OSC_amp';
% titleTxt5 = 'Osc Amp';
% PARAM6 = spkW_LST;
% lblTxt6 = 'Input spike weight';
% saveTxt6 = 'Wspk';
% titleTxt6 = 'W_s_p_k';

%% Quatitatively analysis for each parameter cases ->What is the  biggest different in response when osc 40 ? osc 20 or when amp 1 ? amp 0
% For First Case 
% Different in 

%%
if(SAVE_WORKSPACE)
    save([dirLoc dirFig 'Saved_Workspace_result_' SPECIFIED_BASAL_ACT_CODENAME date '.mat' ],'-v7.3');
end