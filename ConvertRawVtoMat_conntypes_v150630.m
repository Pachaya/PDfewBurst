%%
% New Simulation

% TRIAL_NO_LST_LST =  6:10; 
for TRIAL_NO_LST = 1 : 5

%% PARAMETERS
rTC_LST = [50:10:100];
wmTC_LST = [50:10:100];
% TRIAL_NO_LST = 5;
OSC_F_LST= 40;
OSC_Amp_LST = [ 0 0.5 1];


PARAM1 = rTC_LST;
lblTxt1 = 'Range of connection';
saveTxt1 = 'rTC';
titleTxt1 = 'Range_T_C';
PARAM2 = wmTC_LST;
lblTxt2 = 'Weighting Factor';
saveTxt2 = 'wmTC';
titleTxt2 = 'W_c_o_n_n';
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


% Connection Types
NumCnvrgntTypes = 3; CnvrgntTypes = cell(NumCnvrgntTypes,1);
% CnvrgntTypes{1}.FileName = 'GPmVLmd1_0del_KO2';
% CnvrgntTypes{2}.FileName = 'GPmVLmd1_0del_KO2_avgPWuniformTC';
% CnvrgntTypes{3}.FileName = 'GPmVLmd1_0del_KO2_avgPnegExpWTC' ;
CnvrgntTypes{1}.TitleName = 'G-G Connectivity: Gaussian, Strength: Gaussian';
CnvrgntTypes{2}.TitleName = 'S-U Connectivity: Uniform, Strength: Uniform';
CnvrgntTypes{3}.TitleName = 'S-E Connectivity: Uniform, Strength: Exponential';
CnvrgntTypes{1}.CodeName = 'GG';
CnvrgntTypes{2}.CodeName = 'SU';
CnvrgntTypes{3}.CodeName = 'SE';
CnvrgntTypes{1}.leg = 'P:Gauss, W:Gauss';
CnvrgntTypes{2}.leg = 'P:uniform, W:Uniform';
CnvrgntTypes{3}.leg = 'P:uniform, W:Exp';


PoisInputFr = 20;

PATH = SetPath;
dirLoc = [PATH 'TheoreticalSim/30-Jun-2015/']; % [PATH 'OscInput_Sim/']; % /OscInput_varyTCtype_Sim/
% dirFig = ['CompareConTypes_' num2str( PoisInputFr) 'Hz_' get_Parameters_RangeTxt( PARAMETERS,[1,2,4,5]) '/'];
% mkdir([dirLoc dirFig])
outLoc = ['SimResult\' date '/'];
mkdir(outLoc)
% CheckFile = zeros(3,   4 ,    3,     1,     1,     3);
%%   Files
%InputSpkTrain_GG_InputFR20_OscF40_OscrltAmp0.5_Wscale1e-05_W50_range150_Trial1_T5000
%SomaVolt_SU_InputFR20_OscF40_OscrltAmp0.5_Wscale1e-05_W75_range150_Trial1_T5000

%%

NEVER_LOAD_DATA = 1;
if (NEVER_LOAD_DATA)
    
    for ct_ii = 1 : NumCnvrgntTypes
        tc = tic();
       disp('-----------------------------------------------------------------------------------');
       disp(CnvrgntTypes{ct_ii}.TitleName)
       disp('-----------------------------------------------------------------------------------');
      ACT_Record = cell(ACT_Rec_size);
        for p1_ii = 1 : length(PARAM1)
            for p2_ii = 1 : length(PARAM2)
                for p3_ii = 1 : length(PARAM3)
                    for p4_ii = 1 : length(PARAM4)
                        for p5_ii = 1 : length(PARAM5)
                            t1 = tic();
                            r_ii = p1_ii; wm_ii = p2_ii;
                            TRIAL_NO = TRIAL_NO_LST(p3_ii);
                            of_ii = p4_ii; oa_ii = p5_ii;
                            
                            CnvrgentConnTypeTxt =  CnvrgntTypes{ct_ii}.CodeName;
                            Input_spk_avg_fr = PoisInputFr;
                            OSC_F = OSC_F_LST(of_ii);
                            OSC_rltAmp = OSC_Amp_LST(oa_ii);
                            W_SCALE = 0.00001;
                            weightingFactor=  wmTC_LST(wm_ii);
                            range = rTC_LST(r_ii);
                            
                            tstop = 5000;
                            CtrlCode = '';
                            
                            SimCode = sprintf( '%s_%s_InputFR%g_OscF%g_OscrltAmp%g_Wscale%g_W%g_range%g_Trial%g_T%g', CtrlCode, CnvrgentConnTypeTxt, Input_spk_avg_fr, OSC_F, OSC_rltAmp, W_SCALE, weightingFactor, range,TRIAL_NO,tstop) ;
                            disp(SimCode)
%                             CheckFile(ct_ii,p1_ii, p2_ii, p3_ii, p4_ii,  p5_ii) = CheckFileExist( dirLoc, SimCode);
                            
                            %Layer 1 
                            Read_Input_Spktrain
                            L1.spktrain  = SPKTRAIN;
                            
                            %Layer 2
                            Read_SomaVolt_to_spktrain
                            
                            L2.spktrain  = spktrain;
                            tmp.L1 = L1; tmp.L2 = L2;    tmp.SimulationCode = SimCode;                        
                            ACT_Record{p1_ii, p2_ii, p3_ii, p4_ii,  p5_ii} = tmp;
                            toc(t1);
                            
                        end
                    end
                end
            end
        end
        CnvrgntTypes{ct_ii}.ACT_Record = ACT_Record;
        toc(tc)
    end
end

save([outLoc 'CnvrgntTypes_' num2str( PoisInputFr) 'Hz_' get_Parameters_RangeTxt( PARAMETERS,[1,2,3]) '_' date '.mat'], 'CnvrgntTypes', 'PARAMETERS','-v7.3');
disp('Done & Saved')
disp([ '**************   TRIAL_NO_LST =  ' num2str(TRIAL_NO_LST)])
end