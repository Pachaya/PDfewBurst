% Check spike train
clc 
clear all
close all


simCode  = 'ResFunc';
CoreName = '';

SAVE_FIG = 1 ; Close_Fig_aftr_save = 1;

Nsample = 100;
SynchLvl = 0;
Tcut = 500;
Tmeasure = 5000;
TSTOP = Tcut + Tmeasure;
tstop = TSTOP;
InputFR = 40;
W_SPK = 0.029;
W_VL_M1 = 0.002/Nsample ;


InputFR_LST = [ 10 25:25:450];
SynchLvl_LST = [0];%[0:1/Nsample:1];
W_VL_M1_LST = [0.0001 : 0.0002 : 0.0001*Nsample]; %[0.0001 : 0.0001 : 0.0001*Nsample];

WspkMult_LST = [1]; %[1 1.5 2 2.5 3  5];
IGmean_LST = 0; %[0 -0.1 -0.5 -1 -1.5]; % 0
IGsig_LST = [0];
OSC_F_LST =  [20, 40, 10];
OSC_AMP_LST = [0 0.5 1];

PARAM1 = InputFR_LST;
lblTxt1 = 'Input Frequency';
saveTxt1 = 'InputFR';
titleTxt1 = 'Input Frequency';
PARAM2 = OSC_F_LST;
lblTxt2 = 'Oscilation Frequencyt';
saveTxt2 = 'OscF';
titleTxt2 = 'OSC Freq';
PARAM3 = OSC_AMP_LST;
lblTxt3 = 'Oscilation Amplitude relative to mean input FR';
saveTxt3 = 'OscAmp';
titleTxt3 = 'OSC Amp';
PARAM4 = WspkMult_LST;
lblTxt4 = 'Weighting factor';
saveTxt4 = 'WspkMult';
titleTxt4 = 'W_m_u_l_t';

N_Param = 4;
ACT_Rec_size = zeros(1,N_Param);
for mm = 1 : N_Param
    PARAMETERS{mm}.PARAM = eval(sprintf('PARAM%d',mm)) ;
    PARAMETERS{mm}.lblTxt = eval(sprintf('lblTxt%d',mm));
    PARAMETERS{mm}.titleTxt = eval(sprintf('titleTxt%d',mm));
    PARAMETERS{mm}.saveTxt = eval(sprintf('saveTxt%d',mm));
    ACT_Rec_size(mm) = eval(sprintf('length(PARAM%d)',mm));
end
if(length(ACT_Rec_size) == 1)
    ACT_Rec_size = [ACT_Rec_size 1];
end
ACT_Record = cell(ACT_Rec_size);
Check_Status = zeros(ACT_Rec_size);
coreName = '';

% Directory
PATH = SetPath;
dirLoc = [PATH 'ResFunc_cellType/'];
dirFig = ['DelThisfold/']; %['smallResFunc_' get_Parameters_RangeTxt( PARAMETERS,[1:4]) '/'];
mkdir([dirLoc dirFig])

for p1_ii = 1 : length(PARAM1)
    for p2_ii = 1 : length(PARAM2)
        for p3_ii = 1 : length(PARAM3)
            for p4_ii = 1 : length(PARAM4)
                plist = [p1_ii, p2_ii, p3_ii,p4_ii];
                
                inF_ii = p1_ii; sl_ii = 1; wVM_ii = 1;
                igm_ii = 1; igs_ii = 1; ww_ii = p4_ii;
                of_ii = p2_ii; oa_ii = p3_ii;
                
                W_VL_M1= W_VL_M1_LST(wVM_ii);
                SynchLvl = SynchLvl_LST(sl_ii);
                InputFR = InputFR_LST(inF_ii);
                CUTTIME = 500;
                TSTOP = 5500;
                
                IGmean = IGmean_LST(igm_ii);
                IGsig = IGsig_LST(igs_ii);
                TestW = WspkMult_LST(ww_ii);
                Wscale = 0.001;
                if(IGmean == 0)
                    Wspk = Wscale*TestW;
                else
                    Wspk =  IGmean*-10*Wscale*TestW;
                end
                
                Osc_F = OSC_F_LST(of_ii);
                Osc_amp = OSC_AMP_LST(oa_ii);
                PHASE =0;
                ampTxt = sprintf('%2.2f',Osc_amp);
                %RecordSpkPvalues_VL_Nsample100_TSTOP5500_InputFR10_Wspk0.003_IGmean0_IGsig0_oscF20Hz_amp0.10_phase0
                simCode =[ coreName 'Nsample' num2str(Nsample) '_TSTOP' num2str(TSTOP) '_InputFR' num2str(InputFR)  ...
                    '_Wspk' num2str(Wspk) '_IGmean' num2str(IGmean) '_IGsig' num2str(IGsig) ...
                    '_oscF' num2str(Osc_F) 'Hz_amp' ampTxt '_phase' num2str(PHASE) ];
                disp(simCode)
                

            
          
            Name_postfix  = [ simCode ];
%             Check_Status(plist) = CheckFileExist( dirLoc, Name_postfix  );
            
%             RecordSpkTrain_VL_Nsample100_TSTOP5500_InputFR10_Wspk0.0015_IGmean0_IGsig0_oscF10Hz_amp0.00_phase0
            
              if exist([dirLoc 'RecordSpkTrain_VL_' Name_postfix '.txt'], 'file')
                        disp(' ----- Yes -----')                       
                        found =1;
                    else
                        disp('### Sorry ###')
                        found = 0;
              end
                    Check_Status(plist) = found; 
               
          end
      end
  end
end




%%