% Response Function
clc
close all
clear all


simCode  = 'ResFunc';
CoreName = '';

SAVE_FIG = 0 ; Close_Fig_aftr_save = 1;

Nsample = 100;
SynchLvl = 0;
Tcut = 500;
Tmeasure = 5000;
TSTOP = Tcut + Tmeasure;
tstop = TSTOP;
InputFR = 40;
W_SPK = 0.029;
W_VL_M1 = 0.002/Nsample ;


InputFR_LST = [100];
SynchLvl_LST = [0];%[0:1/Nsample:1];
W_VL_M1_LST = [0.0001 : 0.0002 : 0.0001*Nsample]; %[0.0001 : 0.0001 : 0.0001*Nsample];

WspkMult_LST = [1 1.5 2 2.5 3  5];
IGmean_LST = [0 -0.1 -0.5 -1 -1.5]; % 0
IGsig_LST = [0 0.2 0.5 1 ];

PARAM1 = InputFR_LST;
lblTxt1 = 'Input Frequency';
saveTxt1 = 'InputFR';
titleTxt1 = 'Input Frequency';
PARAM2 = IGmean_LST;
lblTxt2 = 'Mean of tonic Current';
saveTxt2 = 'IGmean';
titleTxt2 = 'IG_m_e_a_n';
PARAM3 = IGsig_LST;
lblTxt3 = 'Sigma of tonic current';
titleTxt3 = 'IGsig';
saveTxt3 = 'IG_s_i_g';
PARAM4 = WspkMult_LST;
lblTxt4 = 'Weighting factor';
titleTxt4 = 'WspkMult';
saveTxt4 = 'W_m_u_l_t';

N_Param = 4;
ACT_Rec_size = zeros(1,N_Param);
for ii = 1 : N_Param
    PARAMETERS{ii}.PARAM = eval(sprintf('PARAM%d',ii)) ;
    PARAMETERS{ii}.lblTxt = eval(sprintf('lblTxt%d',ii));
    PARAMETERS{ii}.titleTxt = eval(sprintf('titleTxt%d',ii));
    PARAMETERS{ii}.saveTxt = eval(sprintf('saveTxt%d',ii));
    ACT_Rec_size(ii) = eval(sprintf('length(PARAM%d)',ii));
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
% dirFig = ['Fig' get_Parameters_RangeTxt( PARAMETERS,[1,2,4,5]) '/'];
% mkdir([dirLoc dirFig])

for p1_ii = 1 : length(PARAM1)
  for p2_ii = 1 : length(PARAM2)
      for p3_ii = 1 : length(PARAM3)
          for p4_ii = 1 : length(PARAM4)
              plist = [p1_ii, p2_ii, p3_ii,p4_ii];
              
            inF_ii = p1_ii; sl_ii = 1; wVM_ii = 1;
            igm_ii = p2_ii; igs_ii = p3_ii; ww_ii = p4_ii;
            
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
            
            simCode =[ coreName 'Nsample' num2str(Nsample) '_TSTOP' num2str(TSTOP) '_InputFR' num2str(InputFR)  '_Wspk' num2str(Wspk) '_IGmean' num2str(IGmean) '_IGsig' num2str(IGsig)];
            disp(simCode)
%                 '_wSPK' num2str(W_SPK) '_wVLM1_' num2str(W_VL_M1)];
%             SomaVolt_WT_VL_Nsample100_TSTOP5500_InputFR15
            
            %WT_VL
            Name_postfix  = ['WT_VL_' simCode ];
            Check_Status(plist) = CheckFileExist( dirLoc, Name_postfix  );
               
          end
      end
  end
end




%%