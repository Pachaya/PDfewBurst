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

rTC_LST =  [ 100 150 200 250 ]; %
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



 PoisInputFr = 50;
 
% Directory

PATH = SetPath;
dirLoc = [PATH 'OscInput_varyTCtype_Sim/']; % [PATH 'OscInput_Sim/']; % /OscInput_varyTCtype_Sim/
dirFig = ['uniformP63WTC_' num2str( PoisInputFr) 'Hz_' get_Parameters_RangeTxt( PARAMETERS,[1,2,4,5,6]) '/'];
mkdir([dirLoc dirFig])

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
                        
                        coreFileName = 'GPmVLmd1_0del_KO2_avgPWuniformTC' ; %% 'GPmVLmd1_0del_KO2'  for Gaussian distribution ,   'GPmVLmd1_0del_KO2_uniformP63WTC' for uniform distribution , 'GPmVLmd1_0del_KO2_negExpWTC' for constant connectivity with random w from negative exponential distribution
                        
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
%%
if(SAVE_BASAL_ACT)
    save([dirLoc dirFig 'Saved_Activity_result_' SPECIFIED_BASAL_ACT_CODENAME date '.mat' ], 'ACT_Record','PARAMETERS','-v7.3');
end



%% 


%% Check Osc F and Amp of one current case
p1_ii = 1; p2_ii =1; p3_ii = 1; p6_ii = 1;
nR = length(PARAM4); nC = length(PARAM5);
fosc = figure;  set(fosc, 'position', [   147         799        2094         471]); set(gcf,'PaperPositionMode','auto');
cnt1 =0;
fg_fft = figure;  set(fg_fft,'position',[ 458         600        1421         568]); set(gcf,'PaperPositionMode','auto');
cntcnt = 0;

for p4_ii = 1 : length(PARAM4)  % OSC F
    for p5_ii = 1 : length(PARAM5) % OSC Amp
        
        tmpTxt = get_Parameters_titleText(PARAMETERS,[4, 5],[p4_ii,p5_ii]);
        %%% for VL
        tmpWTtrain_VL =  ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii,p6_ii}.VL.WT.All.spktrain;
        tmpKOtrain_VL =  ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii,p6_ii}.VL.KO.All.spktrain;
        
        % Check Oscillation Frequency of VL
        
        T = TSTOP;
        tt = 0:1:T;
        
        % Spike Train
        [KOci, KOti] = find(tmpKOtrain_VL);
        [WTci, WTti] = find(tmpWTtrain_VL);
        
        % mean firing rate
        KOfmu = mean(tmpKOtrain_VL, 1)*1000;
        WTfmu = mean(tmpWTtrain_VL, 1)*1000;
        %         fosc = figure;  set(fosc, 'position', [496         763        1042         420]); set(gcf,'PaperPositionMode','auto');
        %         subplot(211);
        %         plot(WTfmu,'k');   hold on;  plot(KOfmu,'r'); legend('WT','KO'); title('Raw instantaneous average firing rate')
        
        gsig = 3;   % filter [ms]
        gfil = exp( -[-round(4*gsig):round(4*gsig)].^2 /2/gsig^2);
        gfil = gfil/sum(gfil(:));
        
        KOfmus = conv(KOfmu, gfil, 'same');
        WTfmus = conv(WTfmu, gfil, 'same');
        %         subplot(212); plot(WTfmus,'k');   hold on;  plot(KOfmus,'r'); legend('WT','KO'); title('Smooth instantaneous average firing rate');
        %         suptitle(tmpTxt)
        cnt1 = cnt1+1;
        figure(fosc); subplot(nR,nC,cnt1);
        plot(WTfmus,'k');   hold on;  plot(KOfmus,'r'); legend('WT','KO','location','best'); title(tmpTxt);  %title('Smooth instantaneous average firing rate');
        
        % FFT
        KOfft = abs(fftshift(fft(KOfmus)));
        WTfft = abs(fftshift(fft(WTfmus)));
        
        
        dHz = 1/T * 1000;
        Hz_ind = [-round(T/2):round(T/2)]*dHz;
        %         Hz_ind = Hz_ind(1:T);
        %         fg_fft = figure;  set(fg_fft,'position',[ 458         600        1421         568]); set(gcf,'PaperPositionMode','auto');
        %         subplot(121);
        figure(fg_fft);  cntcnt = cntcnt +1;
        subplot(nR,nC,cntcnt);
        plot(Hz_ind, WTfft,'k'); hold on;  plot(Hz_ind, KOfft,'r');  legend('WT','KO');
        xlim([1 100]); ylabel('FFT amplitude (a.u.)'); xlabel('Frequency (Hz)'); title(tmpTxt); %title('raw FFT')
        
        [Hind, KO_FFT, WT_FFT] = analyze_fft(Hz_ind, WTfft, KOfft);
        %         subplot(122);
        %         plot(Hind, WT_FFT.FF,'k') ;     hold on; plot(Hind, KO_FFT.FF,'r');  legend('WT','KO');
        %         xlim([1 100]); ylabel('FFT amplitude (a.u.)'); xlabel('Frequency (Hz)'); title('Smooth FFT');
        %         suptitle(tmpTxt);
        
    end
end
tmpTxt = get_Parameters_titleText(PARAMETERS, [1:3, 6], [p1_ii p2_ii p3_ii p6_ii ]);
figure(fosc); suptitle({['Smooth instantaneous average firing rate (Gaussian filter size = ' num2str(gsig) ')'], tmpTxt});
figure(fg_fft); suptitle({'Raw FFT', tmpTxt});
if (SAVE_FIG)
    tmpTxt = get_Parameters_saveText(PARAMETERS, [1:3 6], [p1_ii p2_ii p3_ii p6_ii]);
    ffig = [ dirLoc dirFig 'ChecksampleVLInst_Avg_FR_' tmpTxt ];
    saveas(  fosc, [ffig '.jpg'], 'jpg');        saveas(  fosc, [ffig '.fig'], 'fig');
    ffig = [ dirLoc dirFig 'ChecksampleVL_Osc_fft_' tmpTxt ];
    saveas(  fg_fft, [ffig '.jpg'], 'jpg');        saveas(  fg_fft, [ffig '.fig'], 'fig');
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get info for TC convergenc --> number of VL per M1 , average maxW , average summation of weight
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p3_ii = 1; p4_ii = 1; p5_ii =1;
ssize = 1500; dist = ssize/2; pmax =0.85; W_scale = 1E-05;
fprintf('RangeTC\t\t wmTC\t\t numVL/M1\t\t sumW\t\t maxW\t\n');
numVL_M1_LST = zeros(length(PARAM1),length(PARAM2));
sumW_LST = 		 zeros(length(PARAM1),length(PARAM2));
maxW_LST  = zeros(length(PARAM1),length(PARAM2));

for p1_ii = 1 : length(PARAM1)
    for p2_ii = 1 : length(PARAM2)
        %         disp('##################################################################################################')
        %        disp( [lblTxt1 ' = ' num2str(PARAM1(p1_ii)) ]);
        %          disp( [lblTxt2 ' = ' num2str(PARAM2(p2_ii)) ]);
        rangeTC = PARAM1(p1_ii);      sigTC = rangeTC/sqrt(2);
        wmTC = PARAM2(p2_ii);
        
        Name_postfix_WT  = ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii }.VL.WT.Simulation_Code;  % The connection is same for WT and KO
        %    Name_postfix_KO  = ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii }.VL.KO.Simulation_Code;
        display = 0;
        [TC_basedOnM1_WT, TC_sumW_WT, TC_maxW_WT, TC_numVL_WT ]               = ExtractTC_info( dirLoc, Name_postfix_WT,  display);
        %disp([ num2str(rangeTC) '        ' num2str(mean(TC_numVL_WT)) '        ' num2str(specified_wmTC) '        ' num2str(mean(TC_sumW_WT)) '        ' num2str(mean(TC_maxW_WT))  ])
        fprintf(' %3.0f\t\t%7.4f\t\t%7.4X\t\t%7.4X\t\t%7.4X\n', rangeTC, wmTC, mean(TC_numVL_WT), mean(TC_sumW_WT), mean(TC_maxW_WT));
        numVL_M1_LST(p1_ii,p2_ii) = mean(TC_numVL_WT);
        sumW_LST(p1_ii,p2_ii) = 	mean(TC_sumW_WT);
        maxW_LST(p1_ii,p2_ii)  = mean(TC_maxW_WT);
    end
end
VLperM1 = numVL_M1_LST(:,1);
%%
M1_BaselineActivity_Matrix



%%
if(SAVE_WORKSPACE)
    save([dirLoc dirFig 'Saved_Workspace_result_' SPECIFIED_BASAL_ACT_CODENAME date '.mat' ],'-v7.3');
end