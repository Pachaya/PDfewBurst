
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
TSTOP = 4000;

CUTTIME = 500;
PhotoInjT = 1500;
PhotoStop = 2500;
DelayT = 0;
BurstRange = 100;

% Parameters setting
NUM_TRIAL = 5;
ncells = 1150;
M1_ncells = 166;
TRIAL_LST = 1 : NUM_TRIAL;

rTC_LST = [25 50 75 100 150 200 250 ]; %
wmTC_LST = [20 30 40 50 ]; %[10 25 50 75 100];

LightAmp_LST = [0.5];
LightDur_LST = [1000];
GPmVLw_mean_LST = [ 0.5];
GPmVLw_sig_LST =[0];

OSC_F_LST = [20 40]; 
OSC_phase_LST = 0;


PARAM1 = rTC_LST;
lblTxt1 = 'Range of thalamocortical connection';
saveTxt1 = 'rTC';
titleTxt1 = 'Range_T_C';
PARAM2 = wmTC_LST;
lblTxt2 = 'Weight of thalamocortical connection';
saveTxt2 = 'wmTC';
titleTxt2 = 'W_T_C';
PARAM3 = GPmVLw_mean_LST;
lblTxt3 = 'GPm-VL weight mean';
saveTxt3 = 'GPmVLw_m';
titleTxt3 = 'W_m_e_a_n';
PARAM4 = OSC_F_LST; % TRIAL_LST;
lblTxt4 = 'Frequency of oscillating input (Hz)'; 
saveTxt4 = 'OscF'; %'trial';
titleTxt4 = 'Osc Freq';
PARAM5 = OSC_phase_LST; % TRIAL_LST;
lblTxt5 = 'Phase of oscillating input (Hz)';
saveTxt5 = 'OscPhase'; %'trial';
titleTxt5 = 'Osc Phase';

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
% Check_Status = zeros(ACT_Rec_size);

% Directory

dirLoc = 'Test_Osc_Input/' ;
[s,computername] = dos('ECHO %COMPUTERNAME%');
switch computername
    case 'USER-PC'
        dirLoc = [ 'D:\Pachaya\150313 Code\SimResult\' dirLoc];
    case 'VSLAB-PC'
        dirLoc = ['E:\PDmodelFewBurst\SimResult\' dirLoc];
end
dirFig = ['Fig_OSC_big' get_Parameters_RangeTxt( PARAMETERS,[1,2]) '/'];
mkdir([dirLoc dirFig])




%Start donwnload simulation data
for p1_ii = 1 : length(PARAM1)
    for p2_ii = 1 : length(PARAM2)
        for p3_ii = 1 : length(PARAM3)
            for p4_ii = 1 : length(PARAM4)
                for p5_ii = 1 : length(PARAM5)
                    
                    r_ii = p1_ii; wm_ii = p2_ii;
                    la_ii = 1;  m_ii = p3_ii; s_ii = 1; ld_ii= 1; 
                    of_ii = p4_ii; op_ii = p5_ii;
                    TRIAL_NO = 1;
                    
                    for cell_type = 1 : 2
                        tmp = [];
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
                        
                        rTC = rTC_LST(r_ii);
                        wmTC = wmTC_LST(wm_ii);
                        
                        GPmLight = LightAmp_LST(la_ii);
                        GPm_w_mn = GPmVLw_mean_LST(m_ii);
                        GPm_w_sg = GPmVLw_sig_LST(s_ii);
                        
                        osc_f = OSC_F_LST(of_ii);
                        phase = OSC_phase_LST(op_ii);
                        
                        CUTTIME = 500;
                        PhotoInjT = 1500;
                        PhotoStop = PhotoInjT + GPmLightDur;
                        DelayT = 0;
                        BurstRange = 100;
                        
                        
                        %GPmVLmd1_0del_rTC150_wmTC50_KO_GPmInput_Amp0.5_Dur1000_GPmVLw_m0.5_sig0_InGauss0.2_IGmean0_IGmeanSig0_W0.029_SpecifiedPoisSpk_sig0.00Hz_oscF40Hz_phase0_T4000_trial1
                        
                        Simulation_Code = [coreFileName '_rTC' num2str(rTC) '_wmTC' num2str(wmTC) '_' cTxt '_' 'GPmInput_Amp' num2str(GPmLight) '_Dur' num2str(GPmLightDur) '_GPmVLw_m' num2str(GPm_w_mn) '_sig' num2str(GPm_w_sg) ...
                            '_InGauss' num2str(InGauss_STDEV) '_IGmean' num2str(NoiseMEAN) '_IGmeanSig' num2str(IGmeanSig) '_W' num2str(W_Weight) '_SpecifiedPoisSpk_sig0.00Hz_oscF' num2str(osc_f) 'Hz_phase' num2str(phase) '_T' num2str(TSTOP) '_trial' num2str(TRIAL_NO)];
                        % PDfewBurst_GPmVLmd1_0del_rTC120_wmTC4_KO_GPmInput_Amp0.3_Dur1000_GPmVLw_m0.04_sig0_InGauss0.2_IGmean-0.15_IGmeanSig0_W0.0012_SpecifiedPoisSpk_sig0.00Hz_T4000_trial1
                        disp('==================================================================================================')
                        disp(Simulation_Code)
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
                    
                    
                    
                    Basal_Act.VL = VL;
                    Basal_Act.M1 = M1;
                    
                    % T-Test goes here // Remove in load_and_$$$$ series
                    % for faster runtime
                    ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii } = Basal_Act;
                    
                    clear Basal_Act
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
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    GPM - VL Connnection
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Analyse_GPmVLConnection

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the instantaneous average firing rate of all cases
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SIGMA_fr_All = 30;
Close_Fig_aftr_save = 1;
for p1_ii = 1 : length(PARAM1)
    for p2_ii = 1 : length(PARAM2)
        for p3_ii = 1 : length(PARAM3)
            for p4_ii = 1 : length(PARAM4)
                for p5_ii = 1 : length(PARAM5)
                    
                    %%% for VL
                    tmpWTtrain_VL =  ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii}.VL.WT.All.spktrain;
                    tmpKOtrain_VL =  ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii}.VL.KO.All.spktrain;
                    %Compare WT and KO
                    avgFR_WT_VL = getInstantaneousFiringRate(tmpWTtrain_VL, SIGMA_fr_All, RES);
                    avgFR_KO_VL = getInstantaneousFiringRate(tmpKOtrain_VL, SIGMA_fr_All, RES);
                    
                    %%% for M1
                    tmpWTtrain_M1 =  ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii}.M1.WT.All.spktrain;
                    tmpKOtrain_M1 =  ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii}.M1.KO.All.spktrain;
                    %Compare WT and KO
                    avgFR_WT_M1 = getInstantaneousFiringRate(tmpWTtrain_M1, SIGMA_fr_All, RES);
                    avgFR_KO_M1 = getInstantaneousFiringRate(tmpKOtrain_M1, SIGMA_fr_All, RES);
                    
                    %Plot and save
                    fFR = figure; set(gcf, 'position', [ 624   256   733   722]); set(gcf,'PaperPositionMode','auto');
                    subplot(211); plot(avgFR_WT_M1,'k');  hold on;   plot(avgFR_KO_M1,'m');  hold on;
                    title('M1');   ylabel('Average firing rate (Hz)');            xlabel('Time(ms)'); legend('WT','KO','location','best');
                    
                    subplot(212); plot(avgFR_WT_VL,'k');  hold on;   plot(avgFR_KO_VL,'m');  hold on;
                    title('VL');   ylabel('Average firing rate (Hz)');            xlabel('Time(ms)'); legend('WT','KO','location','best');
                    
                    tmpTxt = get_Parameters_titleText(PARAMETERS, [1:5], [p1_ii p2_ii p3_ii p4_ii p5_ii]);
                    suptitle({'Average firing rate Fr(t)', tmpTxt});
                    
                    
                    % Oscillation 
                    

 
    T = TSTOP;
    tt = 0:1:T;
    
        % Spike Train
        [KOci, KOti] = find(tmpKOtrain_VL);
        [WTci, WTti] = find(tmpWTtrain_VL);
        
        % mean firing rate
        KOfmu = mean(tmpKOtrain_VL, 1)*1000;
        WTfmu = mean(tmpWTtrain_VL, 1)*1000;
        fosc = figure;  set(fosc, 'position', [496         763        1042         420]); set(gcf,'PaperPositionMode','auto');
        subplot(211);
        plot(WTfmu,'k');   hold on;  plot(KOfmu,'r'); legend('WT','KO'); title('Raw instantaneous average firing rate')
         
        gsig = 3;   % filter [ms]
        gfil = exp( -[-round(4*gsig):round(4*gsig)].^2 /2/gsig^2);
        gfil = gfil/sum(gfil(:));
        
        KOfmus = conv(KOfmu, gfil, 'same');
        WTfmus = conv(WTfmu, gfil, 'same');
        subplot(212); plot(WTfmus,'k');   hold on;  plot(KOfmus,'r'); legend('WT','KO'); title('Smooth instantaneous average firing rate');
        suptitle(tmpTxt)
        
        % FFT
        KOfft = abs(fftshift(fft(KOfmus)));
        WTfft = abs(fftshift(fft(WTfmus)));

       
        dHz = 1/T * 1000;
        Hz_ind = [-round(T/2):round(T/2)]*dHz;
%         Hz_ind = Hz_ind(1:T);
         fg_fft = figure;  set(fg_fft,'position',[ 458         600        1421         568]); set(gcf,'PaperPositionMode','auto');
         subplot(121);
         plot(Hz_ind, WTfft,'k'); hold on;  plot(Hz_ind, KOfft,'r');  legend('WT','KO');
         xlim([1 100]); ylabel('FFT amplitude (a.u.)'); xlabel('Frequency (Hz)'); title('raw FFT')

        [Hind, KO_FFT, WT_FFT] = analyze_fft(Hz_ind, WTfft, KOfft);
        subplot(122);
        plot(Hind, WT_FFT.FF,'k') ;     hold on; plot(Hind, KO_FFT.FF,'r');  legend('WT','KO');
        xlim([1 100]); ylabel('FFT amplitude (a.u.)'); xlabel('Frequency (Hz)'); title('Smooth FFT');
        suptitle(tmpTxt);
        
                    
                    
                    
                    
                    
                    if (SAVE_FIG)
                        tmpTxt = get_Parameters_saveText(PARAMETERS, [1:5], [p1_ii p2_ii p3_ii p4_ii p5_ii]);
                        ffig = [ dirLoc dirFig 'Inst_Avg_FR_' tmpTxt ];
                        saveas(  fFR, [ffig '.jpg'], 'jpg');        saveas(  fFR, [ffig '.fig'], 'fig');
                        ffig = [ dirLoc dirFig 'Inst_Avg_FR_VL_' tmpTxt ];
                        saveas(  fosc, [ffig '.jpg'], 'jpg');        saveas(  fosc, [ffig '.fig'], 'fig');
                        ffig = [ dirLoc dirFig 'VL_Osc_fft_' tmpTxt ];
                        saveas(  fg_fft, [ffig '.jpg'], 'jpg');        saveas(  fg_fft, [ffig '.fig'], 'fig');
                        if (Close_Fig_aftr_save)
                            close(fFR);  close(fosc); close(fg_fft); 
                        end
                    end
                    clear  tmpWTtrain_VL  tmpKOtrain_VL tmpWTtrain_M1  tmpKOtrain_M1;
                end
            end
        end
    end
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


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   M1  Activity
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for p3_ii = 1 : length(PARAM3)
    for p4_ii = 1 : length(PARAM4)
        for p5_ii = 1 : length(PARAM5)
            
            CntBurstSpkWT = zeros(length(PARAM1),length(PARAM2));
            CntBurstSpkKO = zeros(length(PARAM1),length(PARAM2));
            M1_burstSpkWT = zeros(length(PARAM1),length(PARAM2));
            M1_burstSpkKO = zeros(length(PARAM1),length(PARAM2));
            
            %
            SAVE_FIG=1; Close_Fig_aftr_save =0;
            
            
            fRB = figure; set(fRB,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
            fBspk = figure; set(fBspk,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
            cnt = 0;
            nR = length(PARAM1); nC = length(PARAM2);
            for p1_ii = 1 : length(PARAM1)
                for p2_ii = 1 : length(PARAM2)
                    
                    BasalAct = ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii}.M1; %% for M1
                    
                    M1_burstSpkWT(p1_ii,p2_ii) = mean(BasalAct.WT.All.BurstSpk)./BurstRange*1000;
                    M1_burstSpkKO(p1_ii,p2_ii) = mean(BasalAct.KO.All.BurstSpk)./BurstRange*1000;
                    
                    simTxt =      get_Parameters_titleText(PARAMETERS, [1:5], [ p1_ii,p2_ii, p3_ii p4_ii p5_ii]);
                    %['VL ' titleTxt3 ' = ' num2str(PARAM3(p3_ii)) ', ' titleTxt4 ' = ' num2str(PARAM4(p4_ii)) ', ' titleTxt5 ' = ' num2str(PARAM5(p5_ii))];
                    
                    spkBin = BasalAct.WT.All.spktrain;
                    baseline_fr_WT = get_avg_baseline(spkBin, CUTTIME, PhotoInjT);
                    
                    spkBin = BasalAct.KO.All.spktrain;
                    baseline_fr_KO = get_avg_baseline(spkBin, CUTTIME, PhotoInjT);
                    
                    WT_expBaselineSpk = baseline_fr_WT/1000*BurstRange;
                    WT_burstSpk = BasalAct.WT.All.BurstSpk; % BurstSpk  = Number of spike during BurstRange
                    CntBurstSpkWT(p1_ii,p2_ii) = mean(WT_burstSpk - WT_expBaselineSpk); % Number of spike during BurstRange - expected number of spike from baseline activity
                    
                    KO_expBaselineSpk = baseline_fr_KO/1000*BurstRange;
                    KO_burstSpk = BasalAct.KO.All.BurstSpk; % BurstSpk  = Number of spike during BurstRange
                    CntBurstSpkKO(p1_ii,p2_ii) = mean(KO_burstSpk - KO_expBaselineSpk); % Number of spike during BurstRange - expected number of spike from baseline activity
                    
                    
                    cnt = cnt+1;
                    
                    % Distribution of first rebound spike timing %%%%%%
                    
                    tmpRebound_WT = zeros(M1_ncells,1); tmpRebound_KO = zeros(M1_ncells,1);
                    for id = 1 : M1_ncells
                        tmp = find(BasalAct.WT.All.BurstSpkTrain(id,:),1,'first');
                        if ~isempty(tmp)
                            tmpRebound_WT(id) = tmp;
                        end
                        
                        tmp = find(BasalAct.KO.All.BurstSpkTrain(id,:),1,'first');
                        if ~isempty(tmp)
                            tmpRebound_KO(id) = tmp;
                        end
                    end
                    tmpRebound_WT = tmpRebound_WT(tmpRebound_WT ~= 0); % pick only the case that there is spike during BurstRange ( if there is no spike at all,then t = 0)
                    tmpRebound_KO = tmpRebound_KO(tmpRebound_KO ~= 0); % Contain data of first spike latency
                    
                    hbin =5;
                    histBin = 1:hbin:size(BasalAct.WT.All.BurstSpkTrain,2);
                    
                    figure(fRB);
                    
                    subplot(nR,nC,cnt);
                    [N,X] = hist(tmpRebound_WT,histBin );
                    pH=bar(X,N/M1_ncells*100,'FaceColor','k','EdgeColor','k'); hold on; %Percent of Total Neuron
                    cH = get(pH,'Children');
                    set(cH,'FaceAlpha',0.5); % 0 = transparent, 1 = opaque.
                    [N,X] = hist(tmpRebound_KO,histBin );
                    pH=bar(X,N/M1_ncells*100,'FaceColor','r','EdgeColor','r'); hold on;
                    cH = get(pH,'Children');
                    set(cH,'FaceAlpha',0.5); % 0 = transparent, 1 = opaque.
                    ylabel('% Neurons'); xlabel('Time(ms)'); xlim([0  BurstRange]);
                    title(get_Parameters_titleText(PARAMETERS, [1,2], [ p1_ii,p2_ii]))
                    
                    % Distribution of spikes in bursting range
                    
                    tmpRebound_WT = sum(BasalAct.WT.All.BurstSpkTrain); %# of spike at particular time
                    tmpRebound_KO = sum(BasalAct.KO.All.BurstSpkTrain);
                    ll =length(tmpRebound_WT);
                    if mod(ll,hbin) == 0
                        % Can use reshape fn.
                        tmpRebound_WT = sum(reshape(tmpRebound_WT,[hbin ll/hbin]));
                        tmpRebound_KO = sum(reshape(tmpRebound_KO,[hbin ll/hbin]));
                        xx = histBin;
                    else %If not, just set bin to 1
                        xx = 1:1:size(BasalAct.KO.All.BurstSpkTrain,2);
                    end
                    figure(fBspk);
                    subplot(nR,nC,cnt);
                    pH=bar(xx,tmpRebound_WT/M1_ncells*100,'FaceColor','k','EdgeColor','k'); hold on;
                    cH = get(pH,'Children');
                    set(cH,'FaceAlpha',0.5); % 0 = transparent, 1 = opaque.
                    pH=bar(xx,tmpRebound_KO/M1_ncells*100,'FaceColor','r','EdgeColor','r');
                    cH = get(pH,'Children');
                    set(cH,'FaceAlpha',0.5); % 0 = transparent, 1 = opaque.
                    ylabel('% Neurons'); xlabel('Time(ms)'); xlim([0  BurstRange]);
                    title(get_Parameters_titleText(PARAMETERS, [1,2], [ p1_ii,p2_ii]))
                    
                    clear BasalAct spkBin
                end
            end
            figure(fRB); suptitle(['M1 : Latency of first Rebound spike : ' get_Parameters_titleText(PARAMETERS, [3:5], [ p3_ii,p4_ii,p5_ii])] );  % Note equal to the latency of peak firing rate --> latency of peak firing rate calculate from instantaneous firing rate
            figure(fBspk); suptitle(['M1 : Percent of neuron that fire during bursting period : ' get_Parameters_titleText(PARAMETERS, [3:5], [ p3_ii,p4_ii,p5_ii]) ]);  % Note equal to the latency of peak firing rate --> latency of peak firing rate calculate from instantaneous firing rate
            
            if (SAVE_FIG)
                tmpTxt = get_Parameters_saveText(PARAMETERS, [1:5], [ p1_ii,p2_ii, p3_ii,p4_ii, p5_ii]);
                ffig = [ dirLoc dirFig 'FirstReboundSpike_M1_'  tmpTxt];
                saveas(fRB, [ffig '.fig'] , 'fig');
                saveas(fRB, [ffig '.jpg'] , 'jpg');
                ffig = [ dirLoc dirFig 'AllBurstSpike_M1_'  tmpTxt]; % Better minus by expected baseline activity
                saveas(fBspk, [ffig '.fig'] , 'fig');
                saveas(fBspk, [ffig '.jpg'] , 'jpg');
                
                if (Close_Fig_aftr_save)
                    close(fRB); close(fBspk);
                end
                
            end
      


%             %% Parameters Matrix  rangeTC vs. weighing Factor for normFR_atBursting
%             %----- consider using Z-score with baseline activity as mean here
%             % %// Sample data
%             % matrix = rand(10, 10);
%             % testData = rand(10, 10);
%             %
%             % %// Obtain mu and sigma
%             % mu = mean(matrix, 1);
%             % sigma = std(matrix, [], 1);
%             % %// or use: [Z, mu, sigma] = zscore(matrix);
%             %
%             % %// Convert into z-scores using precalculated mu and sigma
%             % C = bsxfun(@rdivide, bsxfun(@minus, testData, mu), sigma);


normFR_atBurstingWT= CntBurstSpkWT./BurstRange*1000; %M1's average firing rate during burst period at VL
normFR_atBurstingKO= CntBurstSpkKO./BurstRange*1000; %M1's average firing rate during burst period at VL


SAVE_FIG =1;
figLoc =[  194         333        1760         719];
tt1 = 'WT'; tt2 = 'KO';
xt = PARAM2; yt = PARAM1;
xl = titleTxt2; yl =titleTxt1; %yl =  'Average number of VL per M1 (cell)';
suptitleTxt = {'Average firing rate of M1 during bursting range  - baseline activity', get_Parameters_titleText(PARAMETERS, [3:5], [ p3_ii,p4_ii,p5_ii])} ;

[ fg_M1normFR ] = plotFigure_paramMat( normFR_atBurstingWT,normFR_atBurstingKO, tt1,tt2, xl,yl, xt,yt, suptitleTxt,figLoc );
if(SAVE_FIG)
    fg = fg_M1normFR;
       tmpTxt = get_Parameters_saveText(PARAMETERS, [1:5], [ p1_ii,p2_ii, p3_ii,p4_ii, p5_ii]);
    figname = ['M1normFR_atBursting' tmpTxt];
    saveas(fg,[dirLoc dirFig figname '.fig'],'fig');
    saveas(fg,[dirLoc dirFig figname '.jpg'],'jpg');
end
VLperM1 = numVL_M1_LST(:,1);
xt = PARAM2; yt = round(VLperM1);
xl = titleTxt2; yl =  'Average number of VL per M1 (cell)';
suptitleTxt = {'Average firing rate of M1 during bursting range  - baseline activity',get_Parameters_titleText(PARAMETERS, [3:5], [ p3_ii,p4_ii,p5_ii])};
[ fg_M1normFRnumVL ] = plotFigure_paramMat( normFR_atBurstingWT,normFR_atBurstingKO, tt1,tt2, xl,yl, xt,yt, suptitleTxt,figLoc );
if(SAVE_FIG)
    fg = fg_M1normFRnumVL;
       tmpTxt = get_Parameters_saveText(PARAMETERS, [1:5], [ p1_ii,p2_ii, p3_ii,p4_ii, p5_ii]);
    figname = ['M1normFR_atBursting_numVLM1' tmpTxt];
    saveas(fg,[dirLoc dirFig figname '.fig'],'fig');
    saveas(fg,[dirLoc dirFig figname '.jpg'],'jpg');
end

%% M1 raw activity

M1_burstSpkWT = zeros(length(PARAM1),length(PARAM2));
M1_burstSpkKO = zeros(length(PARAM1),length(PARAM2));
for p1_ii = 1 : length(PARAM1)
    for p2_ii = 1 : length(PARAM2)
        BasalAct = ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii}.M1; %% for M1
        M1_burstSpkWT(p1_ii,p2_ii) = mean(BasalAct.WT.All.BurstSpk)./BurstRange*1000;
        M1_burstSpkKO(p1_ii,p2_ii) = mean(BasalAct.KO.All.BurstSpk)./BurstRange*1000;
    end
end

SAVE_FIG = 1;
xt = PARAM2; yt = PARAM1;
xl = titleTxt2; yl =titleTxt1;
suptitleTxt = {'Average firing rate of M1 during bursting range',get_Parameters_titleText(PARAMETERS, [3:5], [ p3_ii,p4_ii,p5_ii])};
[ fg_M1FR ] = plotFigure_paramMat( M1_burstSpkWT,M1_burstSpkKO, tt1,tt2, xl,yl, xt,yt, suptitleTxt,figLoc );
if(SAVE_FIG)
    fg = fg_M1FR;
       tmpTxt = get_Parameters_saveText(PARAMETERS, [1:5], [ p1_ii,p2_ii, p3_ii,p4_ii, p5_ii]);
    figname = ['M1FR_atBursting' tmpTxt];
    saveas(fg,[dirLoc dirFig figname '.fig'],'fig');
    saveas(fg,[dirLoc dirFig figname '.jpg'],'jpg');
end

xt = PARAM2; yt = round(VLperM1);
xl = titleTxt2; yl =  'Average number of VL per M1 (cell)';
suptitleTxt = {'Average firing rate of M1 during bursting range', get_Parameters_titleText(PARAMETERS, [3:5], [ p3_ii,p4_ii,p5_ii])};
[ fg_M1FRnumVL ] = plotFigure_paramMat( M1_burstSpkWT,M1_burstSpkKO, tt1,tt2, xl,yl, xt,yt, suptitleTxt,figLoc );
if(SAVE_FIG)
    fg = fg_M1FRnumVL;
       tmpTxt = get_Parameters_saveText(PARAMETERS, [1:5], [ p1_ii,p2_ii, p3_ii,p4_ii, p5_ii]);
    figname = ['M1FR_atBurstingNumVL' tmpTxt];
    saveas(fg,[dirLoc dirFig figname '.fig'],'fig');
    saveas(fg,[dirLoc dirFig figname '.jpg'],'jpg');
end

      
        end
    end
end


%%
if(SAVE_WORKSPACE)
    save([dirLoc dirFig 'Saved_Workspace_result_' SPECIFIED_BASAL_ACT_CODENAME date '.mat' ],'-v7.3');
end