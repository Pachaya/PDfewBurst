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
SAVE_BASAL_ACT = 1; SPECIFIED_BASAL_ACT_CODENAME = 'test';
SAVE_FIG = 1; Close_Fig_aftr_save = 0;
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

%PDfewBurst_GPmVLmd1_0del_rTC120_wmTC10_WT_GPmInput_Amp0.3_Dur1000_GPmVLw_m0.06_sig0.01_InGauss0.2_IGmean0_IGmeanSig0_W0.029_SpecifiedPoisSpk_sig0.00Hz_oscF20.00Hz_phase0_T1000_trial1
% Parameters setting
NUM_TRIAL = 5;
ncells = 1150;
M1_ncells = 166;
TRIAL_LST = 1 : NUM_TRIAL;

rTC_LST = [120]; %
wmTC_LST = [10]; %[10 25 50 75 100];

LightAmp_LST = [0.3 ];
GPmVLw_mean_LST = [0.06];
GPmVLw_sig_LST =[ 0.01 ];

OSC_F_LST = [10];

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
PARAM5 = OSC_F_LST;
lblTxt5 = 'Oscillation Frequency';
saveTxt5 = 'OSC_F';
titleTxt5 = 'Osc Freq';
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
dirFig = 'TestOsc';
%dirFig = ['Fig_TCconvergences_big' get_Parameters_RangeTxt( PARAMETERS,[1,2]) '/'];
mkdir([dirLoc dirFig])




%Start donwnload simulation data
for p1_ii = 1 : length(PARAM1)
    for p2_ii = 1 : length(PARAM2)
        for p3_ii = 1 : length(PARAM3)
            for p4_ii = 1 : length(PARAM4)
                for p5_ii = 1 : length(PARAM5)
                    
                    r_ii = p1_ii; wm_ii = p2_ii;
                    la_ii = p3_ii;  m_ii = p4_ii; s_ii =1; of_ii = p5_ii;
                    TRIAL_NO = 1;
                    
                  for cell_type = 1 : 2
               
                        tmp = [];
                        if (cell_type == 1)
                            cTxt = 'WT';
                        elseif (cell_type == 2)
                            cTxt = 'KO';
                        end
                        
                        coreFileName = 'PDfewBurst_GPmVLmd1_test' ;
                        
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
                        phase =0;

                        osc_f = OSC_F_LST(of_ii);
                        Simulation_Code = [coreFileName '_rTC' num2str(rTC) '_wmTC' num2str(wmTC) '_' cTxt '_' 'GPmInput_Amp' num2str(GPmLight) '_Dur' num2str(GPmLightDur) '_GPmVLw_m' num2str(GPm_w_mn) '_sig' num2str(GPm_w_sg) ...
                            '_InGauss' num2str(InGauss_STDEV) '_IGmean' num2str(NoiseMEAN) '_IGmeanSig' num2str(IGmeanSig) '_W' num2str(W_Weight) '_SpecifiedPoisSpk_sig0.00Hz_oscF' num2str(osc_f) '.00Hz_phase0_T' num2str(TSTOP) '_trial' num2str(TRIAL_NO)];

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
%% Single Cell average --> VL # 385
% Mean 20 OSc Hz = 50 
% Mean 10 / Amp 5  
%%  FFT 
SIGMA_fr_All = 10;

 
    T = TSTOP;
    tt = 0:1:T;
    %%% for VL
                    tmpWTtrain_VL =  ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii}.VL.WT.All.spktrain;
                    tmpKOtrain_VL =  ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii}.VL.KO.All.spktrain;
%                     Compare WT and KO
%                     avgFR_WT_VL = getInstantaneousFiringRate(tmpWTtrain_VL, SIGMA_fr_All, RES);
%                     avgFR_KO_VL = getInstantaneousFiringRate(tmpKOtrain_VL, SIGMA_fr_All, RES);

%     Efmu = avgFR_KO_VL;
%     Ifmu = avgFR_WT_VL;
%     figure; subplot(311); plot(Efmu,'r'); hold on; plot(Ifmu,'k');
%     
        % Spike Train
        [KOci, KOti] = find(tmpKOtrain_VL);
        [WTci, WTti] = find(tmpWTtrain_VL);
        
        % mean firing rate
        KOfmu = mean(tmpKOtrain_VL, 1)*1000;
        WTfmu = mean(tmpWTtrain_VL, 1)*1000;
        fosc = figure;  set(fosc, 'position', [496         763        1042         420]);
        subplot(211);
        plot(WTfmu,'k');   hold on;  plot(KOfmu,'r'); legend('WT','KO'); title('Raw instantaneous average firing rate')
         
        gsig = 3;   % filter [ms]
        gfil = exp( -[-round(4*gsig):round(4*gsig)].^2 /2/gsig^2);
        gfil = gfil/sum(gfil(:));
        
        KOfmus = conv(KOfmu, gfil, 'same');
        WTfmus = conv(WTfmu, gfil, 'same');
        subplot(212); plot(WTfmus,'k');   hold on;  plot(KOfmus,'r'); legend('WT','KO'); title('Smooth instantaneous average firing rate')
        % FFT
        KOfft = abs(fftshift(fft(KOfmus)));
        WTfft = abs(fftshift(fft(WTfmus)));
%         subplot(313); plot(WTfft,'k');   hold on;  plot(KOftt,'r'); legend('WT','KO'); title('FFT')
       
        dHz = 1/T * 1000;
        Hz_ind = [-round(T/2):round(T/2)]*dHz;
%         Hz_ind = Hz_ind(1:T);
         figure; plot(Hz_ind, WTfft,'k'); hold on;  plot(Hz_ind, KOfft,'r');  legend('WT','KO');
         xlim([1 100]); ylabel('FFT amplitude (a.u.)'); xlabel('Frequency (Hz)'); title('raw FFT')

        [Hind, KO_FFT, WT_FFT] = analyze_fft(Hz_ind, WTfft, KOfft);
        figure; 
        plot(Hind, WT_FFT.FF,'k') ;     hold on; plot(Hind, KO_FFT.FF,'r');  legend('WT','KO');
        xlim([1 100]); ylabel('FFT amplitude (a.u.)'); xlabel('Frequency (Hz)'); title('Smooth FFT');

        

%% Plot the instantaneous average firing rate


SIGMA_fr_All = 30;
Close_Fig_aftr_save = 0;

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
                    
                    if (SAVE_FIG)
                        tmpTxt = get_Parameters_saveText(PARAMETERS, [1:5], [p1_ii p2_ii p3_ii p4_ii p5_ii]);
                        ffig = [ dirLoc dirFig 'Inst_Avg_FR_' tmpTxt ];
                        saveas(  fFR, [ffig '.jpg'], 'jpg');        saveas(  fFR, [ffig '.fig'], 'fig');
                        if (Close_Fig_aftr_save)
                            close(fFR)
                        end
                    end
                    clear  tmpWTtrain_VL  tmpKOtrain_VL tmpWTtrain_M1  tmpKOtrain_M1;
                end
            end
        end
    end
end


%% Test avg FR



SIGMA_fr_All = 30;
p1_ii =1; p2_ii =1; p3_ii =1; p4_ii =1; p5_ii =1; 
                    % RasterPlot 
                    
                    %%% for VL
                    tmpWTtrain_VL =  ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii}.VL.WT.All.spktrain;
             
                    %Compare WT and KO
                    avgFR_WT_VL = getInstantaneousFiringRate(tmpWTtrain_VL, SIGMA_fr_All, RES);
                   
                    
                    
                    %Plot and save
                    fFR = figure; set(gcf, 'position', [ 624   256   733   722]); set(gcf,'PaperPositionMode','auto');  
                    plot(avgFR_WT_VL,'k');  hold on;   
                    title('VL');   ylabel('Average firing rate (Hz)');            xlabel('Time(ms)'); 
                    
                    tmpTxt = get_Parameters_titleText(PARAMETERS, [1:5], [p1_ii p2_ii p3_ii p4_ii p5_ii]);
                    suptitle({'Average firing rate Fr(t)', tmpTxt});
       %% Raster Plot 
    SAVE_RASTER_FIG = 1;
p1_ii =1; p2_ii =1; p3_ii =1; p4_ii =1; p5_ii =1; 

            
            BasalAct = ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii}.VL;
            
            simTxt = ['VL ' titleTxt3 ' = ' num2str(PARAM3(p3_ii)) ', ' titleTxt4 ' = ' num2str(PARAM4(p4_ii)) ', ' titleTxt5 ' = ' num2str(PARAM5(p5_ii))];
            
            spkBin = BasalAct.WT.All.spktrain;
            baseline_fr_WT = get_avg_baseline(spkBin, CUTTIME, PhotoInjT);
            [freq_all_WT,spkTime_all_WT, fg_handle_WT ] = raster_from_spkbin_BurstRange( spkBin,PhotoStop, Tstop,PhotoInjT, PhotoStop, DelayT, BurstRange, [ simTxt ' : WT']);
            set(fg_handle_WT, 'position',[  449   450   791   528])
            
            spkBin = BasalAct.KO.All.spktrain;
            baseline_fr_KO = get_avg_baseline(spkBin, CUTTIME, PhotoInjT);
            [freq_all_KO,spkTime_all_KO, fg_handle_KO ] = raster_from_spkbin_BurstRange( spkBin,PhotoStop, Tstop,PhotoInjT, PhotoStop, DelayT, BurstRange, [ simTxt ' : KO']);
            set(fg_handle_KO, 'position',[  449   450   791   528])
            
            if (SAVE_RASTER_FIG)
                tmpTxt = get_Parameters_saveText(PARAMETERS, [3:5], [ p3_ii p4_ii p5_ii]);
                ffig = [ dirLoc dirFig 'RasterPlot_' tmpTxt];
                ttt = tic(); saveas( fg_handle_WT, [ffig '_WT.jpg'], 'jpg'); toc(ttt);
                ttt = tic(); saveas( fg_handle_WT, [ffig '_WT.fig'], 'fig'); toc(ttt);
                saveas( fg_handle_KO, [ffig '_KO.jpg'], 'jpg')
                saveas( fg_handle_KO, [ffig '_KO.fig'], 'fig')
                if (Close_Fig_aftr_save)
                   close(fg_handle_WT);                        close(fg_handle_KO);
                end
                
            end
     