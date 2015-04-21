
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
SHOW_RASTER = 0;
SAVE_WORKSPACE = 0;

DoTTest = 1;      DoSampleTtest = 1;
SAVE_NEURON_ACT = 1;
FNAME_SIM = 'Min150324' ;
plotSampleV = 0;
SPECIFIED_RATIO = 0; 
Basal_Act_set  = [];

NoiseMEAN = 0;

SAVE_FIG = 1;

dirLoc = 'Model1_MinParam150324/' ;
rEE_rEI_rIE_rII = [1 1 1 1];
w_MULT = 125;


pulseHz = 1; N_pulse = 1;
Input_I_amp_list = [0];
Input_I_dur_list = [0];
% Input_I_amp = -0.3; %[-1; -2; -3; -4; -5;];% [-1; -5]; %
% Input_I_dur = 5; %[10; 20; 30; 40; 50;];
ReboundPeakAmp_WTlist = zeros(length(Input_I_amp_list),length(Input_I_dur_list));
ReboundPeakDel_WTlist = zeros(length(Input_I_amp_list),length(Input_I_dur_list));
ReboundPeakAmp_KOlist = zeros(length(Input_I_amp_list),length(Input_I_dur_list));
ReboundPeakDel_KOlist = zeros(length(Input_I_amp_list),length(Input_I_dur_list));
Save_BasalAct =  cell(length(Input_I_amp_list),length(Input_I_dur_list));


ADD_I_to_M1 = 0;
ADD_I_to_VL = 0;

FIG_ALL = 0;
dirFig = 'Fig/';
NoiseMEAN_WTKO = [0; 0;]; %[0.0291; 0.0075;]; %[0.0295; 0.016]; %%%%%%%%%%%%%%%%%%%%%%%%%%% [ WT; KO;] for 10 Hz -> [0.0295; 0.016]; for 5 Hz -> [0.0291; 0.0075;]; 5Hz 0.0318 0.3  for WT 5Hz (2) -> [0.0318; 0.0075;]
NoiseSIG_WTKO  =  [0.3; 0.3;]; %[0.8; 0.8;];  %[0.24; 0.3];  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  [WT; KO;] for 10 Hz -> [0.4; 0.5]; for 5 Hz -> [0.24; 0.3];  5Hz 0.0318 0.3  for WT 5Hz (2) -> [0.3; 0.3;]
VARIE_MEAN_InGauss = 1;
NoiseMEANsigma_WTKO = [0; 0;]; %%%%% [WT; KO;]
TSTOP = 4000;

rTC_LST = [120:10:150];
wmTC_LST = [7;];
rTC = 250;
LightAmp_LST = [1.5];
GPmVLw_mean_LST = [1];
GPmVLw_sig_LST =[0.05; 0.1; 0.2; 0.3];

GPmLightDur_LST = [5, 25, 50, 100];

CUTTIME = 500;
BurstRange = 100;
DelayT = 0;
PhotoInjT = 1000;


PARAM1 = LightAmp_LST;
legTxt1 = 'Light Stimulus amplitude (nA)'; %'E_L '; %'soma size';
titleTxt1 = 'I_A_m_p';
saveTxt1 = 'GPmAmp';
PARAM2 = GPmLightDur_LST; %GPmVLw_mean_LST;
legTxt2 = 'Light Stimulus Amplitude (nA)'; % 'GPm-VL weight mean';
saveTxt2 = 'GPmDur'; %'GPmVLw_m';
titleTxt2 = 'I_D_u_r';
PARAM3 = GPmVLw_sig_LST;
legTxt3 = 'GPm-VL weight sigma';
saveTxt3 = 'GPmVLw_sig';
titleTxt3 = 'W sig';

ACT_Record = cell(length(PARAM1),length(PARAM2),length(PARAM3));
N_Param = 3;
for ii = 1 : N_Param
    PARAMETERS{ii}.PARAM = eval(sprintf('PARAM%d',ii)) ;
    PARAMETERS{ii}.legTxt = eval(sprintf('legTxt%d',ii));
    PARAMETERS{ii}.titleTxt = eval(sprintf('titleTxt%d',ii));
    PARAMETERS{ii}.saveTxt = eval(sprintf('saveTxt%d',ii));
end

 
% for p1_ii = 1 : length(PARAM1)
%     for p2_ii = 1 : length(PARAM2)
%         for p3_ii = 1 : length(PARAM3)
%             la_ii = p1_ii; m_ii = 1; ld_ii = p2_ii; s_ii = p3_ii;

NUM_TRIAL = 5;
ACT_Record = cell(NUM_TRIAL,1);
ncells = 1150;
for TRIAL_NO = 1 :NUM_TRIAL
    for cell_type = 1 : 2
        
        if (cell_type == 1)
            cTxt = 'WT';
        elseif (cell_type == 2)
            cTxt = 'KO';
        end
        
        coreFileName = 'PDfewBurst_GPmVLmd1' ;
        
        InGauss_STDEV = 0.2; %0.2;, 0.3
        NoiseMEAN = -0.15;
        IGmeanSig = 0;
        W_Weight = 0.0012;
        W_sig = 0;
        PoisInputFr = 0;
        PoisInputSig = 0;
        TSTOP = 10000;
        %         Rand123Seed_ID1_PDfewBurst_GPmVLmd1_WT_InGauss0.13_IGmean-0.2_IGmeanSig0_W0.001_SpecifiedPoisSpk_sig0.00Hz_T2000_trial4
        %                 GPmLightDur = GPmLightDur_LST(ld_ii);
        %                 PhotoStop = PhotoInjT+GPmLightDur;
        %                 PhotoStop = PhotoStop;
        
        %                 GPmLight = LightAmp_LST(la_ii);
        %                 GPm_w_mn = GPmVLw_mean_LST(m_ii);
        %                 GPm_w_sg = GPmVLw_sig_LST(s_ii);
        
        
        Name_postfix = [coreFileName '_' cTxt '_InGauss' num2str(InGauss_STDEV) '_IGmean' num2str(NoiseMEAN) '_IGmeanSig' num2str(IGmeanSig) '_W' num2str(W_Weight) '_SpecifiedPoisSpk_sig' num2str(PoisInputSig)  '.00Hz_T' num2str(TSTOP) '_trial' num2str(TRIAL_NO) ];
%         Name_postfix = [coreFileName '_' cTxt  '_IGmean' num2str(NoiseMEAN) '_IGmeanSig' num2str(IGmeanSig) '_W' num2str(W_Weight) '_Wsig' num2str(W_sig) '_SpecifiedPoisSpk_sig' num2str(PoisInputSig)  '.00Hz_T' num2str(TSTOP) '_trial' num2str(TRIAL_NO) ];
        % PDfewBurst_GPmVLmd1_KO_IGmean-0.15_IGmeanSig0_W0.0015_SpecifiedPoisSpk_sig0.00Hz_T10000_trial1
        
        
        
        disp('==================================================================================================')
        disp(Name_postfix)
        if (cell_type == 1)
            cTxt = 'WT';
        else
            cTxt = 'KO';
        end
        %                 figNameCode = [ saveTxt1 num2str(PARAM1(p1_ii))  '_' saveTxt2 num2str(PARAM2(p2_ii)) '_' cTxt];
        figNameCode = [cTxt '_Trial_' num2str(TRIAL_NO) ];
        %                 CheckFileExist( dirLoc, Name_postfix  )
        
        VL_BaselineAct
        
        
        tmp.E = E;      tmp.I = I;  tmp.All = All; %tmp.VL = VL; tmp.M1 = M1;
        if (cell_type == 1)
            WT = tmp;
            cTxt = 'WT';
        elseif (cell_type == 2)
            KO = tmp;
            cTxt = 'KO';
        end
        if (plotSampleV)
            %soma's volt
            [NEcell,NIcell, Tstop,somaVall] = get_save_vec_file_BigNet([dirLoc 'SomaVolt_' Name_postfix '.txt'] );
            samE = 385;
            samI = 129;
            samTstr = 1000; samTstp = 1800; %500+ Input_I_dur +200;
            fgSmpl = figure; plot(samTstr:samTstp, somaVall(samE,samTstr:samTstp), 'r'); hold on;
            plot(samTstr:samTstp, somaVall(samE+1,samTstr:samTstp), 'm');
            if(ADD_I_to_VL)
                plot(samTstr:samTstp,somaVall(N_E+samI,samTstr:samTstp),'b');
                plot(samTstr:samTstp,somaVall(N_E+samI+1,samTstr:samTstp),'k');
            end
            title([cTxt ': sample membrane potential' ])
            if(SAVE_FIG)
                saveas(fgSmpl, [dirLoc dirFig 'SamV' figNameCode '.fig'], 'fig')
                saveas(fgSmpl, [dirLoc dirFig 'SamV' figNameCode '.jpg'], 'jpg')
            end
        end
        
    end
    
    Basal_Act.INOISE.MEAN = NoiseMEAN;
    Basal_Act.INOISE.STDEV = InGauss_STDEV;
    
    Basal_Act.WT = WT;
    Basal_Act.KO = KO;
 
    if (DoTTest)
        % two-tailed t-test % during normal
        
        % disp('PoisSPk')
        % [hP,pP] = ttest(WT.PoisSpk.fr_data, KO.PoisSpk.fr_data)
        disp('--------------------------------------------------------------------------------------------------')
        disp('E cells : two-tailed t-test')
        [hE,pE] = ttest2(WT.All.fr_data, KO.All.fr_data, 'Tail','both');
        disp([ ' p-value = ' num2str(pE)])
        if(hE)
            disp('H0 rejected : average normal firing rate of E cells in WT significantly different from KO (p < 0.05)')
        else
            disp('Do not rejected H0 : average normal firing rate of E cells in WT does not significantly different from KO (p > 0.05)')
        end
        Basal_Act.ttest2.hE = hE;
        Basal_Act.ttest2.pE = pE;
      
        
    end
    %% Sample T-Test

    if(DoSampleTtest)
        
        Nsample = 15;
        Nrepeat = 100;
        disp(['===== Sample Test : #sample = ' num2str(Nsample) ' for ' num2str(Nrepeat) ' trials ====='])
        tresult = cell(Nrepeat,1);
        plist = zeros(Nrepeat,1);
        hlist = zeros(Nrepeat,1);
        tresult_t = cell(Nrepeat,1);
        plist_t = zeros(Nrepeat,1);
        hlist_t = zeros(Nrepeat,1);
 
        
        for ii = 1: Nrepeat
            sample = randperm(ncells, Nsample);
            [ht,pt] = ttest2(WT.All.fr_data(sample), KO.All.fr_data(sample));
            [h,p] = ranksum(WT.All.fr_data(sample), KO.All.fr_data(sample));
            
            tresult{ii}.sampleID = sample;
            tresult{ii}.h = h;
            tresult{ii}.p = p;
            %     disp(['h:' num2str(h) ', p:' num2str(p)])
            plist(ii) = p;
            hlist(ii) = h;
            
            tresult_t{ii}.sampleID = sample;
            
            tresult_t{ii}.h = ht;
            tresult_t{ii}.p = pt;
            %     disp(['ht:' num2str(ht) ', p:' num2str(pt)])
            plist_t(ii) = pt;
            hlist_t(ii) = ht;
        end
        
        disp('##Ranksum##')
        %  sum((plist>0.05))
        disp([num2str(sum((plist>0.05))) ' cases do not rejected  H0 : average normal firing rate of all cells in WT does not significantly different from KO (p > 0.05)'])
        disp('##T-test#')
        %  sum((plist_t>0.05))
        disp([num2str(sum((plist_t>0.05))) ' cases do not rejected  H0 : average normal firing rate of all cells in WT does not significantly different from KO (p > 0.05)'])
        % Common cases
        disp('The other cases are rejecting H0 : average normal firing rate of all cells in WT significantly different from KO (p < 0.05)')
        disp('')
        disp('common cases in Ranksum and T-test')
        find((plist>0.05) == (plist_t>0.05))
        Basal_Act.sampleTest.ranksum = tresult;
        Basal_Act.sampleTest.ttest = tresult_t;
        
    end
    
    %             ACT_Record{p1_ii,p2_ii,p3_ii } = Basal_Act;
    ACT_Record{TRIAL_NO} = Basal_Act;
    
end

sample_data = cell(NUM_TRIAL,1);

for TRIAL_NO = 1 : NUM_TRIAL
Basal_Act = ACT_Record{TRIAL_NO};
WT = Basal_Act.WT; KO = Basal_Act.KO;
Nrepeat = size(Basal_Act.sampleTest.ttest,1);
Nsample = length(Basal_Act.sampleTest.ttest{1}.sampleID );
WT_sample = zeros(Nrepeat,2); KO_sample = zeros(Nrepeat,2);

for ii = 1 : Nrepeat
       sample = Basal_Act.sampleTest.ttest{ii}.sampleID;
       WT_sample(ii,:) =  [mean( WT.All.fr_data(sample))  std( WT.All.fr_data(sample))];
       KO_sample(ii,:) =  [mean( KO.All.fr_data(sample)) std( KO.All.fr_data(sample))];
end
sample_data{TRIAL_NO}.WTsample = WT_sample;
sample_data{TRIAL_NO}.KOsample = KO_sample;
end

%%
if (SAVE_NEURON_ACT)
    save([dirLoc dirFig 'Result_ACT_Record_' coreFileName '.mat'], 'ACT_Record','PARAMETERS','-v7.3');
    save(['Result_ACT_Record_' coreFileName '.mat'], 'ACT_Record','PARAMETERS','-v7.3');
end
%%
if(SHOW_RASTER)
    for TRIAL_NO = 1 : NUM_TRIAL
        simTxt = ['TRIAL NO =' num2str(TRIAL_NO)];
        spkBin = ACT_Record{TRIAL_NO}.WT.All.spktrain;
        [freq_all_WT,spkTime_all_WT, fg_handle_WT ] = raster_from_spkbin( spkBin,cuttime, Tstop, [ simTxt ' : WT']);
        set(fg_handle_WT, 'position',[  449   450   791   528])
        spkBin = ACT_Record{TRIAL_NO}.KO.All.spktrain;
        [freq_all_KO,spkTime_all_KO, fg_handle_KO ] = raster_from_spkbin( spkBin,cuttime, Tstop, [ simTxt ' : KO']);
        set(fg_handle_KO, 'position',[  449   450   791   528])
        disp(['########### Trial = ' num2str(TRIAL_NO) ])
        disp(['### WT avg freq = ' num2str(mean(freq_all_WT)) ' Hz, std = ' num2str(std(freq_all_WT))])
        disp(['### KO avg freq = ' num2str(mean(freq_all_KO)) ' Hz, std = ' num2str(std(freq_all_KO))])
        if (SAVE_FIG)
            %ffig = [ dirLoc dirFig 'RasterPlot_' saveTxt1 num2str(PARAM1(p1_ii)) '_' saveTxt2 num2str(PARAM2(p2_ii)) '_' saveTxt3 num2str(PARAM3(p3_ii))];
            ffig = [ dirLoc dirFig 'RasterPlot_Trial' num2str(TRIAL_NO)];
            saveas( fg_handle_WT, [ffig '_WT.jpg'], 'jpg')
            saveas( fg_handle_WT, [ffig '_WT.fig'], 'fig')
            saveas( fg_handle_KO, [ffig '_KO.jpg'], 'jpg')
            saveas( fg_handle_KO, [ffig '_KO.fig'], 'fig')
        end
        figure;
        subplot(121); hist(ACT_Record{TRIAL_NO}.WT.All.fr_data,15); title(['WT : Trial #' num2str(TRIAL_NO)])
        subplot(122); hist(ACT_Record{TRIAL_NO}.KO.All.fr_data,15); title(['KO : Trial #' num2str(TRIAL_NO)])
        
    end
end

for TRIAL_NO = 1:NUM_TRIAL
    ff=figure; set(ff, 'position',[ 725   514   792   298]);  set(ff, 'PaperPositionMode','Auto');
    subplot(121); hist(ACT_Record{TRIAL_NO}.WT.All.fr_data,15); title(['WT : Trial #' num2str(TRIAL_NO)])
    subplot(122); hist(ACT_Record{TRIAL_NO}.KO.All.fr_data,15); title(['KO : Trial #' num2str(TRIAL_NO)])
    suptitle('Distribution of baseline activity (Hz)');
    saveas( ff, [ dirLoc dirFig 'FrHist' '_TRIAL' num2str(TRIAL_NO) '.jpg'], 'jpg')
    saveas( ff, [ dirLoc dirFig 'FrHist' '_TRIAL' num2str(TRIAL_NO) '.fig'], 'fig')    
end


%% compare with experimental data
% Exp Data
WTWT_MS = [26.86	28.03281345
    17.2	26.53013345
    3.665	12.27952048
    5.99	15.57402385
    5.62	22.21885032
    1.16	8.547720921
    19.37	27.93744195
    11.8	25.00476264
    28.4	32.55107443
    5.225	14.83667385
    22.995	28.85713199
    14.775	21.87142359
    4.78	16.33355283
    12.175	23.5040088
    6.32	19.61901407
    14.965	22.83668024
    18.18	24.4532686
    21.63	32.89191899
    7.87	18.27853558
    6.92	15.97953435
    4.13	19.17782766
    9.17	18.7284528
    4.415	18.07618927
    5.74	15.39746847
    7.095	19.35737104
    10.425	29.78641098
    3.435	13.16139916
    9.645	19.07929447
    8.73	19.21023761
    6.93	16.2483261
    13.2962963	25.50678635
    11.58201058	22.93539692
    32.76719577	39.36207525
    11.93650794	23.02569219
    4.835	14.73940064
    0.955	8.514519195
    ];

KOKO_MS = [7.075	16.29041274
    10.605	20.72167099
    3.51	12.66885798
    1.18	7.836789962
    8.945	21.08654916
    12.425	26.52371229
    12.68	22.36247102
    12.02	23.79055167
    16.23	25.20846154
    25.34	33.11225532
    17.15	26.69768896
    6.435	17.14723093
    19.095	31.36405487
    4.305	13.9172981
    21.915	31.42699846
    26.81	35.46521003
    15.015	25.60530666
    8.81	19.84018912
    ];

ff= figure; set(gcf, 'position',[ 725   514   792   298]);  set(gcf, 'PaperPositionMode','Auto');
subplot(121); hist(WTWT_MS(:,1),15); title(['WT : Experimental data'])
subplot(122); hist(KOKO_MS(:,1),15); title(['KO : Experimental data'])        
        
sct_fg = figure; scatter(WTWT_MS(:,1),WTWT_MS(:,2),'k');
hold on; scatter(KOKO_MS(:,1),KOKO_MS(:,2),'k');
xlabel('Mean average firing rate'); ylabel('sigma of mean firing rate');
tmp_MS = [WTWT_MS; KOKO_MS];
 p = polyfit(tmp_MS(:,1),tmp_MS(:,2),2);
    xx = 0:1:35;
    yy = polyval(p,xx);
    figure(sct_fg); hold on;
    plot(xx,yy,'color','k','LineWidth',2)
   disp(['Fitting Curve for experimental data:  y = ' num2str(p(1)) 'x^2 +' num2str(p(2)) 'x +' num2str(p(3)) ]) 
binsize = 25;

ccc = [0.5:(1-0.5)/(NUM_TRIAL-1):1];
polyfit_p = cell(NUM_TRIAL,3);

for TRIAL_NO = 1:NUM_TRIAL
     disp(['Fitting Curve for Trial#' num2str(TRIAL_NO)]) 
    binT = cuttime:binsize:Tstop;
    midT = binT(1:end-1) + binsize/2;
    tmp_spk = ACT_Record{TRIAL_NO}.WT.All.spktrain ; 
    
    [tmpWT_MS,segment_fr_dataWT] = get_mean_std_fr_of_segment_fromSpkTrain( tmp_spk, binT, binsize );
    tmp_spk = ACT_Record{TRIAL_NO}.KO.All.spktrain ;
    [tmpKO_MS, segment_fr_dataKO] = get_mean_std_fr_of_segment_fromSpkTrain( tmp_spk, binT, binsize );
    
    figure(sct_fg);
    % scatter(tmpWT_MS(:,1),tmpWT_MS(:,2),'MarkerEdgeColor',[0 0 ccc(TRIAL_NO)],'MarkerFaceColor',[0 0 ccc(TRIAL_NO)]);
    scatter(tmpWT_MS(:,1),tmpWT_MS(:,2),'MarkerEdgeColor',[0 0 ccc(TRIAL_NO)]);
    p = polyfit(tmpWT_MS(:,1),tmpWT_MS(:,2),2);
    
    yy = polyval(p,xx);
    figure(sct_fg); hold on;
    plot(xx,yy,'color',[0 ccc(TRIAL_NO) 0],'LineWidth',2)
    polyfit_p{TRIAL_NO,1} = p;
    disp(['WT:  y = ' num2str(p(1)) 'x^2 +' num2str(p(2)) 'x +' num2str(p(3)) ]) 
%     k = waitforbuttonpress;
    
    
    % hold on; scatter(tmpKO_MS(:,1),tmpKO_MS(:,2),'MarkerEdgeColor',[ccc(TRIAL_NO) 0 0],'MarkerFaceColor',[ccc(TRIAL_NO) 0 0]);
    hold on; scatter(tmpKO_MS(:,1),tmpKO_MS(:,2),'MarkerEdgeColor',[ccc(TRIAL_NO) 0 0]);    
    p2 = polyfit(tmpKO_MS(:,1),tmpKO_MS(:,2),2); % 2 for power 2 , 3- for power3
   
    yy = polyval(p2,xx);
    figure(sct_fg); hold on;
    plot(xx,yy,'color',[ccc(TRIAL_NO) ccc(TRIAL_NO) 0],'LineWidth',2)
    polyfit_p{TRIAL_NO,2} = p2;
    disp(['KO:  y = ' num2str(p2(1)) 'x^2 +' num2str(p2(2)) 'x +' num2str(p2(3)) ]) 
%     k = waitforbuttonpress; 
    tmp_MS = [tmpWT_MS; tmpKO_MS];
    p3 = polyfit(tmp_MS(:,1),tmp_MS(:,2),2); % 2 for power 2 , 3- for power3
   
    yy = polyval(p3,xx);    
    polyfit_p{TRIAL_NO,3} = p3;
    disp(['ALL: y = ' num2str(p3(1)) 'x^2 +' num2str(p3(2)) 'x +' num2str(p3(3)) ]) 
end
yyy = -0.0241.*xx.*xx + 1.4433.*xx + 4.7504; %parameter search result
plot(xx, yyy,'k:','linewidth',2)
title(['Distribution of Mean and sigma of real experiment data : bin size = ' num2str(binsize)])
saveas( sct_fg, [ dirLoc dirFig 'FrMeanSig3' '_binsize' num2str(binsize) '.jpg'], 'jpg')
saveas( sct_fg, [ dirLoc dirFig 'FrMeanSig3' '_binsize' num2str(binsize) '.fig'], 'fig')    


for TRIAL_NO = 1:NUM_TRIAL
    disp(['Trial#' num2str(TRIAL_NO)])
    tmp = ACT_Record{TRIAL_NO}.WT.All.fr_data;
    disp(['WT: mean = ' num2str(mean(tmp)) ' , std = ' num2str(std(tmp))]);
    tmp = ACT_Record{TRIAL_NO}.KO.All.fr_data;
    disp(['KO: mean = ' num2str(mean(tmp)) ' , std = ' num2str(std(tmp))]);
end

%% Check Input - Output Relationship
 
xx = 0 :1:35;
yfit = 0.0084.*xx.*xx.*xx - 0.2176.*xx.*xx + 24.865.*xx + 39.196; % Fit data
ff = figure; 
plot(xx,yfit,'k.-');
xlabel('Output Firing rate'); ylabel('Input Firing rate');
hold on;
LEG{1} = 'Fitting equation';
for TRIAL_NO = 1:NUM_TRIAL
    title(num2str(TRIAL_NO))
    WT_in = dlmread(['PoisInputFR150324_WT_' num2str(TRIAL_NO) '.txt']); WT_in = WT_in(2:end);
    tmp = ACT_Record{TRIAL_NO}.WT.All.fr_data;
    scatter(tmp,WT_in,'MarkerEdgeColor',[0 0 ccc(TRIAL_NO)]);
    LEG{2*TRIAL_NO} = ['WT, trial#' num2str(TRIAL_NO)];
    k = waitforbuttonpress;
    
    KO_in = dlmread(['PoisInputFR150324_KO_' num2str(TRIAL_NO) '.txt']); KO_in = KO_in(2:end);
    tmp = ACT_Record{TRIAL_NO}.KO.All.fr_data;
    scatter(tmp,KO_in,'MarkerEdgeColor',[ccc(TRIAL_NO) 0 0]);
    LEG{2*TRIAL_NO+1} = ['KO, trial#' num2str(TRIAL_NO)];
    k = waitforbuttonpress;
end
title('Output-Input Relationship')
legend(LEG,'location','best')


%%
if(0)
    
    CntBurstSpkWT = zeros(length(PARAM1),length(PARAM2), length(PARAM3));
    CntBurstSpkKO = zeros(length(PARAM1),length(PARAM2), length(PARAM3));
    ChanceOfBurstingWT = zeros(length(PARAM1),length(PARAM2), length(PARAM3));
    ChanceOfBurstingKO = zeros(length(PARAM1),length(PARAM2), length(PARAM3));
    for p1_ii = 1 : length(PARAM1)
        fLA1 = figure; set(fLA1,'position',[302          49        1290         948]) ;  set(gcf,'PaperPositionMode','auto')
        fLA2 = figure; set(fLA2,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
        fB1 = figure; set(fB1,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
        fB2 = figure; set(fB2,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
        cnt = 0;
        for p2_ii = 1 : length(PARAM2)
            for p3_ii = 1 : length(PARAM3)
                la_ii = p1_ii; m_ii = 1; ld_ii = p2_ii; s_ii = p3_ii;
                GPmLight = LightAmp_LST(la_ii);
                GPmLightDur = GPmLightDur_LST(ld_ii);
                GPm_w_mn = GPmVLw_mean_LST(m_ii);
                GPm_w_sg = GPmVLw_sig_LST(s_ii);
                
                PhotoStop = PhotoInjT+GPmLightDur;
                
                
                simTxt = ['I_A_m_p =' num2str(GPmLight) ', I_d_u_r =' num2str(GPmLightDur) ',GPmVL mean =' num2str(GPm_w_mn) ', sig =' num2str(GPm_w_sg)];
                spkBin = ACT_Record{p1_ii,p2_ii,p3_ii}.WT.All.spktrain;
                [freq_all_WT,spkTime_all_WT, fg_handle_WT ] = raster_from_spkbin_BurstRange( spkBin,PhotoStop, Tstop,PhotoInjT, PhotoStop, DelayT, BurstRange, [ simTxt ' : WT']);
                set(fg_handle_WT, 'position',[  449   450   791   528])
                spkBin = ACT_Record{p1_ii,p2_ii,p3_ii}.KO.All.spktrain;
                [freq_all_KO,spkTime_all_KO, fg_handle_KO ] = raster_from_spkbin_BurstRange( spkBin,PhotoStop, Tstop,PhotoInjT, PhotoStop, DelayT, BurstRange, [ simTxt ' : KO']);
                set(fg_handle_KO, 'position',[  449   450   791   528])
                if (SAVE_FIG)
                    ffig = [ dirLoc dirFig 'RasterPlot_' saveTxt1 num2str(PARAM1(p1_ii)) '_' saveTxt2 num2str(PARAM2(p2_ii)) '_' saveTxt3 num2str(PARAM3(p3_ii))];
                    saveas( fg_handle_WT, [ffig '_WT.jpg'], 'jpg')
                    saveas( fg_handle_WT, [ffig '_WT.fig'], 'fig')
                    saveas( fg_handle_KO, [ffig '_KO.jpg'], 'jpg')
                    saveas( fg_handle_KO, [ffig '_KO.fig'], 'fig')
                end
                CntBurstSpkWT(p1_ii,p2_ii,p3_ii) = mean(ACT_Record{p1_ii,p2_ii,p3_ii}.WT.All.BurstSpk);
                CntBurstSpkKO(p1_ii,p2_ii,p3_ii) = mean(ACT_Record{p1_ii,p2_ii,p3_ii}.KO.All.BurstSpk);
                ChanceOfBurstingWT(p1_ii,p2_ii,p3_ii) =  sum((ACT_Record{p1_ii,p2_ii,p3_ii}.WT.All.BurstSpk > 0))/(nE+nI);
                ChanceOfBurstingKO(p1_ii,p2_ii,p3_ii) =  sum((ACT_Record{p1_ii,p2_ii,p3_ii}.KO.All.BurstSpk > 0))/(nE+nI);
                
                cnt = cnt+1;
                figure(fLA1);
                subplot(length(PARAM2), length(PARAM3),cnt)
                [nb,xb]=hist(ACT_Record{p1_ii,p2_ii,p3_ii}.WT.All.fr_AfterLightOff,10); bh=bar(xb,nb); hold on;
                set(bh,'facecolor',[0 0 0]);
                xlim([0 max(xb)+1])
                title([ titleTxt2 ' = ' num2str(PARAM2(p2_ii)) ', ' titleTxt3 ' = ' num2str(PARAM3(p3_ii)) ])
                figure(fLA2);
                subplot(length(PARAM2), length(PARAM3),cnt)
                [nb,xb]=hist(ACT_Record{p1_ii,p2_ii,p3_ii}.KO.All.fr_AfterLightOff,10); bh=bar(xb,nb); hold on;
                set(bh,'facecolor',[1 0 0]);
                xlim([0 max(xb)+1])
                title([ titleTxt2 ' = ' num2str(PARAM2(p2_ii)) ', ' titleTxt3 ' = ' num2str(PARAM3(p3_ii)) ])
                
                figure(fB1);
                subplot(length(PARAM2), length(PARAM3),cnt)
                [nb,xb]=hist(ACT_Record{p1_ii,p2_ii,p3_ii}.WT.All.BurstSpk); bh=bar(xb,nb); hold on;
                set(bh,'facecolor',[0 0 0]);
                xlim([0 max(xb)+1])
                title([ titleTxt2 ' = ' num2str(PARAM2(p2_ii)) ', ' titleTxt3 ' = ' num2str(PARAM3(p3_ii)) ])
                figure(fB2);
                subplot(length(PARAM2), length(PARAM3),cnt)
                [nb,xb]=hist(ACT_Record{p1_ii,p2_ii,p3_ii}.KO.All.BurstSpk,10); bh=bar(xb,nb); hold on;
                set(bh,'facecolor',[1 0 0]);
                xlim([0 max(xb)+1])
                title([ titleTxt2 ' = ' num2str(PARAM2(p2_ii)) ', ' titleTxt3 ' = ' num2str(PARAM3(p3_ii)) ])
            end
        end
        figure(fLA1); suptitle(['Avg Fr after light off WT: GPmInput Amp =' num2str(GPmLight) ])
        figure(fB1); suptitle(['Number of Burst spike WT: GPmInput Amp =' num2str(GPmLight) ])
        figure(fLA2); suptitle(['Avg Fr after light off KO: GPmInput Amp =' num2str(GPmLight) ])
        figure(fB2);  suptitle(['Number of Burst spike KO: GPmInput Amp =' num2str(GPmLight) ])
        if (SAVE_FIG)
            ffig = [ dirLoc dirFig 'AvgFrLightOff_' saveTxt1 num2str(PARAM1(p1_ii))];
            saveas( fLA1, [ffig '_WT.jpg'], 'jpg')
            saveas( fLA1, [ffig '_WT.fig'], 'fig')
            saveas( fLA2, [ffig '_KO.jpg'], 'jpg')
            saveas( fLA2, [ffig '_KO.fig'], 'fig')
            ffig = [ dirLoc dirFig 'NumBurstSpk_' saveTxt1 num2str(PARAM1(p1_ii))];
            saveas( fB1, [ffig '_WT.jpg'], 'jpg')
            saveas( fB1, [ffig '_WT.fig'], 'fig')
            saveas( fB2, [ffig '_KO.jpg'], 'jpg')
            saveas( fB2, [ffig '_KO.fig'], 'fig')
        end
    end
    
    %
    for p1_ii = 1 : length(PARAM1)
        Bfg = figure;  set(Bfg,'position',[680   283   712   695]); set(gcf,'PaperPositionMode','auto');
        Bfg2 = figure;  set(Bfg2,'position',[680   283   712   695]);  set(gcf,'PaperPositionMode','auto');
        
        for p2_ii = 1 : length(PARAM2)
            figure(Bfg)
            subplot(3,2, p2_ii); hold on;
            plot( PARAM3, squeeze(CntBurstSpkWT(p1_ii,p2_ii,:)),'*-k');
            plot( PARAM3, squeeze(CntBurstSpkKO(p1_ii,p2_ii,:)),'*-r');
            legend('WT','KO','location','Best')
            title([ titleTxt2 ' = ' num2str(PARAM2(p2_ii))]);
            xlabel('Sigma'); ylabel('Avg Burst spike')
            
            figure(Bfg2)
            subplot(3,2, p2_ii); hold on;
            plot( PARAM3, squeeze(ChanceOfBurstingWT(p1_ii,p2_ii,:)),'*-k');
            plot( PARAM3, squeeze(ChanceOfBurstingKO(p1_ii,p2_ii,:)),'*-r');
            legend('WT','KO','location','Best')
            title([ titleTxt2 ' = ' num2str(PARAM2(p2_ii))]);
            xlabel('Sigma'); ylabel('Chance of Burst spike')
            
        end
        figure(Bfg)
        suptitle([ legTxt1 ' = ' num2str(PARAM1(p1_ii))]);
        figure(Bfg2)
        suptitle([ 'Chance of Bursting, ' legTxt1 ' = ' num2str(PARAM1(p1_ii))]);
        
        if (SAVE_FIG)
            ffig = [ dirLoc dirFig 'NumBurstSpk_' saveTxt1 num2str(PARAM1(p1_ii)) ];
            saveas(  Bfg, [ffig '.jpg'], 'jpg');        saveas(  Bfg, [ffig '.fig'], 'fig');
            ffig = [ dirLoc dirFig 'ChanceOfBurst_' saveTxt1 num2str(PARAM1(p1_ii)) ];
            saveas(  Bfg2, [ffig '.jpg'], 'jpg');        saveas(  Bfg2, [ffig '.fig'], 'fig');
        end
    end
    
end
