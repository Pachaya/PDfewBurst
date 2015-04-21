%% For compare WT&KO response of one set of parameters

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
FNAME_SIM = 'PD4' ; %'PreLim_INOISE0.2_InSPK50_EE1_EI2_IE20_II12_mult125_w0.0006';
plotSampleV = 0;
SPECIFIED_RATIO = 0; %TestSim_WT_InGauss0.2_W0.0005_50.00Hz_rEE0_rEI1_rIE2_rII1_Wmult10
Basal_Act_set  = [];

NoiseMEAN = 0;

SAVE_FIG = 1;

dirLoc = 'TestCenterM1_TC_optParam/' ; %'TestSim_BasalAct/';  %'VL_LocalConn/';
rEE_rEI_rIE_rII = [1 1 1 1];
w_MULT = 125;
avgFR_CUTTIME = 500;

pulseHz = 1; N_pulse = 1;
Input_I_amp_list = [0];
Input_I_dur_list = [0];
Input_I_amp = 0; %[-1; -2; -3; -4; -5;];% [-1; -5]; %
Input_I_dur = 500; %[10; 20; 30; 40; 50;];
ReboundPeakAmp_WTlist = zeros(length(Input_I_amp_list),length(Input_I_dur_list));
ReboundPeakDel_WTlist = zeros(length(Input_I_amp_list),length(Input_I_dur_list));
ReboundPeakAmp_KOlist = zeros(length(Input_I_amp_list),length(Input_I_dur_list));
ReboundPeakDel_KOlist = zeros(length(Input_I_amp_list),length(Input_I_dur_list));
Save_BasalAct =  cell(length(Input_I_amp_list),length(Input_I_dur_list));
InjStartT = 1000;
InjStopT = InjStartT+Input_I_dur;
reboundDelay = 0;
burstRange = 100;
cutTime = avgFR_CUTTIME;

ADD_I_to_M1 = 0;
ADD_I_to_VL = 0;

FIG_ALL = 0;
dirFig = 'Fig_150303NoiseOnly/';
NoiseMEAN_WTKO = [0; 0;]; %[0.0291; 0.0075;]; %[0.0295; 0.016]; %%%%%%%%%%%%%%%%%%%%%%%%%%% [ WT; KO;] for 10 Hz -> [0.0295; 0.016]; for 5 Hz -> [0.0291; 0.0075;]; 5Hz 0.0318 0.3  for WT 5Hz (2) -> [0.0318; 0.0075;]
NoiseSIG_WTKO  =  [0; 0;]; %[0.8; 0.8;];  %[0.24; 0.3];  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  [WT; KO;] for 10 Hz -> [0.4; 0.5]; for 5 Hz -> [0.24; 0.3];  5Hz 0.0318 0.3  for WT 5Hz (2) -> [0.3; 0.3;]
VARIE_MEAN_InGauss = 1;
NoiseMEANsigma_WTKO = [0; 0;]; %%%%% [WT; KO;]
TSTOP = 3000; 

rTC_LST = [120:10:150];
wmTC_LST = [10;];

PARAM1 = rTC_LST;
legTxt1 = 'TC range'; %'E_L '; %'soma size';
saveTxt1 = 'rTC';
PARAM2 = wmTC_LST;
legTxt2 = 'TC weight';
saveTxt2 = 'wmTC';


 for p1_ii = 1 : length(PARAM1)
     for p2_ii = 1 : length(PARAM2)
        rTC_ii = p1_ii; wmTC_ii = p2_ii;
for cell_type = 1 : 2
    NoiseMEAN = NoiseMEAN_WTKO(cell_type);  %%%%%%%%%%%%%%%%%%%%%%%%%%% [ WT; KO;]
    NoiseSIG  = NoiseSIG_WTKO(cell_type);  %%%%%%%%%%%%%%[WT; KO;]
    NoiseMEANsigma = NoiseMEANsigma_WTKO(cell_type);
    
    InjStopT = InjStartT+Input_I_dur;
    coreFileName = 'PD4_testM1_TC_' ; %PD4_testM1_TC_1Hz_amp-0.3_dur500_N1_rTC150_wmTC6_KO_IGmean0_IGmeanSig0_W0.001_T3000
    coreFileName = [coreFileName num2str(pulseHz)  'Hz_amp' num2str(Input_I_amp) '_dur' num2str(Input_I_dur) '_N' num2str(N_pulse) '_'];
    
    rTC = rTC_LST(rTC_ii); wmTC = wmTC_LST(wmTC_ii);
    coreFileName = [coreFileName 'rTC' num2str(rTC) '_' ];
    coreFileName = [coreFileName 'wmTC' num2str(wmTC) '_' ];
    
    InGauss_STDEV = NoiseSIG; %0.2;
    W_Weight = 0.001;
    PoisInputFr = 40;
    
    if (cell_type == 1)
        coreFileName = [coreFileName 'WT_'];  % TestInGauss_WT_W0.0005_40.00Hz
    elseif (cell_type == 2)
        coreFileName = [coreFileName 'KO_'];
    end
    
    
    if(InGauss_STDEV == 0)
        if(VARIE_MEAN_InGauss)
            Name_postfix = [coreFileName 'IGmean' num2str(NoiseMEAN) '_IGmeanSig' num2str(NoiseMEANsigma) '_W' num2str(W_Weight) ];
        else
            Name_postfix = [coreFileName 'W' num2str(W_Weight) ];
        end
    else
        if(VARIE_MEAN_InGauss)
            Name_postfix = [coreFileName 'InGauss' num2str(InGauss_STDEV) '_IGmean' num2str(NoiseMEAN) '_IGmeanSig' num2str(NoiseMEANsigma) '_W' num2str(W_Weight) ];
        else
            Name_postfix = [coreFileName 'InGauss' num2str(InGauss_STDEV) '_W' num2str(W_Weight) ];
        end
    end
    if (PoisInputFr ~= 0)
    Name_postfix = [Name_postfix  '_' num2str(PoisInputFr) '.00Hz' ];
    end
    if(SPECIFIED_RATIO)
        rEE = num2str(rEE_rEI_rIE_rII(1));    rEI = num2str(rEE_rEI_rIE_rII(2));
        rIE = num2str(rEE_rEI_rIE_rII(3));    rII = num2str(rEE_rEI_rIE_rII(4));
        WMult = num2str(w_MULT);
        tmpTxt = Name_postfix;
        Name_postfix = [tmpTxt  '_rEE' rEE '_rEI' rEI '_rIE' rIE '_rII' rII '_Wmult' WMult ]; % For first run that did not add InSpk = 50 Hz in Name
    end
    Name_postfix = [Name_postfix '_T' num2str(TSTOP)];
%     PD4_testM1_TC_1Hz_amp-0.3_dur500_N1_rTC150_wmTC7_KO_IGmean0_IGmeanSig0_W0.001_T3000
    
    disp('==================================================================================================')
    disp(Name_postfix)
    if (cell_type == 1)
        cTxt = 'WT';
    else
        cTxt = 'KO';
    end
    figNameCode = [ saveTxt1 num2str(PARAM1(p1_ii))  '_' saveTxt2 num2str(PARAM2(p2_ii)) '_' cTxt];
    %             VL_Basal_ActivityCenter %call function
    %             VL_Basal_Activity
    getCenterPart_M1 = 0; 
    Sim_Neuron_Activity_PD4  
   
   
    
    tmp.E = E;      tmp.I = I;  tmp.All = All; tmp.VL = VL; tmp.M1 = M1;
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
        title([cTxt ': sample membrane potential, Inject I : amp = ' num2str(Input_I_amp) ', dur = ' num2str(Input_I_dur) ])
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
DoTTest = 0;
if (DoTTest)
    % two-tailed t-test % during normal 
    
    % disp('PoisSPk')
    % [hP,pP] = ttest(WT.PoisSpk.fr_data, KO.PoisSpk.fr_data)
    disp('--------------------------------------------------------------------------------------------------')
    disp('E cells : two-tailed t-test')
    [hE,pE] = ttest2(WT.VL.normal.fr_data, KO.VL.normal.fr_data, 'Tail','both');
    disp([ ' p-value = ' num2str(pE)])
    if(hE)
        disp('H0 rejected : average normal firing rate of E cells in WT significantly different from KO (p < 0.05)')
    else
        disp('Do not rejected H0 : average normal firing rate of E cells in WT does not significantly different from KO (p > 0.05)')
    end
     Basal_Act.ttest2.hE = hE;
    Basal_Act.ttest2.pE = pE;
       if(ADD_I_to_VL)
    disp('--------------------------------------------------------------------------------------------------')
    disp('I cells : two-tailed t-test')
    [hI,pI] = ttest2(WT.I.fr_data, KO.I.fr_data);
    disp([ ' p-value = ' num2str(pI)])
    if(hI)
        disp('H0 rejected : average normal firing rate of I cells in WT significantly different from KO (p < 0.05)')
    else
        disp('Do not rejected H0 : average normal firing rate of I cells in WT does not significantly different from KO (p > 0.05)')
    end
    disp('--------------------------------------------------------------------------------------------------')
    disp('All cells : two-tailed t-test')
    [hA,pA]  = ttest2([WT.VL.normal.fr_data; WT.I.fr_data], [KO.VL.normal.fr_data; KO.I.fr_data]);
    disp([ ' p-value = ' num2str(pA)])
    if(hA)
        disp('H0 rejected : average normal firing rate of all cells in WT significantly different from KO (p < 0.05)')
    else
        disp('Do not rejected  H0 : average normal firing rate of all cells in WT does not significantly different from KO (p > 0.05)')
    end
        Basal_Act.ttest2.hI = hI;
    Basal_Act.ttest2.pI = pI;
    Basal_Act.ttest2.hA = hA;
    Basal_Act.ttest2.pA = pA;
       end

end
%% Sample T-Test 
DoSampleTtest = 1; 
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
    [ht,pt] = ttest2(WT.VL.normal.fr_data(sample), KO.VL.normal.fr_data(sample));
    [h,p] = ranksum(WT.VL.normal.fr_data(sample), KO.VL.normal.fr_data(sample));
    
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

     end
 end

%% Raster Plot
% [freq_all_WT,spkTime_all_WT, fg_handle_WT ] = raster_from_spkbin( [WT.E.spktrain; WT.I.spktrain],avgFR_CUTTIME, Tstop, 'center cells : WT');
% [freq_all_KO,spkTime_all_KO, fg_handle_KO ] = raster_from_spkbin( [KO.E.spktrain; KO.I.spktrain],avgFR_CUTTIME, Tstop, 'center cells : KO');
if ( 0) 
%%  Testing Filter Size    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TestFilterSize = 0; 
if(TestFilterSize)
DOWNLOAD_SAVE_BASALACT = 0;
if(DOWNLOAD_SAVE_BASALACT)
    load([dirLoc 'Data_amp-1_-5_dur10_50_Trial1fixed.mat'])
end
if(0)
    for ii = 1 : 10
        figure(ii)
        k = waitforbuttonpress;
    end
end

SIGMA_fr_All_list =[0 1 5:5:21];
zoomx1 = 900; zoomx2 = 1450;
fg_SIGMA = figure; set(fg_SIGMA, 'position',[445    48   993   948]);
subplot(length(SIGMA_fr_All_list)+4,1,1:2);
raster_from_spkbin_cga( [WT.E.spktrain; WT.I.spktrain],avgFR_CUTTIME, Tstop,2);
ylim([80 160]); xlim([zoomx1 zoomx2]); title('WT'); ylabel('cell ID#')
subplot(length(SIGMA_fr_All_list)+4,1,3:4);
raster_from_spkbin_cga( [KO.E.spktrain; KO.I.spktrain],avgFR_CUTTIME, Tstop,2);
ylim([80 160]); xlim([zoomx1 zoomx2]); title('KO'); ylabel('cell ID#')

for S_ii = 1 : length(SIGMA_fr_All_list)
    ncells = NcenterE+NcenterI;
    SIGMA_fr_All = SIGMA_fr_All_list(S_ii); %min(Input_I_dur,40);
    
    if( SIGMA_fr_All == 0)
        sWT_All = sum(WT.All.spktrain(1:end,:));
        sKO_All = sum(KO.All.spktrain(1:end,:));
    else
        %All
        sigma =  SIGMA_fr_All; % ms
        edges =[-3*sigma:1:3*sigma];
        kernel = normpdf(edges,0,sigma);
        kernel = kernel*RES; % Normalized with bin width
        figure; plot(kernel); title(['sigma = ' num2str(sigma)])
        sWT_All = conv(sum(WT.All.spktrain(1:end,:)),kernel);
        sKO_All = conv(sum(KO.All.spktrain(1:end,:)),kernel);
        center = ceil(length(edges)/2);
        sWT_All = sWT_All(center:Tstop+center-1);
        sKO_All = sKO_All(center:Tstop+center-1);
    end
    %     fg_Allcut = figure;  set(fg_Allcut, 'position', [  221         443        1606         475]);
    figure(fg_SIGMA) ; subplot(length(SIGMA_fr_All_list)+4,1,S_ii+4)
    avgFR_WT_All = sWT_All/ncells*1000;
    avgFR_KO_All = sKO_All/ncells*1000;
    plot(avgFR_WT_All,'k','LineWidth',2);  hold on;
    plot(avgFR_KO_All,'m','LineWidth',2);  hold on;
    %     legend('WT','KO')
    % title(['All cell, Inject I : amp = ' num2str(Input_I_amp) ', dur = ' num2str(Input_I_dur) ])
    %     ylabel('Average firing rate (Hz)')
    %     xlabel('Time(ms)')
    
    xlim([zoomx1 zoomx2])
    ylim( [0 max([ avgFR_WT_All avgFR_WT_All])+2])
    %     legend('WT','KO')
    
    % PeakAvgFr_WT_All = max(avgFR_WT_All(zoomx1:zoomx2));
    %     [pks, locs] = findpeaks(avgFR_WT_All); pks =  pks(locs>InjStopT); locs = locs(locs > InjStopT); locs = locs(1); pks =  pks(1); % Find the first peak
    %     peakvalWT_All = pks; peakTimeWT_All = locs;
    %     [pks, locs] = findpeaks(avgFR_KO_All);  pks =  pks(locs>InjStopT); locs = locs(locs > InjStopT); locs = locs(1); pks =  pks(1); % Find the first peak
    %     peakvalKO_All = pks; peakTimeKO_All = locs;
    %     % [peakvalWT_All, peakTimeWT_All] = max(avgFR_WT_All(InjStopT:InjStopT+burstRange)); % Consider finding after stop injection negative current only
    %     % [peakvalKO_All, peakTimeKO_All] = max(avgFR_KO_All(InjStopT:InjStopT+burstRange));
    %     % peakTimeWT_All = peakTimeWT_All + InjStopT ;
    %     % peakTimeKO_All = peakTimeKO_All + InjStopT ;
    %     % figure(17); hold on;
    %     scatter(peakTimeWT_All, peakvalWT_All,150,'*k','LineWidth',2); peakTimeWT_All = peakTimeWT_All-1; hold on; %plot start at time 0
    %     scatter(peakTimeKO_All, peakvalKO_All,150,'*m','LineWidth',2); peakTimeKO_All = peakTimeKO_All -1;
    rectangle('Position',[InjStartT , 0,Input_I_dur,  max([ avgFR_WT_All avgFR_WT_All])+2],'EdgeColor','r');
    set(gca,'box','off')
    title([ 'Sigma = ' num2str( SIGMA_fr_All)])
    if(S_ii == round(length(SIGMA_fr_All_list)/2))
        ylabel('Average firing rate (Hz)')
    end
    %     title({[ 'Sigma = ' num2str( SIGMA_fr_All) ' : Average firing rate Fr(t) cell in WT and KO , Inject I : amp = ' num2str(Input_I_amp) ', dur = ' num2str(Input_I_dur)],['WT rebound peak = ' num2str(peakvalWT_All) ' Hz at ' num2str(peakTimeWT_All) ' ms , KO rebound peak = ' num2str(peakvalKO_All) ' Hz at ' num2str(peakTimeKO_All) ' ms']});
    %%% Raster Plot and Fr(t)
    %     raster_from_spkbin2(  WT.E.spktrain, WT.I.spktrain,avgFR_CUTTIME, Tstop,[ 'WT, Inject I : amp = ' num2str(Input_I_amp) ', dur = ' num2str(Input_I_dur) ],RES, SIGMA_fr_E,SIGMA_fr_I,SIGMA_fr_All)
    %     raster_from_spkbin2(  KO.E.spktrain, KO.I.spktrain,avgFR_CUTTIME, Tstop, ['KO, Inject I : amp = ' num2str(Input_I_amp) ', dur = ' num2str(Input_I_dur) ],RES, SIGMA_fr_E,SIGMA_fr_I,SIGMA_fr_All)
end
figure(fg_SIGMA)

xlabel('Time(ms)')
suptitle(['Inject I : amp = ' num2str(Input_I_amp) ', dur = ' num2str(Input_I_dur)   ])
end
%% Average firing rate
ncells = size(WT.E.spktrain,1) +size(WT.I.spktrain,1);
SIGMA_fr_E = 20;
SIGMA_fr_I = 20;
SIGMA_fr_All = 15; %min(Input_I_dur,40);
%E
sigma =  SIGMA_fr_E; % ms
edges =[-3*sigma:1:3*sigma];
kernel = normpdf(edges,0,sigma);
kernel = kernel*RES; % Normalized with bin width
sWT_E = conv(sum(WT.E.spktrain(1:end,:)),kernel);
sKO_E = conv(sum(KO.E.spktrain(1:end,:)),kernel); length( conv(sum(WT.E.spktrain(1:end,:)),kernel))
center = ceil(length(edges)/2);
sWT_E = sWT_E(center:Tstop+center); %Note time in  spike train start from 0
sKO_E = sKO_E(center:Tstop+center);

% figure; plot(kernel); %just test delete later
% yy = func_gauss(edges, sigma); figure; plot(yy*RES) %----> check peak
% %% fft of overall fr
% T = Tstop; % size(WT.E.spktrain,2);
% fftWT_E = abs(fftshift(fft(sWT_E)));
% fftKO_E = abs(fftshift(fft(sKO_E)));
% dHz = 1/T * 1000;
% Hz_ind = [-round(T/2):round(T/2)]*dHz;
% Hz_ind = Hz_ind(1:T);
%
%     %smooth  FFT
%     gsig = 14;   % filter [Hz] ---------------------------------------> Calculate percent different when changing sigma size
%     gfil2 = exp( -[-round(4*gsig):round(4*gsig)].^2 /2/gsig^2);
%     gfil2 = gfil2/sum(gfil2(:));
%
%     fft_smoothWT_E = conv(fftWT_E, gfil2, 'same');
%     fft_smoothKO_E = conv(fftKO_E, gfil2, 'same');
%     fftLimX1 = 1;   fftLimX2 = 10;
%     figure;
%      hold on
%     plot(Hz_ind, fftWT_E, 'k-')
%     plot(Hz_ind, fftKO_E, 'm-')
%     hold on
%     xlim([fftLimX1 fftLimX2])
%     xlim([1 5])

% %%
PLOT_separate = 0;
if(PLOT_separate)
    fFR_E = figure; set(fFR_E, 'position', [578   185   780   584]); subplot(211);
    plot(sWT_E/NcenterE,'m');  hold on
    plot(sKO_E/NcenterE,'k');  hold on
    legend('WT','KO')
    title(['Average spike density functions of E cell in WT and KO']);
    ylabel('Spike per ms')
    xlabel('Time(ms)')
    subplot(212)
    plot(sWT_E/NcenterE*1000,'m');  hold on;
    plot(sKO_E/NcenterE*1000,'k');  hold on;
    legend('WT','KO')
    title('Average firing rate Fr(t) of E cell in WT and KO')
    ylabel('Average firing rate (Hz)')
    xlabel('Time(ms)')
end

%I
sigma =  SIGMA_fr_I; % ms
edges =[-3*sigma:1:3*sigma];
kernel = normpdf(edges,0,sigma);
kernel = kernel*RES; % Normalized with bin width
sWT_I = conv(sum(WT.I.spktrain(1:end,:)),kernel);
sKO_I = conv(sum(KO.I.spktrain(1:end,:)),kernel);
center = ceil(length(edges)/2);
sWT_I = sWT_I(center:Tstop+center-1);
sKO_I = sKO_I(center:Tstop+center-1);
if(PLOT_separate)
    fFR_I = figure; set(fFR_I, 'position', [578   185   780   584]);  subplot(211)
    plot(sWT_I/NcenterI,'k');  hold on
    plot(sKO_I/NcenterI,'m');  hold on
    legend('WT','KO')
    title(['Average spike density functions of I cell in WT and KO']);
    ylabel('Spike per ms')
    xlabel('Time(ms)')
    subplot(212)
    plot(sWT_I/NcenterI*1000,'k');  hold on;
    plot(sKO_I/NcenterI*1000,'m');  hold on;
    legend('WT','KO')
    title('Average firing rate Fr(t) of I cell in WT and KO')
    ylabel('Average firing rate (Hz)')
    xlabel('Time(ms)')
end


figure;
plot(sum(WT.All.spktrain(1:end,:)),'k');      hold on;
plot(sum(KO.All.spktrain(1:end,:)),'m');

%All
sigma =  SIGMA_fr_All; % ms
edges =[-3*sigma:1:3*sigma];
kernel = normpdf(edges,0,sigma);
kernel = kernel*RES; % Normalized with bin width
% figure; plot(kernel)
sWT_All = conv(sum(WT.All.spktrain(1:end,:)),kernel);
sKO_All = conv(sum(KO.All.spktrain(1:end,:)),kernel);
center = ceil(length(edges)/2);
sWT_All = sWT_All(center:Tstop+center-1);
sKO_All = sKO_All(center:Tstop+center-1);
if(PLOT_separate)
    fFR_All = figure; set(fFR_All, 'position', [578   185   780   584]); subplot(211)
    plot(sWT_All/ncells,'k');  hold on
    plot(sKO_All/ncells,'m');  hold on
    legend('WT','KO')
    title(['Average spike density functions of All cell in WT and KO']);
    ylabel('Spike per ms')
    xlabel('Time(ms)')
    subplot(212)
    plot(sWT_All/ncells*1000,'k');  hold on;
    plot(sKO_All/ncells*1000,'m');  hold on;
    legend('WT','KO')
    title('Average firing rate Fr(t) of All cell in WT and KO')
    ylabel('Average firing rate (Hz)')
    xlabel('Time(ms)')
end

% Combine
if(PLOT_separate)
    fFR_Combine = figure; set(fFR_Combine, 'position', [  589   230   780   767]);  subplot(311)
    plot(sWT_E/NcenterE*1000,'k');  hold on;
    plot(sKO_E/NcenterE*1000,'m');  hold on;
    legend('WT','KO')
    title('E cell')
    ylabel('Average firing rate (Hz)')
    xlabel('Time(ms)')
    subplot(312)
    plot(sWT_I/NcenterI*1000,'k');  hold on;
    plot(sKO_I/NcenterI*1000,'m');  hold on;
    legend('WT','KO')
    title('I cell')
    ylabel('Average firing rate (Hz)')
    xlabel('Time(ms)')
    subplot(313)
    plot(sWT_All/ncells*1000,'k');  hold on;
    plot(sKO_All/ncells*1000,'m');  hold on;
    legend('WT','KO')
    title('All cell')
    ylabel('Average firing rate (Hz)')
    xlabel('Time(ms)')
    suptitle(['Average firing rate Fr(t) cell in WT and KO , Inject I : amp = ' num2str(Input_I_amp) ', dur = ' num2str(Input_I_dur) ])
end
%%

fg_Allcut = figure;  set(fg_Allcut, 'position', [  221         443        1606         475]);
avgFR_WT_All = sWT_All/ncells*1000;
avgFR_KO_All = sKO_All/ncells*1000;
plot(avgFR_WT_All,'k','LineWidth',3);  hold on;
plot(avgFR_KO_All,'m','LineWidth',3);  hold on;
legend('WT','KO')
% title(['All cell, Inject I : amp = ' num2str(Input_I_amp) ', dur = ' num2str(Input_I_dur) ])
ylabel('Average firing rate (Hz)')
xlabel('Time(ms)')
zoomx1 = 100; zoomx2 = 2000;
xlim([zoomx1 zoomx2])
ylim( [0 max([ avgFR_WT_All avgFR_WT_All])+2])
legend('WT','KO')

% PeakAvgFr_WT_All = max(avgFR_WT_All(zoomx1:zoomx2));
[pks, locs] = findpeaks(avgFR_WT_All); pks =  pks(locs>InjStopT); locs = locs(locs > InjStopT); locs = locs(1); pks =  pks(1); % Find the first peak
peakvalWT_All = pks; peakTimeWT_All = locs;
[pks, locs] = findpeaks(avgFR_KO_All);  pks =  pks(locs>InjStopT); locs = locs(locs > InjStopT); locs = locs(1); pks =  pks(1); % Find the first peak
peakvalKO_All = pks; peakTimeKO_All = locs;
% [peakvalWT_All, peakTimeWT_All] = max(avgFR_WT_All(InjStopT:InjStopT+burstRange)); % Consider finding after stop injection negative current only
% [peakvalKO_All, peakTimeKO_All] = max(avgFR_KO_All(InjStopT:InjStopT+burstRange));
% peakTimeWT_All = peakTimeWT_All + InjStopT ;
% peakTimeKO_All = peakTimeKO_All + InjStopT ;
% figure(17); hold on;
scatter(peakTimeWT_All, peakvalWT_All,150,'*k','LineWidth',2); peakTimeWT_All = peakTimeWT_All-1; hold on; %plot start at time 0
scatter(peakTimeKO_All, peakvalKO_All,150,'*m','LineWidth',2); peakTimeKO_All = peakTimeKO_All -1;
rectangle('Position',[InjStartT , 0,Input_I_dur,  max([ avgFR_WT_All avgFR_WT_All])+2],'EdgeColor','r');
set(gca,'box','off')
title({['Average firing rate Fr(t) cell in WT and KO , Inject I : amp = ' num2str(Input_I_amp) ', dur = ' num2str(Input_I_dur)],['WT rebound peak = ' num2str(peakvalWT_All) ' Hz at ' num2str(peakTimeWT_All) ' ms , KO rebound peak = ' num2str(peakvalKO_All) ' Hz at ' num2str(peakTimeKO_All) ' ms']});
%% Raster Plot and Fr(t)
raster_from_spkbin2(  WT.E.spktrain, WT.I.spktrain,avgFR_CUTTIME, Tstop,[ 'WT, Inject I : amp = ' num2str(Input_I_amp) ', dur = ' num2str(Input_I_dur) ],RES, SIGMA_fr_E,SIGMA_fr_I,SIGMA_fr_All)
raster_from_spkbin2(  KO.E.spktrain, KO.I.spktrain,avgFR_CUTTIME, Tstop, ['KO, Inject I : amp = ' num2str(Input_I_amp) ', dur = ' num2str(Input_I_dur) ],RES, SIGMA_fr_E,SIGMA_fr_I,SIGMA_fr_All)
%% save
% ReboundPeakAmp_WTlist(A_ii, D_ii) = peakvalWT_All ; ReboundPeakDel_WTlist(A_ii, D_ii) = peakTimeWT_All;
% ReboundPeakAmp_KOlist(A_ii, D_ii) = peakvalKO_All; ReboundPeakDel_KOlist(A_ii, D_ii) = peakTimeKO_All;
Basal_Act.WT.avgFR = avgFR_WT_All;
Basal_Act.KO.avgFR = avgFR_KO_All;
A_ii = 1; D_ii = 1; 
Save_BasalAct{A_ii, D_ii}  = Basal_Act;
% disp(['( ' num2str(A_ii) ',' num2str(D_ii) ' ) I amp = ' num2str(Input_I_amp_list(A_ii)) ', dur = ' num2str(Input_I_dur_list(D_ii))  ])
disp(['======================  Inject I : amp = ' num2str(Input_I_amp) ', dur = ' num2str(Input_I_dur) '  ======================']);
%     end
% end
%%
%% Testing  Multiple parameter set --> fixed sigma to be 0, 1, 6 , 11
if(0)
SIGMA_List = [ 0 1 5 10]; CASE = 2;
 fftLimX1 = 1;   fftLimX2 = 10;
if(CASE ==1)
    Injection_Param = [ -1 50; -3 50; -5 50;];  Injection_Param_ID = [1 5; 3 5; 5 5;];
end
if(CASE == 2)
    Injection_Param = [ -1 10; -1 20; -1 50;]; Injection_Param_ID = [1 1; 1 2; 1 5;];
end
if(0)
    for ii = 1 : 37
        figure(ii)
        k = waitforbuttonpress;
    end
end


% zoomx1 = 450; zoomx2 = 1000;

fg_SIGMA = figure; set(fg_SIGMA, 'position',[  274          64        1404         915]);
fg_fft = figure; set(fg_fft, 'position',[  274          64        1404         915]);

% figure;
% NRows = size(Injection_Param_ID,1)+4; NCols = length(SIGMA_List);

% for cc = 1 :  length(SIGMA_List)
%
%     WT = Save_BasalAct{Injection_Param_ID(1,1),Injection_Param_ID(1,2)}.WT;
%     KO = Save_BasalAct{Injection_Param_ID(1,1),Injection_Param_ID(1,2)}.KO;
%
% subplot(NRows, NCols ,[cc NCols+cc]);
% raster_from_spkbin_cga( [WT.E.spktrain; WT.I.spktrain],avgFR_CUTTIME, Tstop,2);
% ylim([80 160]); xlim([zoomx1 zoomx2]); title('WT'); ylabel('cell ID#')
% subplot(NRows, NCols ,[2*NCols+cc 3*NCols+cc]);
% raster_from_spkbin_cga( [KO.E.spktrain; KO.I.spktrain],avgFR_CUTTIME, Tstop,2);
% ylim([80 160]); xlim([zoomx1 zoomx2]); title('KO'); ylabel('cell ID#')
% end

NCols = size(Injection_Param_ID,1); NRows = length(SIGMA_List)+4;
for p_ii = 1: size(Injection_Param_ID ,1)
    cc = p_ii;
    WT = Save_BasalAct{Injection_Param_ID(p_ii,1),Injection_Param_ID(p_ii,2)}.WT;
    KO = Save_BasalAct{Injection_Param_ID(p_ii,1),Injection_Param_ID(p_ii,2)}.KO;
    figure( fg_SIGMA );
    subplot(NRows, NCols ,[cc NCols+cc]);
    raster_from_spkbin_cga( [WT.E.spktrain; WT.I.spktrain],avgFR_CUTTIME, Tstop,2);
    ylim([80 160]); xlim([zoomx1 zoomx2]); ylabel('WT cell ID#');
    ttt= title(['Inject I : amp = ' num2str( Injection_Param(p_ii,1)) ', dur = ' num2str( Injection_Param(p_ii,2))   ]);
    set(ttt,'FontSize', 14)
    
    subplot(NRows, NCols ,[2*NCols+cc 3*NCols+cc]);
    raster_from_spkbin_cga( [KO.E.spktrain; KO.I.spktrain],avgFR_CUTTIME, Tstop,2);
    ylim([80 160]); xlim([zoomx1 zoomx2]);  ylabel('KO cell ID#')
    
    for S_ii = 1 : length(SIGMA_List)
        ncells = NcenterE+NcenterI;
        SIGMA_fr_All = SIGMA_List(S_ii); %min(Input_I_dur,40);
        
        if( SIGMA_fr_All == 0)
            sWT_All = sum(WT.All.spktrain(1:end,:));
            sKO_All = sum(KO.All.spktrain(1:end,:));
        else
            %All
            sigma =  SIGMA_fr_All; % ms
            edges =[-3*sigma:1:3*sigma];
            kernel = normpdf(edges,0,sigma);
            kernel = kernel*RES; % Normalized with bin width
            %                     figure; plot(kernel); title(['sigma = ' num2str(sigma)])
            sWT_All = conv(sum(WT.All.spktrain(1:end,:)),kernel);
            sKO_All = conv(sum(KO.All.spktrain(1:end,:)),kernel);
            center = ceil(length(edges)/2);
            sWT_All = sWT_All(center:Tstop+center-1);
            sKO_All = sKO_All(center:Tstop+center-1);
        end
        %     fg_Allcut = figure;  set(fg_Allcut, 'position', [  221         443        1606         475]);
        figure(fg_SIGMA) ;
        subplot(NRows, NCols, (3+S_ii)*NCols +cc)
        avgFR_WT_All = sWT_All/ncells*1000;
        avgFR_KO_All = sKO_All/ncells*1000;
        plot(avgFR_WT_All,'k','LineWidth',2);  hold on;
        plot(avgFR_KO_All,'m','LineWidth',2);  hold on;
        
        xlim([zoomx1 zoomx2])
        ylim( [0 max([ avgFR_WT_All avgFR_WT_All])+2])
        
        rectangle('Position',[InjStartT , 0,Input_I_dur,  max([ avgFR_WT_All avgFR_WT_All])+2],'EdgeColor','r');
        set(gca,'box','off')
        title([ 'sigma = ' num2str( SIGMA_fr_All)]);
        
        if(S_ii == round(length(SIGMA_List)/2))
            figure(fg_SIGMA) ;
            yl=ylabel('Average firing rate (Hz)');
            %                     set(yl,'FontSize', 12)
        end
        % fft
        ttl  =  ['Sigma= ' num2str(SIGMA_fr_All) ' :Inject I : amp = ' num2str( Injection_Param(p_ii,1)) ', dur = ' num2str( Injection_Param(p_ii,2))   ] ;
        [Return_fft] = fft_from_sAvgFR (sWT_All, sKO_All, 5, 50,  1.5,ttl) ;
        figure(fg_fft);
        subplot(length(SIGMA_List), NCols, (S_ii-1)*NCols+cc);              hold on;
        plot(Return_fft.Hz_ind, Return_fft.fft_smoothWT_All_norm, 'k-');
        plot(Return_fft.Hz_ind, Return_fft.fft_smoothKO_All_norm, 'm-');
        cut_Hz_ind = (Return_fft.Hz_ind > 5) & (Return_fft.Hz_ind <50);
        [pks, locs] = max(Return_fft.fft_smoothWT_All_norm(cut_Hz_ind));
        tmpCut_Hz = Return_fft.Hz_ind(cut_Hz_ind); %tmpCut_fft = Return_fft.fft_smoothWT_All_norm(cut_Hz_ind);
        scatter(tmpCut_Hz(locs), pks,'*b');
        ttl_WT = ['WT: ' num2str(tmpCut_Hz(locs)) ' Hz, power = ' num2str(pks) ];
        [pks, locs] = max(Return_fft.fft_smoothKO_All_norm(cut_Hz_ind));
        tmpCut_Hz = Return_fft.Hz_ind(cut_Hz_ind);% tmpCut_fft = Return_fft.fft_smoothKO_All_norm(cut_Hz_ind);
        scatter(tmpCut_Hz(locs), pks,'*r');
        ttl_KO = [' KO: ' num2str(tmpCut_Hz(locs)) ' Hz, power = ' num2str(pks) ];
        xlim([fftLimX1 fftLimX2]); xlabel('Frequency (Hz)'); ylabel('Power (A.U.)')
        if(S_ii == 1)
            title({ttl,ttl_WT, ttl_KO})
        else
            title([ttl_WT ' ' ttl_KO ])
        end
        
    end
    if(p_ii == round(size(Injection_Param_ID ,1)/2))
        figure(fg_SIGMA) ;
        xl=xlabel('Time(ms)');
        set(xl,'FontSize', 12)
    end
end

if(CASE ==1)
    figure(fg_SIGMA)
    suptitle(['Different parameter sets when fixed duration to equal 50 ms']);
    figure(fg_fft);  suptitle(['Different parameter sets when fixed duration to equal 50 ms']);
end
if(CASE ==2)
    figure(fg_SIGMA)
    suptitle(['Different parameter sets when fixed amplitude to equal -1 nA']);
    figure(fg_fft); suptitle(['Different parameter sets when fixed amplitude to equal -1 nA']);
end

end
%%

%% fft of overall fr
fftLimX1 = 1; fftLimX2 = 50;
T = Tstop; % size(WT.E.spktrain,2);
figure; plot(sWT_All,'k'); hold on; plot(sKO_All,'m')
fftWT_All = abs(fftshift(fft(sWT_All)));
fftKO_All = abs(fftshift(fft(sKO_All)));
dHz = 1/T * 1000;
Hz_ind = [-round(T/2):round(T/2)]*dHz;
Hz_ind = Hz_ind(1:T);
figure;
plot(Hz_ind, fftWT_All, 'k-'); hold on;
plot(Hz_ind, fftKO_All, 'm-')
xlim([fftLimX1 fftLimX2])
%             xlim([1 5])
%smooth  FFT
gsig = 1.5;   % filter [Hz] ---------------------------------------> Calculate percent different when changing sigma size
gfil2 = exp( -[-round(4*gsig):round(4*gsig)].^2 /2/gsig^2);
gfil2 = gfil2/sum(gfil2(:));

fft_smoothWT_All = conv(fftWT_All, gfil2, 'same');
fft_smoothKO_All = conv(fftKO_All, gfil2, 'same');

figure;
hold on
plot(Hz_ind, fft_smoothWT_All, 'k-')
plot(Hz_ind, fft_smoothKO_All, 'm-')
hold on
xlim([fftLimX1 fftLimX2])


%%
if(0)
% peak percent different  --> KO = 1
Color_Letter = 'brkmcy';
figure;
Leg = cell( length(Input_I_dur_list),1);
for D_ii = 1 : length(Input_I_dur_list)
    plot(Input_I_amp_list, (ReboundPeakAmp_WTlist(:,D_ii)./ReboundPeakAmp_KOlist(:,D_ii)-1)*100,['*-' Color_Letter(D_ii)]); hold on
    Leg{D_ii} = ['dur = ' num2str(Input_I_dur_list(D_ii))];
end
legend(Leg)
title('percent different of Rebound Peak value ')
xlabel('Input Current amplitude'); ylabel('percent different of Rebound Peak value ');



plot(Input_I_amp_list, (ReboundPeakAmp_WTlist(:,1)./ReboundPeakAmp_KOlist(:,1)-1)*100,'*-'); hold on
plot(Input_I_amp_list, (ReboundPeakAmp_WTlist(:,2)./ReboundPeakAmp_KOlist(:,2)-1)*100,'*-r')
legend(['dur = ' num2str(Input_I_dur_list(1))] ,['dur = ' num2str(Input_I_dur_list(2))])
title('percent different of Rebound Peak value ')
xlabel('Input Current amplitude'); ylabel('percent different of Rebound Peak value ');

%peak different
figure;
plot(Input_I_amp_list, ReboundPeakAmp_WTlist(:,1) - ReboundPeakAmp_KOlist(:,1),'*-'); hold on
plot(Input_I_amp_list, ReboundPeakAmp_WTlist(:,2) - ReboundPeakAmp_KOlist(:,2),'*-r')
legend(['dur = ' num2str(Input_I_dur_list(1))] ,['dur = ' num2str(Input_I_dur_list(2))])
title('Different of rebound peak value')
xlabel('Input Current amplitude'); ylabel('Different of Rebound Peak value ');

% delay percent different  --> KO = 1
figure;
plot(Input_I_amp_list, (ReboundPeakDel_WTlist(:,1)./ReboundPeakDel_KOlist(:,1)-1)*100,'*-'); hold on
plot(Input_I_amp_list, (ReboundPeakDel_WTlist(:,2)./ReboundPeakDel_KOlist(:,2)-1)*100,'*-r')
legend(['dur = ' num2str(Input_I_dur_list(1))] ,['dur = ' num2str(Input_I_dur_list(2))])
title('percent different of Rebound Peak delay ')
xlabel('Input Current amplitude'); ylabel('percent different of Rebound Peak delay');

%delay different
figure;
plot(Input_I_amp_list, ReboundPeakDel_WTlist(:,1) - ReboundPeakDel_KOlist(:,1),'*-'); hold on
plot(Input_I_amp_list, ReboundPeakDel_WTlist(:,2) - ReboundPeakDel_KOlist(:,2),'*-r')
legend(['dur = ' num2str(Input_I_dur_list(1))] ,['dur = ' num2str(Input_I_dur_list(2))])
title('Different of rebound peak delay');
xlabel('Input Current amplitude'); ylabel('Different of Rebound Peak delay');

%Scatter of rebound peak and delay
figure;
scatter(ReboundPeakAmp_WTlist(:,1), ReboundPeakDel_WTlist(:,1),'*k'); hold on;
scatter(ReboundPeakAmp_WTlist(:,2), ReboundPeakDel_WTlist(:,2),'ob'); hold on;
scatter(ReboundPeakAmp_KOlist(:,1), ReboundPeakDel_KOlist(:,1),'*m'); hold on;
scatter(ReboundPeakAmp_KOlist(:,2), ReboundPeakDel_KOlist(:,2),'or'); hold on;
legend(['WT, dur = ' num2str(Input_I_dur_list(1))], ['WT, dur = ' num2str(Input_I_dur_list(2))], ['KO, dur = ' num2str(Input_I_dur_list(1))],['KO, dur = ' num2str(Input_I_dur_list(2))])
xlabel('Rebound Peak'); ylabel('Reboun Delay')
title('Scatter of Rebound spike peak and delay')

%Scatter of rebound peak and delay percent different
figure;
scatter((ReboundPeakAmp_WTlist(:,1)./ReboundPeakAmp_KOlist(:,1)-1)*100, (ReboundPeakDel_WTlist(:,1)./ReboundPeakDel_KOlist(:,1)-1)*100); hold on;
scatter((ReboundPeakAmp_WTlist(:,2)./ReboundPeakAmp_KOlist(:,2)-1)*100, (ReboundPeakDel_WTlist(:,2)./ReboundPeakDel_KOlist(:,2)-1)*100); hold on;
legend(['dur = ' num2str(Input_I_dur_list(1))] ,['dur = ' num2str(Input_I_dur_list(2))])
xlabel('Rebound Peak'); ylabel('Reboun Delay')
title('Scatter of Rebound spike peak and delay percent difference')

end

if(0)
    for ff = 4 : 5: 128
        figure(ff)
        k = waitforbuttonpress;
        figure(ff+1)
        k = waitforbuttonpress;
    end
end

%% Manually Fixed Peaks
FIXED = 0;
if(FIXED)
    LIST = 0;
    SCATTERs = 0;
    if(LIST)
        for A_ii = 1: length(Input_I_amp_list)
            for D_ii = 1: length(Input_I_dur_list)
                disp(['( ' num2str(A_ii) ',' num2str(D_ii) ' ) I amp = ' num2str(Input_I_amp_list(A_ii)) ', dur = ' num2str(Input_I_dur_list(D_ii))  ])
            end
        end
        
    end
    
    A_ii = 4; D_ii = 3;
    Input_I_amp = Input_I_amp_list(A_ii); Input_I_dur = Input_I_dur_list(D_ii);
    disp(['( ' num2str(A_ii) ',' num2str(D_ii) ' ) I amp = ' num2str(Input_I_amp_list(A_ii)) ', dur = ' num2str(Input_I_dur_list(D_ii))  ])
    InjStopT = InjStartT + Input_I_dur;
    avgFR_WT_All = Save_BasalAct{A_ii, D_ii}.WT.avgFR;
    avgFR_KO_All = Save_BasalAct{A_ii, D_ii}.KO.avgFR;
    
    fg_Allcut = figure;  set(fg_Allcut, 'position', [  221         443        1606         475]);
    plot(avgFR_WT_All,'k','LineWidth',3);  hold on;
    plot(avgFR_KO_All,'m','LineWidth',3);  hold on;
    legend('WT','KO')
    % title(['All cell, Inject I : amp = ' num2str(Input_I_amp) ', dur = ' num2str(Input_I_dur) ])
    ylabel('Average firing rate (Hz)')
    xlabel('Time(ms)')
    zoomx1 = 100; zoomx2 = 2000;
    xlim([zoomx1 zoomx2])
    ylim( [0 max([ avgFR_WT_All avgFR_WT_All])+2])
    legend('WT','KO')
    
    % PeakAvgFr_WT_All = max(avgFR_WT_All(zoomx1:zoomx2));
    [pks, locs] = findpeaks(avgFR_WT_All);
    if(SCATTERs)
        scatter(locs,pks)
    end
    pks =  pks(locs>InjStopT); locs = locs(locs > InjStopT); locs = locs(2); pks =  pks(2); % Find the first peak
    peakvalWT_All = pks; peakTimeWT_All = locs;
    [pks, locs] = findpeaks(avgFR_KO_All);
    if(SCATTERs)
        scatter(locs,pks)
    end
    pks =  pks(locs>InjStopT); locs = locs(locs > InjStopT); locs = locs(4); pks =  pks(4); % Find the first peak
    peakvalKO_All = pks; peakTimeKO_All = locs;
    % [peakvalWT_All, peakTimeWT_All] = max(avgFR_WT_All(InjStopT:InjStopT+burstRange)); % Consider finding after stop injection negative current only
    % [peakvalKO_All, peakTimeKO_All] = max(avgFR_KO_All(InjStopT:InjStopT+burstRange));
    % peakTimeWT_All = peakTimeWT_All + InjStopT ;
    % peakTimeKO_All = peakTimeKO_All + InjStopT ;
    % figure(17); hold on;
    scatter(peakTimeWT_All, peakvalWT_All,150,'*k','LineWidth',2); peakTimeWT_All = peakTimeWT_All-1; hold on; %plot start at time 0
    scatter(peakTimeKO_All, peakvalKO_All,150,'*m','LineWidth',2); peakTimeKO_All = peakTimeKO_All -1;
    rectangle('Position',[InjStartT , 0,Input_I_dur,  max([ avgFR_WT_All avgFR_WT_All])+2],'EdgeColor','r');
    set(gca,'box','off')
    title({['Average firing rate Fr(t) cell in WT and KO , Inject I : amp = ' num2str(Input_I_amp) ', dur = ' num2str(Input_I_dur)],['WT rebound peak = ' num2str(peakvalWT_All) ' Hz at ' num2str(peakTimeWT_All) ' ms , KO rebound peak = ' num2str(peakvalKO_All) ' Hz at ' num2str(peakTimeKO_All) ' ms']});
    
    %%
    ReboundPeakAmp_WTlist(A_ii, D_ii) = peakvalWT_All ; ReboundPeakDel_WTlist(A_ii, D_ii) = peakTimeWT_All;
    ReboundPeakAmp_KOlist(A_ii, D_ii) = peakvalKO_All; ReboundPeakDel_KOlist(A_ii, D_ii) = peakTimeKO_All;
    disp(['WRITE DOWN : ( ' num2str(A_ii) ',' num2str(D_ii) ' ) I amp = ' num2str(Input_I_amp_list(A_ii)) ', dur = ' num2str(Input_I_dur_list(D_ii))  ])
    
    
end


if(0)
    save([dirLoc 'Data_amp-1_-5_dur10_50_Trial2.mat'], 'Save_BasalAct', 'ReboundPeakAmp_WTlist','ReboundPeakAmp_KOlist','ReboundPeakDel_WTlist','ReboundPeakDel_KOlist');
end


%% Parameter Matrix

if(0)
tmpTest = zeros(length(Input_I_amp_list),length(Input_I_dur_list));
cnt = 0;
for A_ii = 1: length(Input_I_amp_list)
    for D_ii = 1: length(Input_I_dur_list)
        tmpTest(A_ii,D_ii) = cnt; cnt = cnt +1;
        disp(['Cnt# ' num2str(cnt) ' ( ' num2str(A_ii) ',' num2str(D_ii) ' ) I amp = ' num2str(Input_I_amp_list(A_ii)) ', dur = ' num2str(Input_I_dur_list(D_ii))  ])
    end
end

figure; imagesc(tmpTest); colorbar(); % caxis([0 0.09])
axis equal xy
xlabel('Duration of Current Injection (ms)')
ylabel('Amplitude of Current Injection (nA)')

for ii = 1:length(Input_I_amp_list)
    for jj = 1:length(Input_I_dur_list)
        text(jj-0.25,ii,{['cnt:' num2str(tmpTest(ii,jj))],['Amp:' num2str(Input_I_amp_list(ii))],['Dur:' num2str(Input_I_dur_list(jj))]},'Color','w')
    end
end

PercentPeakDiff = (ReboundPeakAmp_WTlist./ReboundPeakAmp_KOlist-1)*100;
Paramfg = figure; set(Paramfg, 'position',[242          45        1453         952]);
imagesc(PercentPeakDiff); colorbar(); % caxis([0 0.09])
axis equal xy
xlabel('Duration of Current Injection (ms)')
ylabel('Amplitude of Current Injection (nA)')
for ii = 1:length(Input_I_amp_list)
    for jj = 1:length(Input_I_dur_list)
        tmpTxt = {['Peak Diff:' num2str(PercentPeakDiff(ii,jj)) '%'],['Peak:(' num2str(ReboundPeakAmp_WTlist(ii,jj)) ',' num2str(ReboundPeakAmp_KOlist(ii,jj)) ')'],['Del:(' num2str(ReboundPeakDel_WTlist(ii,jj)) ',' num2str(ReboundPeakDel_KOlist(ii,jj)) ')']};
        if(PercentPeakDiff(ii,jj) > 35) && (PercentPeakDiff(ii,jj) < 60)
            text(jj-0.4,ii,tmpTxt,'Color','k')
        else
            text(jj-0.4,ii,tmpTxt,'Color','w')
        end
    end
end
set(gca,'YTick',[1:length(Input_I_amp_list)])
set(gca,'YTickLabel',Input_I_amp_list)
set(gca,'XTick',[1:length(Input_I_dur_list')])
set(gca,'XTickLabel',Input_I_dur_list)
title({['Parameter Matrix : Percent Peak Different (WT/KO*100%)'] ,'Peak = rebounding spikes peak amplitude(Hz), Del = rebounding spikes peak delay(ms), (Val1,Val2) = rebounding peak of (WT,KO)'});
end
%%



%%


















%%
if(0)
    %Use the same length when testing other wise it is hard on comparison
    maxFittingEnd = max(peakTimeWT_All,peakTimeKO_All);
    % WT
    % x_WT=(InjStopT-200:peakTimeWT_All+1); % Explanatory variable
    x_WT=(InjStopT-200:maxFittingEnd+1);
    y_WT = avgFR_WT_All(x_WT);
    [param_WT,stat_WT]=sigm_fit(x_WT,y_WT) % param = [min, max, x50, slope]
    F = param_WT;
    fit_WT = F(1)+(F(2)-F(1))./(1+10.^((F(3)-x_WT)*F(4)));
    %KO
    % x_KO=(InjStopT-200:peakTimeKO_All+1); % Explanatory variable
    x_KO=(InjStopT-200:maxFittingEnd+1);
    y_KO = avgFR_KO_All(x_KO);
    [param_KO,stat_KO]=sigm_fit(x_KO,y_KO) % param = [min, max, x50, slope]
    F = param_KO;
    fit_KO = F(1)+(F(2)-F(1))./(1+10.^((F(3)-x_KO)*F(4)));
    
    % figure;plot(x_WT,y_WT,'.k'); hold on; plot(x_WT,fit_WT,'b','LineWidth',2);
    % plot(x_KO,y_KO,'.m'); hold on; plot(x_KO,fit_KO,'r','LineWidth',2)
    
    % figure; plot(avgFR_WT_All,'k','LineWidth',2); hold on;
    %  plot(avgFR_KO_All,'m','LineWidth',2);
    % hold on; scatter(x_WT,fit_WT,'.b'); scatter(param_WT(3),avgFR_WT_All(round(param_WT(3))),'*b');
    % scatter(x_KO,fit_KO,'.r'); scatter(param_KO(3),avgFR_KO_All(round(param_KO(3))),'*r');
    
    
    fzoomSlope = figure; set(fzoomSlope, 'position', [  500   558   740   420]);
    plot(avgFR_WT_All,'k','LineWidth',2); hold on;
    plot(avgFR_KO_All,'m','LineWidth',2);
    hold on; plot(x_WT,fit_WT,':b','LineWidth',2);  plot(x_KO,fit_KO,':r','LineWidth',2);
    legend('WT','KO', 'fit WT','fit KO'); %legend('WT fr(t)','KO fr(t)', 'fit WT','fit KO')
    scatter(param_WT(3),avgFR_WT_All(round(param_WT(3))),'*b','LineWidth',3); scatter(param_KO(3),avgFR_KO_All(round(param_KO(3))),'*r','LineWidth',3);
    xlim([InjStopT-200 max(peakTimeWT_All+1,peakTimeKO_All+1)+100]);
    text(param_WT(3)- 170, avgFR_WT_All(round(param_WT(3))), [ ' slope = ' num2str(param_WT(4))],'color','b');
    text(param_KO(3)+20, avgFR_KO_All(round(param_KO(3))), [ ' slope = ' num2str(param_KO(4))],'color','r');
    scatter(peakTimeWT_All+1, peakvalWT_All,150,'*k','LineWidth',2); text(peakTimeWT_All-100, peakvalWT_All+1.5,150,['(t=' num2str(peakTimeWT_All) ' ms, Peak =' num2str(peakvalWT_All) ' Hz)'], 'color','k');
    scatter(peakTimeKO_All+1, peakvalKO_All,150,'*m','LineWidth',2); text(peakTimeKO_All-100, peakvalKO_All+1.5,150,['(t=' num2str(peakTimeKO_All) ' ms, Peak =' num2str(peakvalKO_All) ' Hz)'], 'color','m');
    title(['Measuring slope of rebound spikes : WT = '  num2str(param_WT(4)) , ', KO = '  num2str(param_KO(4)) ])
    xlabel('Time (ms)'); ylabel('Average Firing rate (Hz)');
    ylim([-5 25])
    
    
    % % Define function that will be used to fit data
    % % (F is a vector of fitting parameters)
    % f = @(F,x) F(1)+(F(2)-F(1))./(1+10.^((F(3)-x)*F(4))) ; %F(1) + F(2).*x + F(3).*x.^2;
    % F_fitted = nlinfit(x,y,f,[1 1 1 1]);
    % % Display fitted coefficients
    % disp(['F = ',num2str(F_fitted)])
    % % Plot the data and fit
    % figure;
    % plot(x,y,'*',x,f(F_fitted,x),'g');
    % legend('data','fit')
    
    %  y=1/(1+e^(-x)).
    % param(1)+(param(2)-param(1))./(1+10.^((param(3)-xval)*param(4)))
end

%% Autocorrelation During Burstng

       %  Bursting 
    boxRange = 40; T1 = InjStopT + reboundDelay ; T2 = T1 + burstRange;
    boxRangeM1 = 10; 
    histbox = -boxRange:1:boxRange;
    INJECT = 1;
    BURST= 0;
    if(INJECT)
        [box_WT_VL, Npeak_WT_VL,sumISI_WT_VL]=AutocorrSpikeInterval2(WT.VL.inject.spkMatCut,boxRange, 'WT VL ',0);
    [box_WT_M1, Npeak_WT_M1,sumISI_WT_M1]=AutocorrSpikeInterval2(WT.M1.inject.spkMatCut,boxRangeM1, 'WT M1 ',0);
    [box_KO_VL, Npeak_KO_VL,sumISI_KO_VL]=AutocorrSpikeInterval2(KO.VL.inject.spkMatCut,boxRange, 'KO VL ',0);
    [box_KO_M1, Npeak_KO_M1,sumISI_KO_M1]=AutocorrSpikeInterval2(KO.M1.inject.spkMatCut,boxRangeM1, 'KO M1 ',0);
    end
    
    if(BURST)
    [box_WT_VL, Npeak_WT_VL,sumISI_WT_VL]=AutocorrSpikeInterval2(WT.VL.burst.spkMatCut,boxRange, 'WT VL ',0);
    [box_WT_M1, Npeak_WT_M1,sumISI_WT_M1]=AutocorrSpikeInterval2(WT.M1.burst.spkMatCut,boxRangeM1, 'WT M1 ',0);
    [box_KO_VL, Npeak_KO_VL,sumISI_KO_VL]=AutocorrSpikeInterval2(KO.VL.burst.spkMatCut,boxRange, 'KO VL ',0);
    [box_KO_M1, Npeak_KO_M1,sumISI_KO_M1]=AutocorrSpikeInterval2(KO.M1.burst.spkMatCut,boxRangeM1, 'KO M1 ',0);
    end
    %VL
    [meu_WT_VL,sig_WT_VL] = normfit(sumISI_WT_VL);
    Npeak_WT_VL=hist(sumISI_WT_VL,histbox);
    [meu_KO_VL,sig_KO_VL] = normfit(sumISI_KO_VL);
    Npeak_KO_VL=hist(sumISI_KO_VL,histbox);
    dd_WT_VL= -3.5*sig_WT_VL:sig_WT_VL*3.5;
    dd_KO_VL= -3.5*sig_KO_VL:sig_KO_VL*3.5;
    
    ffE = figure; set(ffE, 'position',[ 68         208        1822         686]);
    
    %VL    
    ffE_sub1 = subplot(121);
    hist(sumISI_WT_VL,histbox); set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on;
    % plot(dd_WT_VL,max(Npeak_WT_VL)*func_gauss(dd_WT_VL, sig_WT_VL),'r','LineWidth',2); hold on;
    width_WT_VL = 2.355*sig_WT_VL; %width at half of maximum, normal fit
    hist(sumISI_KO_VL,histbox,'m'); hold on; tmp = get(gca,'child'); set(tmp(1),'FaceColor','m','EdgeColor','k');
    % plot(dd_KO_VL,max(Npeak_KO_VL)*func_gauss(dd_KO_VL, sig_KO_VL),'b','LineWidth',2); hold on;
    width_KO_VL = 2.355*sig_KO_VL; %width at half of maximum, normal fit
    % title({['WT: [mu,sig] = [' num2str(meu_WT_VL) ',' num2str(sig_WT_VL) '], width = ' num2str(width_WT_VL)],[ ' KO: [mu,sig] = [' num2str(meu_KO_VL) ',' num2str(sig_KO_VL) '], width = ' num2str(width_KO_VL) ]}); xlabel('time (ms)'); ylabel('Number of spikes');
    title('VL'); set(gca,'box', 'off')
    
    %M1
    [meu_WT_M1,sig_WT_M1] = normfit(sumISI_WT_M1);
    Npeak_WT_M1=hist(sumISI_WT_M1,histbox);
    [meu_KO_M1,sig_KO_M1] = normfit(sumISI_KO_M1);
    Npeak_KO_M1=hist(sumISI_KO_M1,histbox);
    dd_WT_M1= -3.5*sig_WT_M1:sig_WT_M1*3.5;
    dd_KO_M1= -3.5*sig_KO_M1:sig_KO_M1*3.5;
    
    %M1    
    ffE_sub1 = subplot(122);
    hist(sumISI_WT_M1,histbox); set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on;
    % plot(dd_WT_M1,max(Npeak_WT_M1)*func_gauss(dd_WT_M1, sig_WT_M1),'r','LineWidth',2); hold on;
    width_WT_M1 = 2.355*sig_WT_M1; %width at half of maximum, normal fit
    hist(sumISI_KO_M1,histbox,'m'); hold on; tmp = get(gca,'child'); set(tmp(1),'FaceColor','m','EdgeColor','k');
    % plot(dd_KO_M1,max(Npeak_KO_M1)*func_gauss(dd_KO_M1, sig_KO_M1),'b','LineWidth',2); hold on;
    width_KO_M1 = 2.355*sig_KO_M1; %width at half of maximum, normal fit
    % title({['WT: [mu,sig] = [' num2str(meu_WT_M1) ',' num2str(sig_WT_M1) '], width = ' num2str(width_WT_M1)],[ ' KO: [mu,sig] = [' num2str(meu_KO_M1) ',' num2str(sig_KO_M1) '], width = ' num2str(width_KO_M1) ]}); xlabel('time (ms)'); ylabel('Number of spikes');
    title('M1'); set(gca,'box', 'off')
    
    % normalized with reference spikes at center 
%       normISIpeak = Npeak./Npeak(ceil(length(box)/2)); 
    normISIpeak_WT_VL = Npeak_WT_VL./Npeak_WT_VL(ceil(length(box_WT_VL)/2)); 
    figure; bar(box_WT_VL,normISIpeak_WT_VL); title([' Normalized'])
    xlabel('time (ms)'); ylabel('Spikes Ratio');
    
    
 %%   % Norm + Envelope
    
    ffNN = figure; set(ffNN, 'position',[ 68         208        1822         686]);
    %VL    
    ffNN_sub1 = subplot(121);
    [peak,Time] = hist(sumISI_WT_VL,histbox);
    bar(Time,peak / trapz(Time,peak));
    set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on;
    % plot(dd_WT_VL,max(Npeak_WT_VL)*func_gauss(dd_WT_VL, sig_WT_VL),'r','LineWidth',2); hold on;
    width_WT_VL = 2.355*sig_WT_VL; %width at half of maximum, normal fit
    [peak,Time] = hist(sumISI_KO_VL,histbox,'m');
    bar(Time,peak / trapz(Time,peak));
    hold on; tmp = get(gca,'child'); set(tmp(1),'FaceColor','m','EdgeColor','k');
    % plot(dd_KO_VL,max(Npeak_KO_VL)*func_gauss(dd_KO_VL, sig_KO_VL),'b','LineWidth',2); hold on;
    width_KO_VL = 2.355*sig_KO_VL; %width at half of maximum, normal fit
    % title({['WT: [mu,sig] = [' num2str(meu_WT_E) ',' num2str(sig_WT_E) '], width = ' num2str(width_WT_E)],[ ' KO: [mu,sig] = [' num2str(meu_KO_E) ',' num2str(sig_KO_E) '], width = ' num2str(width_KO_E) ]}); xlabel('time (ms)'); ylabel('Number of spikes');
    [peak,Time] = hist(sumISI_WT_VL,histbox); hold on;
    plot(Time,peak / trapz(Time,peak),'k','LineWidth',3);
    [peak,Time] = hist(sumISI_KO_VL,histbox);
    plot(Time,peak / trapz(Time,peak),'r','LineWidth',3);
    title('VL'); set(gca,'box', 'off')
    
    ffNN_sub1 = subplot(122);
    [peak,Time] = hist(sumISI_WT_M1,histbox);
    bar(Time,peak / trapz(Time,peak));
    set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on;
    % plot(dd_WT_M1,max(Npeak_WT_M1)*func_gauss(dd_WT_M1, sig_WT_M1),'r','LineWidth',2); hold on;
    width_WT_M1 = 2.355*sig_WT_M1; %width at half of maximum, normal fit
    [peak,Time] = hist(sumISI_KO_M1,histbox,'m');
    bar(Time,peak / trapz(Time,peak));
    hold on; tmp = get(gca,'child'); set(tmp(1),'FaceColor','m','EdgeColor','k');
    % plot(dd_KO_M1,max(Npeak_KO_M1)*func_gauss(dd_KO_M1, sig_KO_M1),'b','LineWidth',2); hold on;
    width_KO_M1 = 2.355*sig_KO_M1; %width at half of maximum, normal fit
    % title({['WT: [mu,sig] = [' num2str(meu_WT_E) ',' num2str(sig_WT_E) '], width = ' num2str(width_WT_E)],[ ' KO: [mu,sig] = [' num2str(meu_KO_E) ',' num2str(sig_KO_E) '], width = ' num2str(width_KO_E) ]}); xlabel('time (ms)'); ylabel('Number of spikes');
    [peak,Time] = hist(sumISI_WT_M1,histbox); hold on;
    plot(Time,peak / trapz(Time,peak),'k','LineWidth',3);
    [peak,Time] = hist(sumISI_KO_M1,histbox);
    plot(Time,peak / trapz(Time,peak),'r','LineWidth',3);
    title('M1'); set(gca,'box', 'off')
    
%%     %Normalize with area under curve Envelope
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ffE = figure; set(ffE, 'position',[ 68         208        1822         686]);
    %VL    
    ffE_sub1 = subplot(121);
    [peakWT_VL,TimeWT_VL] = hist(sumISI_WT_VL,histbox);
    plot(TimeWT_VL,peakWT_VL / trapz(TimeWT_VL,peakWT_VL),'k');
    % set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on;
    % plot(dd_WT_VL,max(Npeak_WT_VL)*func_gauss(dd_WT_VL, sig_WT_VL),'r','LineWidth',2); hold on;
    width_WT_VL = 2.355*sig_WT_VL; %width at half of maximum, normal fit
    [peakKO_VL,TimeKO_VL] = hist(sumISI_KO_VL,histbox,'m'); hold on;
    plot(TimeKO_VL,peakKO_VL / trapz(TimeKO_VL,peakKO_VL),'m');
    % hold on; tmp = get(gca,'child'); set(tmp(1),'FaceColor','m','EdgeColor','k');
    % plot(dd_KO_VL,max(Npeak_KO_VL)*func_gauss(dd_KO_VL, sig_KO_VL),'b','LineWidth',2); hold on;
    width_KO_VL = 2.355*sig_KO_VL; %width at half of maximum, normal fit
    % title({['WT: [mu,sig] = [' num2str(meu_WT_VL) ',' num2str(sig_WT_VL) '], width = ' num2str(width_WT_E)],[ ' KO: [mu,sig] = [' num2str(meu_KO_E) ',' num2str(sig_KO_E) '], width = ' num2str(width_KO_E) ]}); xlabel('time (ms)'); ylabel('Number of spikes');
    title('VL'); set(gca,'box', 'off')
    
    %M1
    ffE_sub2 = subplot(122);
    [peakWT_M1,TimeWT_M1] = hist(sumISI_WT_M1,histbox);
    plot(TimeWT_M1,peakWT_M1 / trapz(TimeWT_M1,peakWT_M1),'k');
    % set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on;
    % plot(dd_WT_M1,max(Npeak_WT_M1)*func_gauss(dd_WT_M1, sig_WT_M1),'r','LineWidth',2); hold on;
    width_WT_M1 = 2.355*sig_WT_M1; %width at half of maximum, normal fit
    
    [peak_KO_M1,Time_KO_M1] = hist(sumISI_KO_M1,histbox,'m'); hold on;
    plot(Time_KO_M1,peak_KO_M1 / trapz(Time_KO_M1,peak_KO_M1),'m');
    % tmp = get(gca,'child'); set(tmp(1),'FaceColor','m','EdgeColor','k');
    % plot(dd_KO_M1,max(Npeak_KO_M1)*func_gauss(dd_KO_M1, sig_KO_M1),'b','LineWidth',2); hold on;
    
    width_KO_M1 = 2.355*sig_KO_M1; %width at half of maximum, normal fit
    % title({['WT: [mu,sig] = [' num2str(meu_WT_M1) ',' num2str(sig_WT_M1) '], width = ' num2str(width_WT_M1)],[ ' KO: [mu,sig] = [' num2str(meu_KO_M1) ',' num2str(sig_KO_M1) '], width = ' num2str(width_KO_M1) ]}); xlabel('time (ms)'); ylabel('Number of spikes');
    title('M1'); set(gca,'box', 'off')
    
%     suptitle(['Rebound : T1 = ' num2str(T1) ', T2 = ' num2str(T2) ', histbox Range = ' num2str(boxRange)]);
%     disp(['Rebound : T1 = ' num2str(T1) ', T2 = ' num2str(T2) ', histbox Range = ' num2str(boxRange)]);
    
    % Only All case
    %
    % figure; plot(kernel); %just test delete later
    % yy = func_gauss(edges, sigma); figure; plot(yy*RES) %----> check peak
   
    %% fft of overall fr
    T = size(peakWT_VL,2); % size(WT.E.spktrain,2);
    fftWT_VL = abs(fftshift(fft(peakWT_VL)));
    fftKO_VL = abs(fftshift(fft(peakKO_VL)));
    dHz = 1/T * 1000;
    Hz_ind = [-round(T/2):round(T/2)]*dHz;
    Hz_ind = Hz_ind(1:T);
    
    %smooth  FFT
    gsig = 14;   % filter [Hz] ---------------------------------------> Calculate percent different when changing sigma size
    gfil2 = exp( -[-round(4*gsig):round(4*gsig)].^2 /2/gsig^2);
    gfil2 = gfil2/sum(gfil2(:));
    
    fft_smoothWT_VL = conv(fftWT_VL, gfil2, 'same');
    fft_smoothKO_VL = conv(fftKO_VL, gfil2, 'same');
    fftLimX1 = 1;   fftLimX2 = 10;
    figure;
    hold on
    plot(Hz_ind, fftWT_VL, 'k-')
    plot(Hz_ind, fftKO_VL, 'm-')
    hold on
    %     xlim([fftLimX1 fftLimX2])
    %     xlim([1 5])
    
    %All
    sigma =  SIGMA_fr_All; % ms
    edges =[-3*sigma:1:3*sigma];
    kernel = normpdf(edges,0,sigma);
    kernel = kernel*RES; % Normalized with bin width
    % figure; plot(kernel)
    sWT_All = conv(sum(WT.All.spktrain(1:end,:)),kernel);
    sKO_All = conv(sum(KO.All.spktrain(1:end,:)),kernel);
    
 
    %%
%% fit WT
if(1)
    %Use the same length when testing other wise it is hard on comparison
%     maxFittingEnd = max(peakTimeWT_All,peakTimeKO_All);
    maxFittingEnd = 1200;
    fitX1 = 1040;  fitX2 = maxFittingEnd+1;  %InjStopT+reboundDelay
    % WT
    % x_WT=(InjStopT-200:peakTimeWT_All+1); % Explanatory variable
    x_WT=(fitX1:maxFittingEnd+1);
    y_WT = avgFR_WT_All(x_WT);
    [param_WT,stat_WT]=sigm_fit(x_WT,y_WT) % param = [min, max, x50, slope]
    F = param_WT;
    fit_WT = F(1)+(F(2)-F(1))./(1+10.^((F(3)-x_WT)*F(4)));
    %KO
    % x_KO=(InjStopT-200:peakTimeKO_All+1); % Explanatory variable
     x_WT=(fitX1:maxFittingEnd+1);
    y_KO = avgFR_KO_All(x_KO);
    [param_KO,stat_KO]=sigm_fit(x_KO,y_KO) % param = [min, max, x50, slope]
    F = param_KO;
    fit_KO = F(1)+(F(2)-F(1))./(1+10.^((F(3)-x_KO)*F(4)));
    
    % figure;plot(x_WT,y_WT,'.k'); hold on; plot(x_WT,fit_WT,'b','LineWidth',2);
    % plot(x_KO,y_KO,'.m'); hold on; plot(x_KO,fit_KO,'r','LineWidth',2)
    
    % figure; plot(avgFR_WT_All,'k','LineWidth',2); hold on;
    %  plot(avgFR_KO_All,'m','LineWidth',2);
    % hold on; scatter(x_WT,fit_WT,'.b'); scatter(param_WT(3),avgFR_WT_All(round(param_WT(3))),'*b');
    % scatter(x_KO,fit_KO,'.r'); scatter(param_KO(3),avgFR_KO_All(round(param_KO(3))),'*r');
    
    
    fzoomSlope = figure; set(fzoomSlope, 'position', [  500   558   740   420]);
    plot(avgFR_WT_All,'k','LineWidth',2); hold on;
    plot(avgFR_KO_All,'m','LineWidth',2);
    hold on; plot(x_WT,fit_WT,':b','LineWidth',2);  plot(x_KO,fit_KO,':r','LineWidth',2);
    legend('WT','KO', 'fit WT','fit KO'); %legend('WT fr(t)','KO fr(t)', 'fit WT','fit KO')
    scatter(param_WT(3),avgFR_WT_All(round(param_WT(3))),'*b','LineWidth',3); scatter(param_KO(3),avgFR_KO_All(round(param_KO(3))),'*r','LineWidth',3);
    xlim([fitX1 fitX2]);
    text(param_WT(3)- 170, avgFR_WT_All(round(param_WT(3))), [ ' slope = ' num2str(param_WT(4))],'color','b');
    text(param_KO(3)+20, avgFR_KO_All(round(param_KO(3))), [ ' slope = ' num2str(param_KO(4))],'color','r');
    scatter(peakTimeWT_All+1, peakvalWT_All,150,'*k','LineWidth',2); text(peakTimeWT_All-100, peakvalWT_All+1.5,150,['(t=' num2str(peakTimeWT_All) ' ms, Peak =' num2str(peakvalWT_All) ' Hz)'], 'color','k');
    scatter(peakTimeKO_All+1, peakvalKO_All,150,'*m','LineWidth',2); text(peakTimeKO_All-100, peakvalKO_All+1.5,150,['(t=' num2str(peakTimeKO_All) ' ms, Peak =' num2str(peakvalKO_All) ' Hz)'], 'color','m');
    title(['Measuring slope of rebound spikes : WT = '  num2str(param_WT(4)) , ', KO = '  num2str(param_KO(4)) ])
    xlabel('Time (ms)'); ylabel('Average Firing rate (Hz)');
%     ylim([-5 25])
    
    
    % % Define function that will be used to fit data
    % % (F is a vector of fitting parameters)
    % f = @(F,x) F(1)+(F(2)-F(1))./(1+10.^((F(3)-x)*F(4))) ; %F(1) + F(2).*x + F(3).*x.^2;
    % F_fitted = nlinfit(x,y,f,[1 1 1 1]);
    % % Display fitted coefficients
    % disp(['F = ',num2str(F_fitted)])
    % % Plot the data and fit
    % figure;
    % plot(x,y,'*',x,f(F_fitted,x),'g');
    % legend('data','fit')
    
    %  y=1/(1+e^(-x)).
    % param(1)+(param(2)-param(1))./(1+10.^((param(3)-xval)*param(4)))
end
%%
if(0)
    %% Autocorrelation During Burstng (500ms after stop) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Caution !  - Do not use these results as it is not finalized yet.
%     InjStartT = 1000;
%     InjStopT = InjStartT+Input_I_dur;
%     reboundDelay = 70;
%     burstRange = 100;

    %  Bursting 
    boxRange = 20; T1 = InjStopT + reboundDelay ; T2 = T1 + burstRange;
    histbox = -boxRange:2:boxRange;
    
    
    [Npeak_WT_VL,width_WT_VL,sumISI_WT_VL]=AutocorrSpikeInterval(WT.VL.burst.spkMatCut,boxRange, 'WT VL ',0);
    [Npeak_WT_M1,width_WT_M1,sumISI_WT_M1]=AutocorrSpikeInterval(WT.M1.burst.spkMatCut,boxRange, 'WT M1 ',0);
    [Npeak_KO_VL,width_KO_VL,sumISI_KO_VL]=AutocorrSpikeInterval(KO.VL.burst.spkMatCut,boxRange, 'KO VL ',0);
    [Npeak_KO_M1,width_KO_M1,sumISI_KO_M1]=AutocorrSpikeInterval(KO.M1.burst.spkMatCut,boxRange, 'KO M1 ',0);
 
    %Raw
    
    
    % close all
    boxRange = 30; T1 = 1000; T2 = 1200;
    histbox = -boxRange:2:boxRange;
    
    % [sumISI,histbox] = AutocorrSpikeInterval(spktrainVec, Trange, ttl,cutTime )
    
    [Npeak_WT_E,width_WT_E,sumISI_WT_E]=AutocorrSpikeInterval(WT.E.spktrain(:,T1:T2),boxRange, 'WT E',0);
    [Npeak_WT_I, width_WT_I,sumISI_WT_I]=AutocorrSpikeInterval(WT.I.spktrain(:,T1:T2),boxRange, 'WT I',0);
    [Npeak_KO_E, width_KO_E,sumISI_KO_E]=AutocorrSpikeInterval(KO.E.spktrain(:,T1:T2),boxRange, 'KO E',0);
    [Npeak_KO_I, width_KO_I,sumISI_KO_I]=AutocorrSpikeInterval(KO.I.spktrain(:,T1:T2),boxRange, 'KO I',0);
    [Npeak_WT_All, width_WT_All,sumISI_WT_All]=AutocorrSpikeInterval(WT.All.spktrain(:,T1:T2),boxRange, 'WT All',0);
    [Npeak_KO_All, width_KO_All,sumISI_KO_All]=AutocorrSpikeInterval(KO.All.spktrain(:,T1:T2),boxRange, 'KO All',0);
    
    %E
    [meu_WT_E,sig_WT_E] = normfit(sumISI_WT_E);
    Npeak_WT_E=hist(sumISI_WT_E,histbox);
    [meu_KO_E,sig_KO_E] = normfit(sumISI_KO_E);
    Npeak_KO_E=hist(sumISI_KO_E,histbox);
    dd_WT_E= -3.5*sig_WT_E:sig_WT_E*3.5;
    dd_KO_E= -3.5*sig_KO_E:sig_KO_E*3.5;
    %I
    [meu_WT_I,sig_WT_I] = normfit(sumISI_WT_I);
    Npeak_WT_I=hist(sumISI_WT_I,histbox);
    [meu_KO_I,sig_KO_I] = normfit(sumISI_KO_I);
    Npeak_KO_I=hist(sumISI_KO_I,histbox);
    dd_WT_I= -3.5*sig_WT_I:sig_WT_I*3.5;
    dd_KO_I= -3.5*sig_KO_I:sig_KO_I*3.5;
    %All
    [meu_WT_All,sig_WT_All] = normfit(sumISI_WT_All);
    Npeak_WT_All=hist(sumISI_WT_All,histbox);
    [meu_KO_All,sig_KO_All] = normfit(sumISI_KO_All);
    Npeak_KO_All=hist(sumISI_KO_All,histbox);
    dd_WT_All= -3.5*sig_WT_All:sig_WT_All*3.5;
    dd_KO_All= -3.5*sig_KO_All:sig_KO_All*3.5;
    
    
    ffE = figure; set(ffE, 'position',[ 68         208        1822         686]);
    %E
    
    ffE_sub1 = subplot(131);
    hist(sumISI_WT_E,histbox); set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on;
    % plot(dd_WT_E,max(Npeak_WT_E)*func_gauss(dd_WT_E, sig_WT_E),'r','LineWidth',2); hold on;
    width_WT_E = 2.355*sig_WT_E; %width at half of maximum, normal fit
    hist(sumISI_KO_E,histbox,'m'); hold on; tmp = get(gca,'child'); set(tmp(1),'FaceColor','m','EdgeColor','k');
    % plot(dd_KO_E,max(Npeak_KO_E)*func_gauss(dd_KO_E, sig_KO_E),'b','LineWidth',2); hold on;
    width_KO_E = 2.355*sig_KO_E; %width at half of maximum, normal fit
    % title({['WT: [mu,sig] = [' num2str(meu_WT_E) ',' num2str(sig_WT_E) '], width = ' num2str(width_WT_E)],[ ' KO: [mu,sig] = [' num2str(meu_KO_E) ',' num2str(sig_KO_E) '], width = ' num2str(width_KO_E) ]}); xlabel('time (ms)'); ylabel('Number of spikes');
    title('E'); set(gca,'box', 'off')
    
    %I
    ffE_sub2 = subplot(132);
    hist(sumISI_WT_I,histbox); set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on;
    % plot(dd_WT_I,max(Npeak_WT_I)*func_gauss(dd_WT_I, sig_WT_I),'r','LineWidth',2); hold on;
    width_WT_I = 2.355*sig_WT_I; %width at half of maximum, normal fit
    
    hist(sumISI_KO_I,histbox,'m'); hold on; tmp = get(gca,'child'); set(tmp(1),'FaceColor','m','EdgeColor','k');
    % plot(dd_KO_I,max(Npeak_KO_I)*func_gauss(dd_KO_I, sig_KO_I),'b','LineWidth',2); hold on;
    
    width_KO_I = 2.355*sig_KO_I; %width at half of maximum, normal fit
    % title({['WT: [mu,sig] = [' num2str(meu_WT_I) ',' num2str(sig_WT_I) '], width = ' num2str(width_WT_I)],[ ' KO: [mu,sig] = [' num2str(meu_KO_I) ',' num2str(sig_KO_I) '], width = ' num2str(width_KO_I) ]}); xlabel('time (ms)'); ylabel('Number of spikes');
    title('I'); set(gca,'box', 'off')
    
    %All
    ffE_sub3 = subplot(133);
    hist(sumISI_WT_All,histbox); set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on;
    % plot(dd_WT_All,max(Npeak_WT_All)*func_gauss(dd_WT_All, sig_WT_All),'r','LineWidth',2); hold on;
    width_WT_All = 2.355*sig_WT_All; %width at half of maximum, normal fit
    
    hist(sumISI_KO_All,histbox,'m'); hold on; tmp = get(gca,'child'); set(tmp(1),'FaceColor','m','EdgeColor','k');
    % plot(dd_KO_All,max(Npeak_KO_All)*func_gauss(dd_KO_All, sig_KO_All),'b','LineWidth',2); hold on;
    
    width_KO_All = 2.355*sig_KO_All; %width at half of maximum, normal fit
    % title({['WT: [mu,sig] = [' num2str(meu_WT_All) ',' num2str(sig_WT_All) '], width = ' num2str(width_WT_All)],[ ' KO: [mu,sig] = [' num2str(meu_KO_All) ',' num2str(sig_KO_All) '], width = ' num2str(width_KO_All) ]}); xlabel('time (ms)'); ylabel('Number of spikes');
    title('All'); set(gca,'box', 'off')
    
    suptitle(['Rebound : T1 = ' num2str(T1) ', T2 = ' num2str(T2) ', histbox Range = ' num2str(boxRange)]);
    disp(['Rebound : T1 = ' num2str(T1) ', T2 = ' num2str(T2) ', histbox Range = ' num2str(boxRange)]);
    %%  Normalized
    
    ffE = figure; set(ffE, 'position',[ 68         208        1822         686]);
    %E
    
    ffE_sub1 = subplot(131);
    [peak,Time] = hist(sumISI_WT_E,histbox);
    bar(Time,peak / trapz(Time,peak));
    set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on;
    % plot(dd_WT_E,max(Npeak_WT_E)*func_gauss(dd_WT_E, sig_WT_E),'r','LineWidth',2); hold on;
    width_WT_E = 2.355*sig_WT_E; %width at half of maximum, normal fit
    [peak,Time] = hist(sumISI_KO_E,histbox,'m');
    bar(Time,peak / trapz(Time,peak));
    hold on; tmp = get(gca,'child'); set(tmp(1),'FaceColor','m','EdgeColor','k');
    % plot(dd_KO_E,max(Npeak_KO_E)*func_gauss(dd_KO_E, sig_KO_E),'b','LineWidth',2); hold on;
    width_KO_E = 2.355*sig_KO_E; %width at half of maximum, normal fit
    % title({['WT: [mu,sig] = [' num2str(meu_WT_E) ',' num2str(sig_WT_E) '], width = ' num2str(width_WT_E)],[ ' KO: [mu,sig] = [' num2str(meu_KO_E) ',' num2str(sig_KO_E) '], width = ' num2str(width_KO_E) ]}); xlabel('time (ms)'); ylabel('Number of spikes');
    title('E'); set(gca,'box', 'off')
    
    %I
    ffE_sub2 = subplot(132);
    [peak,Time] = hist(sumISI_WT_I,histbox);
    bar(Time,peak / trapz(Time,peak));
    set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on;
    % plot(dd_WT_I,max(Npeak_WT_I)*func_gauss(dd_WT_I, sig_WT_I),'r','LineWidth',2); hold on;
    width_WT_I = 2.355*sig_WT_I; %width at half of maximum, normal fit
    
    [peak,Time] = hist(sumISI_KO_I,histbox,'m');
    bar(Time,peak / trapz(Time,peak));
    hold on; tmp = get(gca,'child'); set(tmp(1),'FaceColor','m','EdgeColor','k');
    % plot(dd_KO_I,max(Npeak_KO_I)*func_gauss(dd_KO_I, sig_KO_I),'b','LineWidth',2); hold on;
    
    width_KO_I = 2.355*sig_KO_I; %width at half of maximum, normal fit
    % title({['WT: [mu,sig] = [' num2str(meu_WT_I) ',' num2str(sig_WT_I) '], width = ' num2str(width_WT_I)],[ ' KO: [mu,sig] = [' num2str(meu_KO_I) ',' num2str(sig_KO_I) '], width = ' num2str(width_KO_I) ]}); xlabel('time (ms)'); ylabel('Number of spikes');
    title('I'); set(gca,'box', 'off')
    
    %All
    ffE_sub3 = subplot(133);
    [peak,Time] = hist(sumISI_WT_All,histbox);
    bar(Time,peak / trapz(Time,peak));
    set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on;
    % plot(dd_WT_All,max(Npeak_WT_All)*func_gauss(dd_WT_All, sig_WT_All),'r','LineWidth',2); hold on;
    width_WT_All = 2.355*sig_WT_All; %width at half of maximum, normal fit
    
    [peak,Time] = hist(sumISI_KO_All,histbox,'m');
    bar(Time,peak / trapz(Time,peak));
    hold on; tmp = get(gca,'child'); set(tmp(1),'FaceColor','m','EdgeColor','k');
    % plot(dd_KO_All,max(Npeak_KO_All)*func_gauss(dd_KO_All, sig_KO_All),'b','LineWidth',2); hold on;
    
    width_KO_All = 2.355*sig_KO_All; %width at half of maximum, normal fit
    % title({['WT: [mu,sig] = [' num2str(meu_WT_All) ',' num2str(sig_WT_All) '], width = ' num2str(width_WT_All)],[ ' KO: [mu,sig] = [' num2str(meu_KO_All) ',' num2str(sig_KO_All) '], width = ' num2str(width_KO_All) ]}); xlabel('time (ms)'); ylabel('Number of spikes');
    title('All'); set(gca,'box', 'off')
    
    suptitle(['Rebound : T1 = ' num2str(T1) ', T2 = ' num2str(T2) ', histbox Range = ' num2str(boxRange)]);
    disp(['Rebound : T1 = ' num2str(T1) ', T2 = ' num2str(T2) ', histbox Range = ' num2str(boxRange)]);
    %% Norm - KT
    
    ffE = figure; set(ffE, 'position',[ 68         208        1822         686]);
    %E
    
    ffE_sub1 = subplot(131);
    [peak,Time] = hist(sumISI_WT_E,histbox);
    bar(Time,peak / trapz(Time,peak));
    set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on;
    % plot(dd_WT_E,max(Npeak_WT_E)*func_gauss(dd_WT_E, sig_WT_E),'r','LineWidth',2); hold on;
    width_WT_E = 2.355*sig_WT_E; %width at half of maximum, normal fit
    [peak,Time] = hist(sumISI_KO_E,histbox,'m');
    bar(Time,peak / trapz(Time,peak));
    hold on; tmp = get(gca,'child'); set(tmp(1),'FaceColor','m','EdgeColor','k');
    % plot(dd_KO_E,max(Npeak_KO_E)*func_gauss(dd_KO_E, sig_KO_E),'b','LineWidth',2); hold on;
    width_KO_E = 2.355*sig_KO_E; %width at half of maximum, normal fit
    % title({['WT: [mu,sig] = [' num2str(meu_WT_E) ',' num2str(sig_WT_E) '], width = ' num2str(width_WT_E)],[ ' KO: [mu,sig] = [' num2str(meu_KO_E) ',' num2str(sig_KO_E) '], width = ' num2str(width_KO_E) ]}); xlabel('time (ms)'); ylabel('Number of spikes');
    title('E'); set(gca,'box', 'off')
    
    %I
    ffE_sub2 = subplot(132);
    [peak,Time] = hist(sumISI_WT_I,histbox);
    bar(Time,peak / trapz(Time,peak));
    set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on;
    % plot(dd_WT_I,max(Npeak_WT_I)*func_gauss(dd_WT_I, sig_WT_I),'r','LineWidth',2); hold on;
    width_WT_I = 2.355*sig_WT_I; %width at half of maximum, normal fit
    
    [peak,Time] = hist(sumISI_KO_I,histbox,'m');
    bar(Time,peak / trapz(Time,peak));
    hold on; tmp = get(gca,'child'); set(tmp(1),'FaceColor','m','EdgeColor','k');
    % plot(dd_KO_I,max(Npeak_KO_I)*func_gauss(dd_KO_I, sig_KO_I),'b','LineWidth',2); hold on;
    
    width_KO_I = 2.355*sig_KO_I; %width at half of maximum, normal fit
    % title({['WT: [mu,sig] = [' num2str(meu_WT_I) ',' num2str(sig_WT_I) '], width = ' num2str(width_WT_I)],[ ' KO: [mu,sig] = [' num2str(meu_KO_I) ',' num2str(sig_KO_I) '], width = ' num2str(width_KO_I) ]}); xlabel('time (ms)'); ylabel('Number of spikes');
    title('I'); set(gca,'box', 'off')
    
    %All
    ffE_sub3 = subplot(133);
    [peak,Time] = hist(sumISI_WT_All,histbox);
    bar(Time,peak / trapz(Time,peak));
    set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on;
    % plot(dd_WT_All,max(Npeak_WT_All)*func_gauss(dd_WT_All, sig_WT_All),'r','LineWidth',2); hold on;
    width_WT_All = 2.355*sig_WT_All; %width at half of maximum, normal fit
    
    [peak,Time] = hist(sumISI_KO_All,histbox,'m');
    bar(Time,peak / trapz(Time,peak));
    hold on; tmp = get(gca,'child'); set(tmp(1),'FaceColor','m','EdgeColor','k');
    % plot(dd_KO_All,max(Npeak_KO_All)*func_gauss(dd_KO_All, sig_KO_All),'b','LineWidth',2); hold on;
    
    width_KO_All = 2.355*sig_KO_All; %width at half of maximum, normal fit
    % title({['WT: [mu,sig] = [' num2str(meu_WT_All) ',' num2str(sig_WT_All) '], width = ' num2str(width_WT_All)],[ ' KO: [mu,sig] = [' num2str(meu_KO_All) ',' num2str(sig_KO_All) '], width = ' num2str(width_KO_All) ]}); xlabel('time (ms)'); ylabel('Number of spikes');
    title('All'); set(gca,'box', 'off')
    
    suptitle(['Rebound : T1 = ' num2str(T1) ', T2 = ' num2str(T2) ', histbox Range = ' num2str(boxRange)]);
    disp(['Rebound : T1 = ' num2str(T1) ', T2 = ' num2str(T2) ', histbox Range = ' num2str(boxRange)]);
    %% Normalize Envelope
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ffE = figure; set(ffE, 'position',[ 68         208        1822         686]);
    %E
    
    
    ffE_sub1 = subplot(131);
    [peak,Time] = hist(sumISI_WT_E,histbox);
    plot(Time,peak / trapz(Time,peak),'k');
    % set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on;
    % plot(dd_WT_E,max(Npeak_WT_E)*func_gauss(dd_WT_E, sig_WT_E),'r','LineWidth',2); hold on;
    width_WT_E = 2.355*sig_WT_E; %width at half of maximum, normal fit
    [peak,Time] = hist(sumISI_KO_E,histbox,'m'); hold on;
    plot(Time,peak / trapz(Time,peak),'m');
    % hold on; tmp = get(gca,'child'); set(tmp(1),'FaceColor','m','EdgeColor','k');
    % plot(dd_KO_E,max(Npeak_KO_E)*func_gauss(dd_KO_E, sig_KO_E),'b','LineWidth',2); hold on;
    width_KO_E = 2.355*sig_KO_E; %width at half of maximum, normal fit
    % title({['WT: [mu,sig] = [' num2str(meu_WT_E) ',' num2str(sig_WT_E) '], width = ' num2str(width_WT_E)],[ ' KO: [mu,sig] = [' num2str(meu_KO_E) ',' num2str(sig_KO_E) '], width = ' num2str(width_KO_E) ]}); xlabel('time (ms)'); ylabel('Number of spikes');
    title('E'); set(gca,'box', 'off')
    
    %I
    ffE_sub2 = subplot(132);
    [peak,Time] = hist(sumISI_WT_I,histbox);
    plot(Time,peak / trapz(Time,peak),'k');
    % set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on;
    % plot(dd_WT_I,max(Npeak_WT_I)*func_gauss(dd_WT_I, sig_WT_I),'r','LineWidth',2); hold on;
    width_WT_I = 2.355*sig_WT_I; %width at half of maximum, normal fit
    
    [peak,Time] = hist(sumISI_KO_I,histbox,'m'); hold on;
    plot(Time,peak / trapz(Time,peak),'m');
    % tmp = get(gca,'child'); set(tmp(1),'FaceColor','m','EdgeColor','k');
    % plot(dd_KO_I,max(Npeak_KO_I)*func_gauss(dd_KO_I, sig_KO_I),'b','LineWidth',2); hold on;
    
    width_KO_I = 2.355*sig_KO_I; %width at half of maximum, normal fit
    % title({['WT: [mu,sig] = [' num2str(meu_WT_I) ',' num2str(sig_WT_I) '], width = ' num2str(width_WT_I)],[ ' KO: [mu,sig] = [' num2str(meu_KO_I) ',' num2str(sig_KO_I) '], width = ' num2str(width_KO_I) ]}); xlabel('time (ms)'); ylabel('Number of spikes');
    title('I'); set(gca,'box', 'off')
    
    %All
    ffE_sub3 = subplot(133);
    [peakWT,TimeWT] = hist(sumISI_WT_All,histbox);
    plot(TimeWT,peakWT / trapz(TimeWT,peakWT),'k');
    % set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on;
    % plot(dd_WT_All,max(Npeak_WT_All)*func_gauss(dd_WT_All, sig_WT_All),'r','LineWidth',2); hold on;
    width_WT_All = 2.355*sig_WT_All; %width at half of maximum, normal fit
    
    [peakKO,TimeKO] = hist(sumISI_KO_All,histbox,'m');
    plot(Time,peakKO / trapz(Time,peakKO),'m');
    % hold on; tmp = get(gca,'child'); set(tmp(1),'FaceColor','m','EdgeColor','k');
    % plot(dd_KO_All,max(Npeak_KO_All)*func_gauss(dd_KO_All, sig_KO_All),'b','LineWidth',2); hold on;
    
    width_KO_All = 2.355*sig_KO_All; %width at half of maximum, normal fit
    % title({['WT: [mu,sig] = [' num2str(meu_WT_All) ',' num2str(sig_WT_All) '], width = ' num2str(width_WT_All)],[ ' KO: [mu,sig] = [' num2str(meu_KO_All) ',' num2str(sig_KO_All) '], width = ' num2str(width_KO_All) ]}); xlabel('time (ms)'); ylabel('Number of spikes');
    title('All'); set(gca,'box', 'off')
    
    suptitle(['Rebound : T1 = ' num2str(T1) ', T2 = ' num2str(T2) ', histbox Range = ' num2str(boxRange)]);
    disp(['Rebound : T1 = ' num2str(T1) ', T2 = ' num2str(T2) ', histbox Range = ' num2str(boxRange)]);
    
    % Only All case
    %
    % figure; plot(kernel); %just test delete later
    % yy = func_gauss(edges, sigma); figure; plot(yy*RES) %----> check peak
    %% fft of overall fr
    T = size(peakWT,2); % size(WT.E.spktrain,2);
    fftWT_E = abs(fftshift(fft(peakWT)));
    fftKO_E = abs(fftshift(fft(peakKO)));
    dHz = 1/T * 1000;
    Hz_ind = [-round(T/2):round(T/2)]*dHz;
    Hz_ind = Hz_ind(1:T);
    
    %smooth  FFT
    gsig = 14;   % filter [Hz] ---------------------------------------> Calculate percent different when changing sigma size
    gfil2 = exp( -[-round(4*gsig):round(4*gsig)].^2 /2/gsig^2);
    gfil2 = gfil2/sum(gfil2(:));
    
    fft_smoothWT_E = conv(fftWT_E, gfil2, 'same');
    fft_smoothKO_E = conv(fftKO_E, gfil2, 'same');
    fftLimX1 = 1;   fftLimX2 = 10;
    figure;
    hold on
    plot(Hz_ind, fftWT_E, 'k-')
    plot(Hz_ind, fftKO_E, 'm-')
    hold on
    %     xlim([fftLimX1 fftLimX2])
    %     xlim([1 5])
    
    %All
    sigma =  SIGMA_fr_All; % ms
    edges =[-3*sigma:1:3*sigma];
    kernel = normpdf(edges,0,sigma);
    kernel = kernel*RES; % Normalized with bin width
    % figure; plot(kernel)
    sWT_All = conv(sum(WT.All.spktrain(1:end,:)),kernel);
    sKO_All = conv(sum(KO.All.spktrain(1:end,:)),kernel);
    
    %% Norm + Envelope
    
    ffE = figure; set(ffE, 'position',[ 68         208        1822         686]);
    %E
    
    ffE_sub1 = subplot(131);
    [peak,Time] = hist(sumISI_WT_E,histbox);
    bar(Time,peak / trapz(Time,peak));
    set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on;
    % plot(dd_WT_E,max(Npeak_WT_E)*func_gauss(dd_WT_E, sig_WT_E),'r','LineWidth',2); hold on;
    width_WT_E = 2.355*sig_WT_E; %width at half of maximum, normal fit
    [peak,Time] = hist(sumISI_KO_E,histbox,'m');
    bar(Time,peak / trapz(Time,peak));
    hold on; tmp = get(gca,'child'); set(tmp(1),'FaceColor','m','EdgeColor','k');
    % plot(dd_KO_E,max(Npeak_KO_E)*func_gauss(dd_KO_E, sig_KO_E),'b','LineWidth',2); hold on;
    width_KO_E = 2.355*sig_KO_E; %width at half of maximum, normal fit
    % title({['WT: [mu,sig] = [' num2str(meu_WT_E) ',' num2str(sig_WT_E) '], width = ' num2str(width_WT_E)],[ ' KO: [mu,sig] = [' num2str(meu_KO_E) ',' num2str(sig_KO_E) '], width = ' num2str(width_KO_E) ]}); xlabel('time (ms)'); ylabel('Number of spikes');
    [peak,Time] = hist(sumISI_WT_E,histbox); hold on;
    plot(Time,peak / trapz(Time,peak),'k','LineWidth',3);
    [peak,Time] = hist(sumISI_KO_E,histbox);
    plot(Time,peak / trapz(Time,peak),'r','LineWidth',3);
    title('E'); set(gca,'box', 'off')
    
    %I
    ffE_sub2 = subplot(132);
    [peak,Time] = hist(sumISI_WT_I,histbox);
    bar(Time,peak / trapz(Time,peak));
    set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on;
    % plot(dd_WT_I,max(Npeak_WT_I)*func_gauss(dd_WT_I, sig_WT_I),'r','LineWidth',2); hold on;
    width_WT_I = 2.355*sig_WT_I; %width at half of maximum, normal fit
    
    [peak,Time] = hist(sumISI_KO_I,histbox,'m');
    bar(Time,peak / trapz(Time,peak));
    hold on; tmp = get(gca,'child'); set(tmp(1),'FaceColor','m','EdgeColor','k');
    % plot(dd_KO_I,max(Npeak_KO_I)*func_gauss(dd_KO_I, sig_KO_I),'b','LineWidth',2); hold on;
    
    width_KO_I = 2.355*sig_KO_I; %width at half of maximum, normal fit
    % title({['WT: [mu,sig] = [' num2str(meu_WT_I) ',' num2str(sig_WT_I) '], width = ' num2str(width_WT_I)],[ ' KO: [mu,sig] = [' num2str(meu_KO_I) ',' num2str(sig_KO_I) '], width = ' num2str(width_KO_I) ]}); xlabel('time (ms)'); ylabel('Number of spikes');
    [peak,Time] = hist(sumISI_WT_I,histbox); hold on;
    plot(Time,peak / trapz(Time,peak),'b','LineWidth',3);
    [peak,Time] = hist(sumISI_KO_I,histbox);
    plot(Time,peak / trapz(Time,peak),'r','LineWidth',3);
    title('I'); set(gca,'box', 'off')
    
    %All
    ffE_sub3 = subplot(133);
    [peak,Time] = hist(sumISI_WT_All,histbox);
    bar(Time,peak / trapz(Time,peak));
    set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on;
    % plot(dd_WT_All,max(Npeak_WT_All)*func_gauss(dd_WT_All, sig_WT_All),'r','LineWidth',2); hold on;
    width_WT_All = 2.355*sig_WT_All; %width at half of maximum, normal fit
    
    [peak,Time] = hist(sumISI_KO_All,histbox,'m');
    bar(Time,peak / trapz(Time,peak));
    hold on; tmp = get(gca,'child'); set(tmp(1),'FaceColor','m','EdgeColor','k');
    % plot(dd_KO_All,max(Npeak_KO_All)*func_gauss(dd_KO_All, sig_KO_All),'b','LineWidth',2); hold on;
    
    width_KO_All = 2.355*sig_KO_All; %width at half of maximum, normal fit
    % title({['WT: [mu,sig] = [' num2str(meu_WT_All) ',' num2str(sig_WT_All) '], width = ' num2str(width_WT_All)],[ ' KO: [mu,sig] = [' num2str(meu_KO_All) ',' num2str(sig_KO_All) '], width = ' num2str(width_KO_All) ]}); xlabel('time (ms)'); ylabel('Number of spikes');
    title('All'); set(gca,'box', 'off')
    
    suptitle(['Rebound : T1 = ' num2str(T1) ', T2 = ' num2str(T2) ', histbox Range = ' num2str(boxRange)]);
    disp(['Rebound : T1 = ' num2str(T1) ', T2 = ' num2str(T2) ', histbox Range = ' num2str(boxRange)]);
    %%
    
    % Plot the number of spikes with the ISI interval
    % %%
    %Test with All
    fplotISI= figure; set(fplotISI, 'position',[ 68         208        1822         686]);
    %E
    subplot(131)
    plot(histbox, Npeak_WT_E,'k'); hold on;
    plot(histbox, Npeak_KO_E,'m'); hold on;
    line([histbox(1) histbox(end)], [min(Npeak_WT_E) min(Npeak_WT_E)]);
    line([histbox(1) histbox(end)], [min(Npeak_KO_E) min(Npeak_KO_E)]);
    title(['E : min WT = ' num2str(min(Npeak_WT_E)) ' , min KO = ' num2str(min(Npeak_KO_E))])
    %I
    subplot(132)
    plot(histbox, Npeak_WT_I,'k'); hold on;
    plot(histbox, Npeak_KO_I,'m'); hold on;
    line([histbox(1) histbox(end)], [min(Npeak_WT_I) min(Npeak_WT_I)]);
    line([histbox(1) histbox(end)], [min(Npeak_KO_I) min(Npeak_KO_I)]);
    title(['I : min WT = ' num2str(min(Npeak_WT_I)) ' , min KO = ' num2str(min(Npeak_KO_I))])
    %All
    subplot(133)
    plot(histbox, Npeak_WT_All,'k'); hold on;
    plot(histbox, Npeak_KO_All,'m'); hold on;
    line([histbox(1) histbox(end)], [min(Npeak_WT_All) min(Npeak_WT_All)]);
    line([histbox(1) histbox(end)], [min(Npeak_KO_All) min(Npeak_KO_All)]);
    title(['All : min WT = ' num2str(min(Npeak_WT_All)) ' , min KO = ' num2str(min(Npeak_KO_All))])
    %%
    %%
    
    % Triangle fit
    %All
    %WT
    [minVal_WT_All,minTime_WT_All] = min(Npeak_WT_All);
    stepTomid_WT_All =  round(length(histbox)/2) - minTime_WT_All;
    cutNpeak_WT_All = Npeak_WT_All(minTime_WT_All: round(length(histbox)/2)+stepTomid_WT_All); cutBox_WT_All = histbox((minTime_WT_All: round(length(histbox)/2)+stepTomid_WT_All));
    [maxVal_WT_All, maxTime_WT_All] = max(cutNpeak_WT_All);
    Height_WT_All = maxVal_WT_All - minVal_WT_All;
    halfHeight_WT_All = Height_WT_All/2;  % get Width at the half height position
    [mv_WT_All, mi_WT_All]= min(abs(cutNpeak_WT_All - (halfHeight_WT_All+minVal_WT_All)));
    minusMI_WT_All = cutNpeak_WT_All(mi_WT_All) - (halfHeight_WT_All+minVal_WT_All);
    [ mi2_WT_All,wid_WT_All, width_WT_All] = get_mi2_Width_triangle(mi_WT_All, maxTime_WT_All, minusMI_WT_All, cutNpeak_WT_All, cutBox_WT_All)  ;
    
    %KO
    [minVal_KO_All,minTime_KO_All] = min(Npeak_KO_All);
    stepTomid_KO_All =  round(length(histbox)/2) - minTime_KO_All;
    cutNpeak_KO_All = Npeak_KO_All(minTime_KO_All: round(length(histbox)/2)+stepTomid_KO_All); cutBox_KO_All = histbox((minTime_KO_All: round(length(histbox)/2)+stepTomid_KO_All));
    [maxVal_KO_All, maxTime_KO_All] = max(cutNpeak_KO_All);
    Height_KO_All = maxVal_KO_All - minVal_KO_All;
    halfHeight_KO_All = Height_KO_All/2;  % get Width at the half height position
    [mv_KO_All, mi_KO_All]= min(abs(cutNpeak_KO_All - (halfHeight_KO_All+minVal_KO_All)));
    minusMI_KO_All = cutNpeak_KO_All(mi_KO_All) - (halfHeight_KO_All+minVal_KO_All);
    [ mi2_KO_All,wid_KO_All, width_KO_All] = get_mi2_Width_triangle(mi_KO_All, maxTime_KO_All, minusMI_KO_All, cutNpeak_KO_All, cutBox_KO_All)  ;
    
    
    figure;
    plot(cutBox_WT_All, cutNpeak_WT_All,'k'); hold on; scatter(cutBox_WT_All, cutNpeak_WT_All,'ok');
    line([cutBox_WT_All(1) cutBox_WT_All(end)], [halfHeight_WT_All+minVal_WT_All halfHeight_WT_All+minVal_WT_All],'linestyle',':', 'color', 'b');
    plot(cutBox_KO_All, cutNpeak_KO_All,'m'); hold on; scatter(cutBox_KO_All, cutNpeak_KO_All,'om');
    line([cutBox_KO_All(1) cutBox_KO_All(end)], [halfHeight_KO_All+minVal_KO_All halfHeight_KO_All+minVal_KO_All],'linestyle',':', 'color', 'r');
    
    % title(['All: WT''s all cell autocorrelaton width = ' num2str(width_WT_All) ', KO''s all cell autocorrelaton width = '  num2str(width_KO_All) ])
    title('All')
    suptitle(['Level of Synchronization from width of autocorrelation'])
    % Note : This method has problem because the triangle assumption is too
    % width fo r the real data.
    
    % %%
    % Triangle fit  full all case
    if(0)
        %E
        %WT
        [minVal_WT_E,minTime_WT_E] = min(Npeak_WT_E);
        stepTomid_WT_E =  round(length(histbox)/2) - minTime_WT_E;
        cutNpeak_WT_E = Npeak_WT_E(minTime_WT_E: round(length(histbox)/2)+stepTomid_WT_E); cutBox_WT_E = histbox((minTime_WT_E: round(length(histbox)/2)+stepTomid_WT_E));
        [maxVal_WT_E, maxTime_WT_E] = max(cutNpeak_WT_E);
        Height_WT_E = maxVal_WT_E - minVal_WT_E;
        halfHeight_WT_E = Height_WT_E/2;  % get Width at the half height position
        [mv_WT_E, mi_WT_E]= min(abs(cutNpeak_WT_E - (halfHeight_WT_E+minVal_WT_E)));
        minusMI_WT_E = cutNpeak_WT_E(mi_WT_E) - (halfHeight_WT_E+minVal_WT_E);
        [ mi2_WT_E,wid_WT_E, width_WT_E] = get_mi2_Width_triangle(mi_WT_E, maxTime_WT_E, minusMI_WT_E, cutNpeak_WT_E, cutBox_WT_All)  ;
        
        %KO
        [minVal_KO_E,minTime_KO_E] = min(Npeak_KO_E);
        stepTomid_KO_E =  round(length(histbox)/2) - minTime_KO_E;
        cutNpeak_KO_E = Npeak_KO_E(minTime_KO_E: round(length(histbox)/2)+stepTomid_KO_E); cutBox_KO_E = histbox((minTime_KO_E: round(length(histbox)/2)+stepTomid_KO_E));
        [maxVal_KO_E, maxTime_KO_E] = max(cutNpeak_KO_E);
        Height_KO_E = maxVal_KO_E - minVal_KO_E;
        halfHeight_KO_E = Height_KO_E/2;  % get Width at the half height position
        [mv_KO_E, mi_KO_E]= min(abs(cutNpeak_KO_E - (halfHeight_KO_E+minVal_KO_E)));
        minusMI_KO_E = cutNpeak_KO_E(mi_KO_E) - (halfHeight_KO_E+minVal_KO_E);
        [ mi2_KO_E,wid_KO_E, width_KO_E] = get_mi2_Width_triangle(mi_KO_E, maxTime_KO_E, minusMI_KO_E, cutNpeak_KO_E, cutBox_KO_E)  ;
        
        %I
        %WT
        [minVal_WT_I,minTime_WT_I] = min(Npeak_WT_I);
        stepTomid_WT_I =  round(length(histbox)/2) - minTime_WT_I;
        cutNpeak_WT_I = Npeak_WT_I(minTime_WT_I: round(length(histbox)/2)+stepTomid_WT_I); cutBox_WT_I = histbox((minTime_WT_I: round(length(histbox)/2)+stepTomid_WT_I));
        [maxVal_WT_I, maxTime_WT_I] = max(cutNpeak_WT_I);
        Height_WT_I = maxVal_WT_I - minVal_WT_I;
        halfHeight_WT_I = Height_WT_I/2;  % get Width at the half height position
        [mv_WT_I, mi_WT_I]= min(abs(cutNpeak_WT_I - (halfHeight_WT_I+minVal_WT_I)));
        minusMI_WT_I = cutNpeak_WT_I(mi_WT_I) - (halfHeight_WT_I+minVal_WT_I);
        [ mi2_WT_I,wid_WT_I, width_WT_I] = get_mi2_Width_triangle(mi_WT_I, maxTime_WT_I, minusMI_WT_I, cutNpeak_WT_I, cutBox_WT_All)  ;
        
        %KO
        [minVal_KO_I,minTime_KO_I] = min(Npeak_KO_I);
        stepTomid_KO_I =  round(length(histbox)/2) - minTime_KO_I;
        cutNpeak_KO_I = Npeak_KO_I(minTime_KO_I: round(length(histbox)/2)+stepTomid_KO_I); cutBox_KO_I = histbox((minTime_KO_I: round(length(histbox)/2)+stepTomid_KO_I));
        [maxVal_KO_I, maxTime_KO_I] = max(cutNpeak_KO_I);
        Height_KO_I = maxVal_KO_I - minVal_KO_I;
        halfHeight_KO_I = Height_KO_I/2;  % get Width at the half height position
        [mv_KO_I, mi_KO_I]= min(abs(cutNpeak_KO_I - (halfHeight_KO_I+minVal_KO_I)));
        minusMI_KO_I = cutNpeak_KO_I(mi_KO_I) - (halfHeight_KO_I+minVal_KO_I);
        [ mi2_KO_I,wid_KO_I, width_KO_I] = get_mi2_Width_triangle(mi_KO_I, maxTime_KO_I, minusMI_KO_I, cutNpeak_KO_I, cutBox_KO_I)  ;
        
        
        %All
        %WT
        [minVal_WT_All,minTime_WT_All] = min(Npeak_WT_All);
        stepTomid_WT_All =  round(length(histbox)/2) - minTime_WT_All;
        cutNpeak_WT_All = Npeak_WT_All(minTime_WT_All: round(length(histbox)/2)+stepTomid_WT_All); cutBox_WT_All = histbox((minTime_WT_All: round(length(histbox)/2)+stepTomid_WT_All));
        [maxVal_WT_All, maxTime_WT_All] = max(cutNpeak_WT_All);
        Height_WT_All = maxVal_WT_All - minVal_WT_All;
        halfHeight_WT_All = Height_WT_All/2;  % get Width at the half height position
        [mv_WT_All, mi_WT_All]= min(abs(cutNpeak_WT_All - (halfHeight_WT_All+minVal_WT_All)));
        minusMI_WT_All = cutNpeak_WT_All(mi_WT_All) - (halfHeight_WT_All+minVal_WT_All);
        [ mi2_WT_All,wid_WT_All, width_WT_All] = get_mi2_Width_triangle(mi_WT_All, maxTime_WT_All, minusMI_WT_All, cutNpeak_WT_All, cutBox_WT_All)  ;
        
        %KO
        [minVal_KO_All,minTime_KO_All] = min(Npeak_KO_All);
        stepTomid_KO_All =  round(length(histbox)/2) - minTime_KO_All;
        cutNpeak_KO_All = Npeak_KO_All(minTime_KO_All: round(length(histbox)/2)+stepTomid_KO_All); cutBox_KO_All = histbox((minTime_KO_All: round(length(histbox)/2)+stepTomid_KO_All));
        [maxVal_KO_All, maxTime_KO_All] = max(cutNpeak_KO_All);
        Height_KO_All = maxVal_KO_All - minVal_KO_All;
        halfHeight_KO_All = Height_KO_All/2;  % get Width at the half height position
        [mv_KO_All, mi_KO_All]= min(abs(cutNpeak_KO_All - (halfHeight_KO_All+minVal_KO_All)));
        minusMI_KO_All = cutNpeak_KO_All(mi_KO_All) - (halfHeight_KO_All+minVal_KO_All);
        [ mi2_KO_All,wid_KO_All, width_KO_All] = get_mi2_Width_triangle(mi_KO_All, maxTime_KO_All, minusMI_KO_All, cutNpeak_KO_All, cutBox_KO_All)  ;
        
        
        fWidthISI = figure; set(fWidthISI, 'position', [ 68         208        1822         686]);
        subplot(131)
        plot(cutBox_WT_E, cutNpeak_WT_E,'k'); hold on; scatter(cutBox_WT_E, cutNpeak_WT_E,'ok');
        line([cutBox_WT_E(1) cutBox_WT_E(end)], [halfHeight_WT_E+minVal_WT_E halfHeight_WT_E+minVal_WT_E],'linestyle',':', 'color', 'b');
        plot(cutBox_KO_E, cutNpeak_KO_E,'m'); hold on; scatter(cutBox_KO_E, cutNpeak_KO_E,'om');
        line([cutBox_KO_E(1) cutBox_KO_E(end)], [halfHeight_KO_E+minVal_KO_E halfHeight_KO_E+minVal_KO_E],'linestyle',':', 'color', 'r');
        % title('E cell')
        disp(['[E] WT''s all cell : width = ' num2str(width_WT_E) ', KO''s all cell : width = '  num2str(width_KO_E) ])
        
        title(['[E] WT''s all cell : width = ' num2str(width_WT_E) ', KO''s all cell : width = '  num2str(width_KO_E) ])
        subplot(132)
        plot(cutBox_WT_I, cutNpeak_WT_I,'k'); hold on; scatter(cutBox_WT_I, cutNpeak_WT_I,'ok');
        line([cutBox_WT_I(1) cutBox_WT_I(end)], [halfHeight_WT_I+minVal_WT_I halfHeight_WT_I+minVal_WT_I],'linestyle',':', 'color', 'b');
        plot(cutBox_KO_I, cutNpeak_KO_I,'m'); hold on; scatter(cutBox_KO_I, cutNpeak_KO_I,'om');
        line([cutBox_KO_I(1) cutBox_KO_I(end)], [halfHeight_KO_I+minVal_KO_I halfHeight_KO_I+minVal_KO_I],'linestyle',':', 'color', 'r');
        % title('I cell')
        disp(['[I] WT''s all cell : width = ' num2str(width_WT_I) ', KO''s all cell : width = '  num2str(width_KO_I) ])
        title(['[I] WT''s all cell : width = ' num2str(width_WT_I) ', KO''s all cell : width = '  num2str(width_KO_I) ])
        
        subplot(133)
        plot(cutBox_WT_All, cutNpeak_WT_All,'k'); hold on; scatter(cutBox_WT_All, cutNpeak_WT_All,'ok');
        line([cutBox_WT_All(1) cutBox_WT_All(end)], [halfHeight_WT_All+minVal_WT_All halfHeight_WT_All+minVal_WT_All],'linestyle',':', 'color', 'b');
        plot(cutBox_KO_All, cutNpeak_KO_All,'m'); hold on; scatter(cutBox_KO_All, cutNpeak_KO_All,'om');
        line([cutBox_KO_All(1) cutBox_KO_All(end)], [halfHeight_KO_All+minVal_KO_All halfHeight_KO_All+minVal_KO_All],'linestyle',':', 'color', 'r');
        % title('All')
        % title(['All: WT''s all cell autocorrelaton width = ' num2str(width_WT_All) ', KO''s all cell autocorrelaton width = '  num2str(width_KO_All) ])
        disp(['[All] WT''s all cell : width = ' num2str(width_WT_All) ', KO''s all cell : width = '  num2str(width_KO_All) ])
        title(['[All] WT''s all cell : width = ' num2str(width_WT_All) ', KO''s all cell : width = '  num2str(width_KO_All) ])
        
        suptitle(['Level of Synchronization from width of autocorrelation'])
        % Note : This method has problem because the triangle assumption is too
        % width fo r the real data.
    end
    %%
    % Roughly
    figure;
    plot(histbox, Npeak_WT_All,'k'); hold on;
    % plot(histbox, Npeak_KO_All,'m'); hold on;
    line([histbox(1) histbox(end)], [min(cutNpeak_WT_All) min(cutNpeak_WT_All)]);
    % line([histbox(1) histbox(end)], [min(Npeak_KO_All) min(Npeak_KO_All)]);
    hold on; scatter(cutBox_WT_All(mi_WT_All),cutNpeak_WT_All(mi_WT_All))
    widthr = 2* abs(cutBox_WT_All(mi_WT_All)-cutBox_WT_All(maxTime_WT_All)); %rough
   
    

    
    
    %%  ISI histogram of all simulation. (after 200)
    % Analysis the ISI during rebounding
    T1 = InjStopT; T2 = InjStopT+200; step = 2; stepE = 3;
    %E
    ReboundSpikeTrain_WT_E =  WT.E.spktrain(:,T1:stepE:T2);
    [ISI_combine_WT_E,ISI_WT_E ]  = getISIfromSpikeTrain(ReboundSpikeTrain_WT_E);
    %I
    ReboundSpikeTrain_WT_I =  WT.I.spktrain(:,T1:step:T2);
    [ISI_combine_WT_I,ISI_WT_I ]  = getISIfromSpikeTrain(ReboundSpikeTrain_WT_I);
    %All
    ReboundSpikeTrain_WT_All =  [WT.E.spktrain(:,T1:step:T2); WT.I.spktrain(:,T1:step:T2)];
    [ISI_combine_WT_All,ISI_WT_All ]  = getISIfromSpikeTrain(ReboundSpikeTrain_WT_All);
    % KO
    %E
    ReboundSpikeTrain_KO_E =  KO.E.spktrain(:,T1:stepE:T2);
    [ISI_combine_KO_E,ISI_KO_E ]  = getISIfromSpikeTrain(ReboundSpikeTrain_KO_E);
    %I
    ReboundSpikeTrain_KO_I =  KO.I.spktrain(:,T1:step:T2);
    [ISI_combine_KO_I,ISI_KO_I ]  = getISIfromSpikeTrain(ReboundSpikeTrain_KO_I);
    %All
    ReboundSpikeTrain_KO_All =  [KO.E.spktrain(:,T1:step:T2); KO.I.spktrain(:,T1:step:T2)];
    [ISI_combine_KO_All,ISI_KO_All ]  = getISIfromSpikeTrain(ReboundSpikeTrain_KO_All);
    
    
    fg_ISIdist = figure;  set(fg_ISIdist , 'position',[ 68         208        1822         686]);
    subplot(131)
    hist(ISI_combine_WT_E,length(min(ISI_combine_WT_E):stepE:max(ISI_combine_WT_E))); set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on;
    hist(ISI_combine_KO_E,length(min(ISI_combine_KO_E):stepE:max(ISI_combine_KO_E))); tmp = get(gca,'child');
    set(tmp(1),'FaceColor','m','EdgeColor','m');
    title('E'); xlabel('Inter-Spike-Interval (ms)'); ylabel('Number of spikes');
    subplot(132)
    hist(ISI_combine_WT_I,length(min(ISI_combine_WT_I):step:max(ISI_combine_WT_I))); set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on; tmp = get(gca,'child');
    hist(ISI_combine_KO_I,length(min(ISI_combine_KO_I):step:max(ISI_combine_KO_I))); tmp = get(gca,'child');
    set(tmp(1),'FaceColor','m','EdgeColor','m');
    title('I'); xlabel('Inter-Spike-Interval (ms)'); ylabel('Number of spikes');
    subplot(133)
    hist(ISI_combine_WT_All,length(min(ISI_combine_WT_All):step:max(ISI_combine_WT_All))); set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on; tmp = get(gca,'child');
    hist(ISI_combine_KO_All,length(min(ISI_combine_KO_All):step:max(ISI_combine_KO_All))); tmp = get(gca,'child');
    set(tmp(1),'FaceColor','m','EdgeColor','m');
    title('All'); xlabel('Inter-Spike-Interval (ms)'); ylabel('Number of spikes');
    
    suptitle('Distribution of Inter-Spike-Interval in the neural network during recover from inhibition input')
    
    %All only
    figure;
    hist(ISI_combine_WT_All,length(min(ISI_combine_WT_All):step:max(ISI_combine_WT_All))); set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on; tmp = get(gca,'child');
    hist(ISI_combine_KO_All,length(min(ISI_combine_KO_All):step:max(ISI_combine_KO_All))); tmp = get(gca,'child');
    set(tmp(1),'FaceColor','m','EdgeColor','m');
    title('All'); xlabel('Inter-Spike-Interval (ms)'); ylabel('Number of spikes');
    mean_ISI_WT = mean(ISI_combine_WT_All); mode_ISI_WT = mode(ISI_combine_WT_All);
    mean_ISI_KO = mean(ISI_combine_KO_All); mode_ISI_KO = mode(ISI_combine_KO_All);
    scatter(mean_ISI_WT,0,'*b','LineWidth',3)
    scatter(mean_ISI_KO,0,'*r', 'LineWidth',3)
    legend('WT', 'KO','mean WT', 'mean KO')
    %% Distribution of ISI ---> Use fixed vector for bin center
    T1 = InjStopT; T2 = InjStopT+200; step = 2; stepE = 2;
    %E
    ReboundSpikeTrain_WT_E =  WT.E.spktrain(:,T1:stepE:T2);
    [ISI_combine_WT_E,ISI_WT_E ]  = getISIfromSpikeTrain(ReboundSpikeTrain_WT_E);
    %I
    ReboundSpikeTrain_WT_I =  WT.I.spktrain(:,T1:step:T2);
    [ISI_combine_WT_I,ISI_WT_I ]  = getISIfromSpikeTrain(ReboundSpikeTrain_WT_I);
    %All
    ReboundSpikeTrain_WT_All =  [WT.E.spktrain(:,T1:step:T2); WT.I.spktrain(:,T1:step:T2)];
    [ISI_combine_WT_All,ISI_WT_All ]  = getISIfromSpikeTrain(ReboundSpikeTrain_WT_All);
    % KO
    %E
    ReboundSpikeTrain_KO_E =  KO.E.spktrain(:,T1:stepE:T2);
    [ISI_combine_KO_E,ISI_KO_E ]  = getISIfromSpikeTrain(ReboundSpikeTrain_KO_E);
    %I
    ReboundSpikeTrain_KO_I =  KO.I.spktrain(:,T1:step:T2);
    [ISI_combine_KO_I,ISI_KO_I ]  = getISIfromSpikeTrain(ReboundSpikeTrain_KO_I);
    %All
    ReboundSpikeTrain_KO_All =  [KO.E.spktrain(:,T1:step:T2); KO.I.spktrain(:,T1:step:T2)];
    [ISI_combine_KO_All,ISI_KO_All ]  = getISIfromSpikeTrain(ReboundSpikeTrain_KO_All);
    
    
    fg_ISIdist = figure;  set(fg_ISIdist , 'position',[ 68         208        1822         686]);
    subplot(131)
    Mhist_E = min(min(ISI_combine_WT_E),min(ISI_combine_KO_E)):stepE:max(max(ISI_combine_WT_E),max(ISI_combine_KO_E));
    hist(ISI_combine_WT_E,Mhist_E); set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on;
    hist(ISI_combine_KO_E,Mhist_E); tmp = get(gca,'child');
    set(tmp(1),'FaceColor','m','EdgeColor','m');
    title('E'); xlabel('Inter-Spike-Interval (ms)'); ylabel('Number of spikes');
    subplot(132)
    Mhist_I = min(min(ISI_combine_WT_I),min(ISI_combine_KO_I)):stepE:max(max(ISI_combine_WT_I),max(ISI_combine_KO_I));
    hist(ISI_combine_WT_I,Mhist_I); set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on; tmp = get(gca,'child');
    hist(ISI_combine_KO_I,Mhist_I); tmp = get(gca,'child');
    set(tmp(1),'FaceColor','m','EdgeColor','m');
    title('I'); xlabel('Inter-Spike-Interval (ms)'); ylabel('Number of spikes');
    subplot(133)
    Mhist_All = min(min(ISI_combine_WT_All),min(ISI_combine_KO_All)):stepE:max(max(ISI_combine_WT_All),max(ISI_combine_KO_All));
    
    hist(ISI_combine_WT_All,Mhist_All); set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on; tmp = get(gca,'child');
    hist(ISI_combine_KO_All,Mhist_All); tmp = get(gca,'child');
    set(tmp(1),'FaceColor','m','EdgeColor','m');
    title('All'); xlabel('Inter-Spike-Interval (ms)'); ylabel('Number of spikes');
    
    suptitle('Distribution of Inter-Spike-Interval in the neural network during recover from inhibition input')
    
    %% Contour
    
    fg_ISIcontour = figure;  set(fg_ISIcontour , 'position',[ 68         208        1822         686]);
    subplot(131)
    Mhist_E = min(min(ISI_combine_WT_E),min(ISI_combine_KO_E)):stepE:max(max(ISI_combine_WT_E),max(ISI_combine_KO_E));
    hist(ISI_combine_WT_E,Mhist_E);  set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on;
    hist(ISI_combine_KO_E,Mhist_E); tmp = get(gca,'child');
    set(tmp(1),'FaceColor','m','EdgeColor','m');
    title('E'); xlabel('Inter-Spike-Interval (ms)'); ylabel('Number of spikes');
    subplot(132)
    Mhist_I = min(min(ISI_combine_WT_I),min(ISI_combine_KO_I)):stepE:max(max(ISI_combine_WT_I),max(ISI_combine_KO_I));
    hist(ISI_combine_WT_I,Mhist_I); set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on; tmp = get(gca,'child');
    hist(ISI_combine_KO_I,Mhist_I); tmp = get(gca,'child');
    set(tmp(1),'FaceColor','m','EdgeColor','m');
    title('I'); xlabel('Inter-Spike-Interval (ms)'); ylabel('Number of spikes');
    subplot(133)
    Mhist_All = min(min(ISI_combine_WT_All),min(ISI_combine_KO_All)):stepE:max(max(ISI_combine_WT_All),max(ISI_combine_KO_All));
    
    hist(ISI_combine_WT_All,Mhist_All); set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    hold on; tmp = get(gca,'child');
    hist(ISI_combine_KO_All,Mhist_All); tmp = get(gca,'child');
    set(tmp(1),'FaceColor','m','EdgeColor','m');
    title('All'); xlabel('Inter-Spike-Interval (ms)'); ylabel('Number of spikes');
    
    suptitle('Distribution of Inter-Spike-Interval in the neural network during recover from inhibition input')
    
    %% CV versus mean  ISI
    
    %All WT
    ISI_MEAN_WT_All  = zeros(1,length(ISI_WT_All ));
    ISI_STD_WT_All = zeros(1,length(ISI_WT_All ));
    ISI_CV_WT_All = zeros(1,length(ISI_WT_All ));
    for cID = 1 : length(ISI_WT_All)
        ISI_MEAN_WT_All(cID)  = mean(ISI_WT_All{cID});
        ISI_STD_WT_All(cID) = std(ISI_WT_All{cID});
        ISI_CV_WT_All(cID) = std(ISI_WT_All{cID})/ mean(ISI_WT_All{cID});
    end
    %All KO
    ISI_MEAN_KO_All  = zeros(1,length(ISI_KO_All ));
    ISI_STD_KO_All = zeros(1,length(ISI_KO_All ));
    ISI_CV_KO_All = zeros(1,length(ISI_KO_All ));
    for cID = 1 : length(ISI_KO_All)
        ISI_MEAN_KO_All(cID)  = mean(ISI_KO_All{cID});
        ISI_STD_KO_All(cID) = std(ISI_KO_All{cID});
        ISI_CV_KO_All(cID) = std(ISI_KO_All{cID})/ mean(ISI_KO_All{cID});
    end
    
    valISI_MEAN_WT_A = ISI_MEAN_WT_All(~(isnan(ISI_MEAN_WT_All) | isinf(ISI_MEAN_WT_All) )); % 64 Data points
    valISI_CV_WT_A = ISI_CV_WT_All(~(isnan(ISI_CV_WT_All) | isinf(ISI_CV_WT_All)));
    valISI_MEAN_KO_A = ISI_MEAN_KO_All(~(isnan(ISI_MEAN_KO_All) | isinf(ISI_MEAN_KO_All) )); % 56 Data point
    valISI_CV_KO_A = ISI_CV_KO_All(~(isnan(ISI_CV_KO_All) | isinf(ISI_CV_KO_All)));
    %1 regular case
    figure;
    scatter(ISI_MEAN_WT_All,ISI_CV_WT_All,20,'ok'); hold on; scatter(ISI_MEAN_WT_All,ISI_CV_WT_All,'.k');
    scatter(mean(ISI_MEAN_WT_All(~isnan(ISI_MEAN_WT_All))),mean(ISI_CV_WT_All(~isnan(ISI_CV_WT_All))),50,'.b');
    scatter(ISI_MEAN_KO_All,ISI_CV_KO_All,20,'om'); hold on; scatter(ISI_MEAN_KO_All,ISI_CV_KO_All,'.m');
    scatter(mean(ISI_MEAN_KO_All(~isnan(ISI_MEAN_KO_All))),mean(ISI_CV_KO_All(~isnan(ISI_CV_KO_All))),50,'.r');
    xlabel('Mean ISI (ms)')
    ylabel('CV')
    
    %2 remove NaN data
    figure;
    scatter(valISI_MEAN_WT_A ,valISI_CV_WT_A ,20,'ok'); hold on;
    scatter(mean(valISI_MEAN_WT_A ),mean(valISI_CV_WT_A) ,50,'*b'); hold on;
    scatter(valISI_MEAN_KO_A ,valISI_CV_KO_A ,20,'om'); hold on;
    scatter(mean(valISI_MEAN_KO_A ),mean(valISI_CV_KO_A) ,50,'*r'); hold on;
    
    % remove cells that have only one spike and cause std = 0, hence CV = 0;
    figure;
    % scatter(valISI_MEAN_WT_A ,valISI_CV_WT_A ,20,'ok'); hold on;
    % scatter(mean(valISI_MEAN_WT_A ),mean(valISI_CV_WT_A) ,50,'xb'); hold on;
    % scatter(valISI_MEAN_KO_A ,valISI_CV_KO_A ,20,'om'); hold on;
    % scatter(mean(valISI_MEAN_KO_A ),mean(valISI_CV_KO_A) ,50,'xr'); hold on;
    
    scatter(valISI_MEAN_WT_A(valISI_CV_WT_A ~=0) ,valISI_CV_WT_A(valISI_CV_WT_A ~=0)  ,'.k'); hold on;
    scatter(mean(valISI_MEAN_WT_A(valISI_CV_WT_A ~=0) ),mean(valISI_CV_WT_A(valISI_CV_WT_A ~=0) ) ,100,'*b'); hold on;
    scatter(valISI_MEAN_KO_A(valISI_CV_KO_A ~=0) ,valISI_CV_KO_A(valISI_CV_KO_A ~=0) ,'.m'); hold on;
    scatter(mean(valISI_MEAN_KO_A(valISI_CV_KO_A ~=0) ),mean(valISI_CV_KO_A(valISI_CV_KO_A ~=0)) ,100,'*r'); hold on;
    % legend('WT all','Mean of all WT','KO all','Mean of all KO','WT nonzero','Mean of nonzero WT','KO nonzero','Mean of nonzero KO')
    legend('WT nonzero','Mean of nonzero WT','KO nonzero','Mean of nonzero KO')
    scatter(valISI_MEAN_WT_A(valISI_CV_WT_A ~=0) ,valISI_CV_WT_A(valISI_CV_WT_A ~=0),20  ,'ok'); hold on;
    scatter(valISI_MEAN_KO_A(valISI_CV_KO_A ~=0) ,valISI_CV_KO_A(valISI_CV_KO_A ~=0) ,20,'om'); hold on;
    
    xlabel('Mean ISI (ms)')
    ylabel('CV')
    
    mean_point_WT_All = [mean(valISI_MEAN_WT_A ) mean(valISI_CV_WT_A)];
    mean_point_KO_All = [mean(valISI_MEAN_KO_A ) mean(valISI_CV_KO_A)];
    meanDiff_All =  pdist([mean_point_WT_All; mean_point_KO_All ]);
    
    mean_point_nonzero_WT_All = [mean(valISI_MEAN_WT_A(valISI_CV_WT_A ~=0)) mean(valISI_CV_WT_A(valISI_CV_WT_A ~=0))];
    mean_point_nonzero_KO_All = [mean(valISI_MEAN_KO_A(valISI_CV_KO_A ~=0)) mean(valISI_CV_KO_A(valISI_CV_KO_A ~=0))];
    meanDiffnonZero_All =  pdist([mean_point_nonzero_WT_All; mean_point_nonzero_KO_All]);
    
    % title(['Different of mean = ' num2str(meanDiff_All) ', when consider only nonzero values mean = ' num2str(meanDiffnonZero_All)])
    title(['Different of mean = ' num2str(meanDiffnonZero_All)])
    
    nonZeroISI_MEAN_WT_A = valISI_MEAN_WT_A(valISI_CV_WT_A ~=0); % 52 data points
    nonZeroISI_CV_WT_A = valISI_CV_WT_A(valISI_CV_WT_A ~=0);
    
    nonZeroISI_MEAN_KO_A = valISI_MEAN_KO_A(valISI_CV_KO_A ~=0); % 33 data points
    nonZeroISI_CV_KO_A = valISI_CV_KO_A(valISI_CV_KO_A ~=0);
    
    %% Cumulative ISI
    figure;
    [f_WT,x_WT] = ecdf(ISI_combine_WT_All);
    [f_KO,x_KO] = ecdf(ISI_combine_KO_All);
    plot(x_WT,f_WT,'k'); hold on;
    plot(x_KO,f_KO,'m'); hold on;
    legend('WT','KO')
    ylabel('ratio'); xlabel('ISI(ms)')
    title('Cumulative distribution function of ISI')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% analysis with all values ()
    
    % %Just for Ref / delete later
    % valISI_MEAN_WT_A = ISI_MEAN_WT_All(~(isnan(ISI_MEAN_WT_All) | isinf(ISI_MEAN_WT_All) )); % 64 Data points
    % valISI_CV_WT_A = ISI_CV_WT_All(~(isnan(ISI_CV_WT_All) | isinf(ISI_CV_WT_All)));
    % valISI_MEAN_KO_A = ISI_MEAN_KO_All(~(isnan(ISI_MEAN_KO_All) | isinf(ISI_MEAN_KO_All) )); % 56 Data point
    % valISI_CV_KO_A = ISI_CV_KO_All(~(isnan(ISI_CV_KO_All) | isinf(ISI_CV_KO_All)));
    
    N_trial = 1000;
    randID = 1:length(valISI_MEAN_WT_A)+length(valISI_MEAN_KO_A);
    totalData = length(valISI_MEAN_WT_A)+length(valISI_MEAN_KO_A);
    CombinedData_Mean = [valISI_MEAN_WT_A valISI_MEAN_KO_A];
    CombinedData_CV = [valISI_CV_WT_A valISI_CV_KO_A];
    surrogateISI_MEANdiff = zeros(1,N_trial );
    % pick all cells but classified as WT or KO randomly  (Surrogate data)
    N1 = 30; %length(valISI_MEAN_WT_A);
    N2 = 20; %length(valISI_MEAN_KO_A);
    fshow = figure;
    scatter(mean(valISI_MEAN_WT_A ),mean(valISI_CV_WT_A) ,50,'*b'); hold on;
    scatter(mean(valISI_MEAN_KO_A ),mean(valISI_CV_KO_A) ,50,'*r'); hold on;
    xlim([0 80]); ylim([0 1.4]);
    for trial = 1: N_trial
        clf
        scatter(mean(valISI_MEAN_WT_A ),mean(valISI_CV_WT_A) ,50,'*b'); hold on;
        scatter(mean(valISI_MEAN_KO_A ),mean(valISI_CV_KO_A) ,50,'*r'); hold on;
        xlim([0 80]); ylim([0 1.4]);
        
        tmp = randperm(totalData);
        tmpDataGroup1 = [ mean(CombinedData_Mean(tmp(1:N1))) mean(CombinedData_CV(tmp(1:N1)))];
        tmpDataGroup2 = [ mean(CombinedData_Mean(tmp(N1+1:N1+N2))) mean(CombinedData_CV(tmp(N1+1:N1+N2)))];
        scatter(tmpDataGroup1(1), tmpDataGroup1(2),'xk');
        scatter(tmpDataGroup2(1), tmpDataGroup2(2),'xm');
        surrogateISI_MEANdiff(trial) = pdist([tmpDataGroup1 ; tmpDataGroup2;]);
        title([ '#' num2str(trial) '    ' num2str(surrogateISI_MEANdiff(trial)-meanDiff_All) ])
        pause(0.001);
        %     xlim([20 30])
    end
    sum((surrogateISI_MEANdiff > meanDiff_All))
    %% Nonzero
    N_trial = 1000;
    
    nonZerototalData = length(nonZeroISI_MEAN_WT_A)+length(nonZeroISI_MEAN_KO_A);
    nonZeroCombinedData_Mean = [nonZeroISI_MEAN_WT_A nonZeroISI_MEAN_KO_A];
    nonZeroCombinedData_CV = [nonZeroISI_CV_WT_A nonZeroISI_CV_KO_A];
    nonZerosurrogateISI_MEANdiff = zeros(1,N_trial );
    % pick all cells but classified as WT or KO randomly  (Surrogate data)
    N1 = 30; % length(nonZeroISI_MEAN_WT_A);
    N2 = 20; %length(nonZeroISI_MEAN_KO_A);
    fshow = figure;
    scatter(mean(nonZeroISI_MEAN_WT_A ),mean(nonZeroISI_CV_WT_A) ,50,'*b'); hold on;
    scatter(mean(nonZeroISI_MEAN_KO_A ),mean(nonZeroISI_CV_KO_A) ,50,'*r'); hold on;
    xlim([0 80]); ylim([0 1.4]);
    for trial = 1: N_trial
        clf
        scatter(mean(nonZeroISI_MEAN_WT_A ),mean(nonZeroISI_CV_WT_A) ,50,'*b'); hold on;
        scatter(mean(nonZeroISI_MEAN_KO_A ),mean(nonZeroISI_CV_KO_A) ,50,'*r'); hold on;
        xlim([0 80]); ylim([0 1.4]);
        
        tmp = randperm(nonZerototalData);
        tmpDataGroup1 = [ mean(nonZeroCombinedData_Mean(tmp(1:N1))) mean(nonZeroCombinedData_CV(tmp(1:N1)))];
        tmpDataGroup2 = [ mean(nonZeroCombinedData_Mean(tmp(N1+1:N1+N2))) mean(nonZeroCombinedData_CV(tmp(N1+1:N1+N2)))];
        scatter(tmpDataGroup1(1), tmpDataGroup1(2),'xk');
        scatter(tmpDataGroup2(1), tmpDataGroup2(2),'xm');
        nonZerosurrogateISI_MEANdiff(trial) = pdist([tmpDataGroup1 ; tmpDataGroup2;]);
        title([ '#' num2str(trial) '    ' num2str(nonZerosurrogateISI_MEANdiff(trial)-meanDiffnonZero_All) ])
        pause(0.001);
        %     xlim([20 30])
    end
    sum((nonZerosurrogateISI_MEANdiff > meanDiffnonZero_All))
    
    %% Z-score
    
    ISI_Mean_Combine = [nonZeroISI_MEAN_WT_A nonZeroISI_MEAN_KO_A];
    ISI_CV_Combine = [nonZeroISI_CV_WT_A  nonZeroISI_CV_KO_A ];
    Z_ISI_MEAN_WT = (nonZeroISI_MEAN_WT_A - mean(ISI_Mean_Combine))./std(ISI_Mean_Combine); % Note: There is built-in function zscore(X), for which z-scores = (X-MEAN(X)) ./ STD(X)
    Z_ISI_CV_WT =(nonZeroISI_CV_WT_A - mean(ISI_CV_Combine))./std(ISI_CV_Combine) ;
    Z_ISI_MEAN_KO = (nonZeroISI_MEAN_KO_A - mean(ISI_Mean_Combine))./std(ISI_Mean_Combine);
    Z_ISI_CV_KO = (nonZeroISI_CV_KO_A- mean(ISI_CV_Combine))./std(ISI_CV_Combine);
    
    Z_ISI_Mean_Combine = [Z_ISI_MEAN_WT Z_ISI_MEAN_KO];
    Z_ISI_CV_Combine = [Z_ISI_CV_WT Z_ISI_CV_KO];
    
    mean_point_Z_ISI_WT =  [mean(Z_ISI_MEAN_WT) mean(Z_ISI_CV_WT)];
    mean_point_Z_ISI_KO =  [mean(Z_ISI_MEAN_KO) mean(Z_ISI_CV_KO)];
    meanDiff_Z_ISI = pdist([mean_point_Z_ISI_WT; mean_point_Z_ISI_KO]);
    
    
    figure;
    subplot(121)
    scatter(nonZeroISI_MEAN_WT_A,nonZeroISI_CV_WT_A,'ok' ); hold on; % Normal
    scatter(nonZeroISI_MEAN_KO_A,nonZeroISI_CV_KO_A,'or' ); % Normal
    subplot(122)
    scatter(Z_ISI_MEAN_WT,Z_ISI_CV_WT,'ok' );hold on; % Checking ---> the z-score vals are correct.
    scatter(Z_ISI_MEAN_KO,Z_ISI_CV_KO,'or' ); %
    scatter( mean(Z_ISI_MEAN_WT),mean(Z_ISI_CV_WT),'xk')
    scatter( mean(Z_ISI_MEAN_KO),mean(Z_ISI_CV_KO),'xr')
    
    
    N_trial = 1000;
    
    totalData_Z_ISI = length(Z_ISI_Mean_Combine);
    
    surrogateISI_Z_MEANdiff = zeros(1,N_trial);
    surrogateISI_Z_MEANval_WT = zeros(N_trial,2);
    surrogateISI_Z_MEANval_KO = zeros(N_trial,2);
    % pick all cells but classified as WT or KO randomly  (Surrogate data)
    N1 = length(Z_ISI_MEAN_WT);
    N2 = length(Z_ISI_MEAN_KO);
    fshow = figure;
    scatter( mean(Z_ISI_MEAN_WT),mean(Z_ISI_CV_WT),50,'*k'); hold on;
    scatter( mean(Z_ISI_MEAN_KO),mean(Z_ISI_CV_KO),50,'*r'); hold on;
    xlim([-2 3]); ylim([-2 3]);
    for trial = 1: N_trial
        clf
        scatter( mean(Z_ISI_MEAN_WT),mean(Z_ISI_CV_WT),50,'*k'); hold on;
        scatter( mean(Z_ISI_MEAN_KO),mean(Z_ISI_CV_KO),50,'*r'); hold on;
        xlim([-2 3]); ylim([-2 3]);
        
        tmp = randperm(totalData_Z_ISI);
        tmpDataGroup1 = [ mean(Z_ISI_Mean_Combine(tmp(1:N1))) mean(Z_ISI_CV_Combine(tmp(1:N1)))];
        tmpDataGroup2 = [ mean(Z_ISI_Mean_Combine(tmp(N1+1:N1+N2))) mean(Z_ISI_CV_Combine(tmp(N1+1:N1+N2)))];
        scatter(tmpDataGroup1(1), tmpDataGroup1(2),'xk');
        scatter(tmpDataGroup2(1), tmpDataGroup2(2),'xm');
        surrogateISI_Z_MEANval_WT(trial,:) = tmpDataGroup1;
        surrogateISI_Z_MEANval_KO(trial,:) = tmpDataGroup2;
        surrogateISI_Z_MEANdiff(trial) = pdist([tmpDataGroup1 ; tmpDataGroup2;]);
        title([ '#' num2str(trial) '    ' num2str(surrogateISI_Z_MEANdiff(trial)-meanDiff_Z_ISI) ])
        pause(0.00001);
        %     xlim([20 30])
    end
    eventCnt = sum((surrogateISI_Z_MEANdiff > meanDiff_Z_ISI));
    pval = eventCnt / N_trial;
    
    fz = figure; set(fz, 'position', [  508   214   806   628]);
    scatter(Z_ISI_MEAN_WT,Z_ISI_CV_WT,20,'ok' );hold on;
    scatter(Z_ISI_MEAN_KO,testZ_ISI_CV_KO,20,'om');
    scatter( surrogateISI_Z_MEANval_WT(:,1),surrogateISI_Z_MEANval_WT(:,2),'xk'); hold on;
    scatter( surrogateISI_Z_MEANval_KO(:,1),surrogateISI_Z_MEANval_KO(:,2),'xm'); hold on;
    scatter( mean(Z_ISI_MEAN_WT),mean(Z_ISI_CV_WT),80,'*b','LineWidth',2); hold on;
    scatter( mean(Z_ISI_MEAN_KO),mean(Z_ISI_CV_KO),80,'*r','LineWidth',2); hold on;
    legend('data WT', 'data KO', 'controls'' average WT','controls'' average KO','data average WT','data average KO')
    scatter(Z_ISI_MEAN_WT,Z_ISI_CV_WT,'.k' );hold on;
    scatter(Z_ISI_MEAN_KO,testZ_ISI_CV_KO,'.m');
    xlim([-2 3]); ylim([-2 3]);
    xlabel('mean ISI (Z-score)')
    ylabel('CV (Z-score)')
    title(['Mean value of WT and KO in mean-CV space of ISI is separated by ' num2str(meanDiff_Z_ISI) ' (z-score) with p-value = ' num2str(pval) ' (' num2str(N_trial) ' controlled samples)'])
    
    figure; hist(surrogateISI_Z_MEANdiff,100); [hm, xn] = hist(surrogateISI_Z_MEANdiff,100);
    line([meanDiff_Z_ISI meanDiff_Z_ISI],[1 max(hm)],'Color','r')
    legend('Control Mean different', 'Real Data Mean Different')
    title(['Distribution of mean different of control(surrogate data) and real data ( ' num2str(meanDiff_Z_ISI) ')']);
    %% spike correlation of the cells during rebounding and during basal activity
    % InjStartT = 500;
    % InjStopT = 1000;
    % burstRange = 500;
    % cutTime = avgFR_CUTTIME;
    
    % corr_light_on = calc_corr_interval(MUAlight, LFPlight, MUAlightT, [0 10]);
    %
    % function yy = calc_corr_interval(mua0, lfp0, tt, t_interval)
    %
    % mua = mua0(:, t_interval(1) <= tt & tt < t_interval(2));
    % lfp = lfp0(:, t_interval(1) <= tt & tt < t_interval(2));
    %
    % yy = [];
    %
    % maxlag = round(100/1000/diff(tt(1:2)));     % [ms]
    %
    % for ii=1:size(mua, 1)
    %     [cc0, lags0] = xcorr(mua(ii, :), lfp(ii, :), maxlag, 'coeff');
    % %     [cc0, lags0] = xcorr(mua(ii, :), lfp(ii, :), 'coeff');
    %
    %     cc0 = cc0(lags0<=0);
    %     lags0 = lags0(lags0<=0);
    %
    % %     cind = (abs(cc0(:)) == max(abs(cc0(:))));
    %     cind = (cc0(:) == max(cc0(:)));
    %
    %     cc = cc0(cind);
    %     lags = lags0(cind);
    %
    %     yy = [yy; cc lags*diff(tt(1:2))];
    % end
    
end
%%
if(SAVE_BASAL_ACT)
    Date = clock();
    save([num2str(Date(1))  '_' num2str(Date(2)) '_' num2str(Date(3)) '_' 'VLBasalAct_' FNAME_SIM '.mat'], 'Basal_Act');
    disp(['Save Basal_Act to  ' num2str(Date(1))  '_' num2str(Date(2)) '_' num2str(Date(3)) '_' 'VLBasalAct_' FNAME_SIM '.mat']);% save('140715VL_Basal_Act.mat', 'Basal_Act');
end
end