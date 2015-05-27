
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

dirLoc = 'Model1_MinParam150403_newExpdata/' ;

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
dirFig = 'DeleteThis'; %  'Fig_leastW_100s/';
mkdir([dirLoc dirFig])
InFR_CODENAME = 'PoisInputFr150403_leastW_newExp'; 

NoiseMEAN_WTKO = [0; 0;]; %[0.0291; 0.0075;]; %[0.0295; 0.016]; %%%%%%%%%%%%%%%%%%%%%%%%%%% [ WT; KO;] for 10 Hz -> [0.0295; 0.016]; for 5 Hz -> [0.0291; 0.0075;]; 5Hz 0.0318 0.3  for WT 5Hz (2) -> [0.0318; 0.0075;]
NoiseSIG_WTKO  =  [0.3; 0.3;]; %[0.8; 0.8;];  %[0.24; 0.3];  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  [WT; KO;] for 10 Hz -> [0.4; 0.5]; for 5 Hz -> [0.24; 0.3];  5Hz 0.0318 0.3  for WT 5Hz (2) -> [0.3; 0.3;]
VARIE_MEAN_InGauss = 1;
NoiseMEANsigma_WTKO = [0; 0;]; %%%%% [WT; KO;]
TSTOP = 3000;

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
        %PDfewBurst_GPmVLmd1_KO_InGauss0.2_IGmean0_IGmeanSig0_W0.029_Specifi                                                         edPoisSpk_sig0.00Hz_T100000_trial5

        InGauss_STDEV = 0.2; %0.2;, 0.3
        NoiseMEAN = 0;
        IGmeanSig = 0;
        W_Weight = 0.029;
        W_sig = 0;
        PoisInputFr = 0;
        PoisInputSig = 0;
        TSTOP = 100000;
        %         Rand123Seed_ID1_PDfewBurst_GPmVLmd1_WT_InGauss0.13_IGmean-0.2_IGmeanSig0_W0.001_SpecifiedPoisSpk_sig0.00Hz_T2000_trial4
        %                 GPmLightDur = GPmLightDur_LST(ld_ii);
        %                 PhotoStop = PhotoInjT+GPmLightDur;
        %                 PhotoStop = PhotoStop;
        
        %                 GPmLight = LightAmp_LST(la_ii);
        %                 GPm_w_mn = GPmVLw_mean_LST(m_ii);
        %                 GPm_w_sg = GPmVLw_sig_LST(s_ii);
        
        
        Name_postfix = [coreFileName '_' cTxt '_InGauss' num2str(InGauss_STDEV) '_IGmean' num2str(NoiseMEAN) '_IGmeanSig' num2str(IGmeanSig) '_W' num2str(W_Weight) '_SpecifiedPoisSpk_sig' num2str(PoisInputSig)  '.00Hz_T' num2str(TSTOP) '_trial' num2str(TRIAL_NO) ];
%         Name_postfix = [coreFileName '_' cTxt  '_IGmean' num2str(NoiseMEAN) '_IGmeanSig' num2str(IGmeanSig) '_W' num2str(W_Weight) '_Wsig' num2str(W_sig) '_SpecifiedPoisSpk_sig' num2str(PoisInputSig)  '.00Hz_T' num2str(TSTOP) '_trial' num2str(TRIAL_NO) ];
%    PDfewBurst_GPmVLmd1_KO_InGauss0.2_IGmean0_IGmeanSig0_W0.029_Specifi                                                         edPoisSpk_sig0.00Hz_T100000_trial5

        
        
        
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
        
        Nsample = 18;
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
getExpData ; % WTWT_MS, KOKO_MS, WTWT, KOKO

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
% yyy = 64.78.* exp(0.0147.*xx) - 64.67 ; %parameter search result
% plot(xx, yyy,'k:','linewidth',2)
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
yfit = 64.78.* exp(0.0147.*xx) - 64.67 ; %parameter search result


ff = figure; set(ff,'position',[ 1000         827         647         511]);  set(gcf, 'PaperPositionMode','Auto');
plot(xx,yfit,'k.-');
xlabel('Output Firing rate'); ylabel('Input Firing rate');
hold on;
LEG{1} = 'Fitting equation';
for TRIAL_NO = 1:NUM_TRIAL
    title(num2str(TRIAL_NO))
    WT_in = dlmread([InFR_CODENAME '_WT_' num2str(TRIAL_NO) '.txt']); WT_in = WT_in(2:end);
    tmp = ACT_Record{TRIAL_NO}.WT.All.fr_data;
    scatter(tmp,WT_in,'MarkerEdgeColor',[0 0 ccc(TRIAL_NO)]);
    LEG{2*TRIAL_NO} = ['WT, trial#' num2str(TRIAL_NO)];
    k = waitforbuttonpress;
    
    KO_in = dlmread([InFR_CODENAME '_KO_' num2str(TRIAL_NO) '.txt']); KO_in = KO_in(2:end);
    tmp = ACT_Record{TRIAL_NO}.KO.All.fr_data;
    scatter(tmp,KO_in,'MarkerEdgeColor',[ccc(TRIAL_NO) 0 0]);
    LEG{2*TRIAL_NO+1} = ['KO, trial#' num2str(TRIAL_NO)];
    k = waitforbuttonpress;
end
title('Output-Input Relationship')
legend(LEG,'location','best')
 k = waitforbuttonpress;

saveas( ff, [ dirLoc dirFig 'OutInPlot'  '.jpg'], 'jpg')
saveas( ff, [ dirLoc dirFig 'OutInPlot'  '.fig'], 'fig')   

%%
% 
%   
%  Least W : 100s
% Fitting Curve for experimental data:  y = -0.018795x^2 +1.3976x +8.8904
% Fitting Curve for Trial#1
% WT:  y = -0.024444x^2 +1.4128x +6.9444
% KO:  y = -0.027293x^2 +1.5041x +6.405
% ALL: y = -0.025823x^2 +1.458x +6.6828
% Fitting Curve for Trial#2
% WT:  y = -0.024925x^2 +1.4178x +6.9788
% KO:  y = -0.026941x^2 +1.4895x +6.5069
% ALL: y = -0.025993x^2 +1.4557x +6.7409
% Fitting Curve for Trial#3
% WT:  y = -0.024658x^2 +1.4177x +6.9388
% KO:  y = -0.026855x^2 +1.4861x +6.5498
% ALL: y = -0.025695x^2 +1.4509x +6.7525
% Fitting Curve for Trial#4
% WT:  y = -0.024805x^2 +1.4181x +6.9464
% KO:  y = -0.027578x^2 +1.5129x +6.3327
% ALL: y = -0.02623x^2 +1.4672x +6.6412
% Fitting Curve for Trial#5
% WT:  y = -0.023639x^2 +1.3869x +7.1187
% KO:  y = -0.025988x^2 +1.4578x +6.741
% ALL: y = -0.024731x^2 +1.4209x +6.9388
% Trial#1
% WT: mean = 10.7213 , std = 7.0942
% KO: mean = 12.0213 , std = 6.8073
% Trial#2
% WT: mean = 10.4892 , std = 6.8102
% KO: mean = 11.9884 , std = 6.66
% Trial#3
% WT: mean = 10.4175 , std = 6.8868
% KO: mean = 12.0728 , std = 6.8148
% Trial#4
% WT: mean = 10.2943 , std = 6.8599
% KO: mean = 11.8143 , std = 6.7192
% Trial#5
% WT: mean = 10.6923 , std = 6.8911
% KO: mean = 12.1863 , std = 6.7355
% 
