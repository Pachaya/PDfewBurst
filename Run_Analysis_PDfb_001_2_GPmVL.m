
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
SAVE_NEURON_ACT = 1;
FNAME_SIM = 'GPmVL' ; %'PreLim_INOISE0.2_InSPK50_EE1_EI2_IE20_II12_mult125_w0.0006';
plotSampleV = 0;
SPECIFIED_RATIO = 0; %TestSim_WT_InGauss0.2_W0.0005_50.00Hz_rEE0_rEI1_rIE2_rII1_Wmult10
Basal_Act_set  = [];

NoiseMEAN = 0;

SAVE_FIG = 1;

dirLoc = 'PulseTest/' ; %'TestSim_BasalAct/';  %'VL_LocalConn/';
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


for p1_ii = 1 : length(PARAM1)
    for p2_ii = 1 : length(PARAM2)
        for p3_ii = 1 : length(PARAM3)
            la_ii = p1_ii; m_ii = 1; ld_ii = p2_ii; s_ii = p3_ii;
            
            for cell_type = 1 : 2
                
                if (cell_type == 1)
                    cTxt = 'WT';
                elseif (cell_type == 2)
                    cTxt = 'KO';
                end
                
                coreFileName = 'PDfewBurst_GPmVLmd1' ;
                
                InGauss_STDEV = 0.25; %0.2;, 0.3
                NoiseMEAN = 0;
                IGmeanSig = 0;
                W_Weight = 0.001;
                PoisInputFr = 0;
                TSTOP = 3000;
                
                GPmLightDur = GPmLightDur_LST(ld_ii);                
                PhotoStop = PhotoInjT+GPmLightDur;
                PhotoStop = PhotoStop;


            
                GPmLight = LightAmp_LST(la_ii);
                GPm_w_mn = GPmVLw_mean_LST(m_ii);
                GPm_w_sg = GPmVLw_sig_LST(s_ii);
                
                
                Name_postfix = [coreFileName '_' cTxt '_' 'GPmInput_Amp' num2str(GPmLight) '_Dur' num2str(GPmLightDur) '_GPmVLw_m' num2str(GPm_w_mn) '_sig' num2str(GPm_w_sg) ...
                    '_InGauss' num2str(InGauss_STDEV) '_IGmean' num2str(NoiseMEAN) '_IGmeanSig' num2str(IGmeanSig) '_W' num2str(W_Weight) '_' num2str(PoisInputFr) '.00Hz_T' num2str(TSTOP) ];
                
                % PDfewBurst_GPmVLmd1_KO_GPmInput_Amp0.4_Dur1000_GPmVLw_m0.1_sig0.03_InGauss0.2_IGmean0_IGmeanSig0_W0.001_0.00Hz_T3000
                
                
                
                disp('==================================================================================================')
                disp(Name_postfix)
                if (cell_type == 1)
                    cTxt = 'WT';
                else
                    cTxt = 'KO';
                end
                figNameCode = [ saveTxt1 num2str(PARAM1(p1_ii))  '_' saveTxt2 num2str(PARAM2(p2_ii)) '_' cTxt];
                
%                 CheckFileExist( dirLoc, Name_postfix  )
                
                VL_PhotoactivationAll
                
%                 if(0)
                
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
            DoSampleTtest = 0;
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
            
            ACT_Record{p1_ii,p2_ii,p3_ii } = Basal_Act;
            
        end
    end
end


%% 
if (SAVE_NEURON_ACT)
    save([dirLoc dirFig 'Result_ACT_Record_' coreFileName '.mat'], 'ACT_Record','-v7.3');
end
%% 


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


    