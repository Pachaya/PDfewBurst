
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

SAVE_WORKSPACE = 0;

DoTTest = 1; DoSampleTtest = 1;
SAVE_BASAL_ACT = 1;
plotSampleV = 1;
SPECIFIED_RATIO = 0; %TestSim_WT_InGauss0.2_W0.0005_50.00Hz_rEE0_rEI1_rIE2_rII1_Wmult10


NoiseMEAN = 0;

SAVE_FIG = 1; Close_Fig_aftr_save = 1;

dirLoc = 'Model1_MinParam150403_fullsim/' ;
rEE_rEI_rIE_rII = [1 1 1 1];
w_MULT = 125;
avgFR_CUTTIME = 500;

pulseHz = 1; N_pulse = 1;
Input_I_amp_list = [0];
Input_I_dur_list = [0];
Input_I_amp = -0.3; %[-1; -2; -3; -4; -5;];% [-1; -5]; %
Input_I_dur = 500; %[10; 20; 30; 40; 50;];
ReboundPeakAmp_WTlist = zeros(length(Input_I_amp_list),length(Input_I_dur_list));
ReboundPeakDel_WTlist = zeros(length(Input_I_amp_list),length(Input_I_dur_list));
ReboundPeakAmp_KOlist = zeros(length(Input_I_amp_list),length(Input_I_dur_list));
ReboundPeakDel_KOlist = zeros(length(Input_I_amp_list),length(Input_I_dur_list));
Save_BasalAct =  cell(length(Input_I_amp_list),length(Input_I_dur_list));
InjStartT = 1500;
InjStopT = InjStartT+Input_I_dur;
reboundDelay = 0;
burstRange = 100;
cutTime = avgFR_CUTTIME;

ADD_I_to_M1 = 0;
ADD_I_to_VL = 0;

FIG_ALL = 0;
dirFig = 'Fig_fullSim_GPmVL/';
mkdir([dirLoc dirFig])
TSTOP = 4000;

rTC_LST = [120; 150];
wmTC_LST = [50];



LightAmp_LST = [0.4];
GPmVLw_mean_LST = [0.2 0.3 0.4 0.5];    
GPmVLw_sig_LST =[0 0.05 0.1]; 

% Candidates Model 1 :
% I amp = 0.3, mean 0.06, sigma = 0.01, 0.005 ?0.03
% I amp = 0.5, mean 0.04, sigma = 0.03, 0.01 ?0.03
%PDfewBurst_GPmVLmd1_rTC120_wmTC10_WT_GPmInput_Amp0.3_Dur1000_GPmVLw_m0.06_sig0.01_InGauss0.2_IGmean-0.15_IGmeanSig0_W0.0015_SpecifiedPoisSpk_sig0.00Hz_T4000_trial3
NUM_TRIAL = 5;
ncells = 1150;
TRIAL_LST = 1 : NUM_TRIAL;

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
PARAM5 = GPmVLw_sig_LST;
lblTxt5 = 'GPm-VL weight sigma';
saveTxt5 = 'GPmVLw_sig';
titleTxt5 = 'W_s_i_g';



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

CUTTIME = 500;
PhotoInjT = 1500;
PhotoStop = 2500;
DelayT = 0;
BurstRange = 100;

%CheckSimulationByName_Combine_GpmVL_VLM1

for p1_ii = 1 : length(PARAM1)
    for p2_ii = 1 : length(PARAM2)
        for p3_ii = 1 : length(PARAM3)
            for p4_ii = 1 : length(PARAM4)
                for p5_ii = 1 : length(PARAM5)
                    
                    r_ii = p1_ii; wm_ii = p2_ii;
                    la_ii = p3_ii;  m_ii = p4_ii; s_ii = p5_ii;
                    TRIAL_NO = 1;
                    
                    for cell_type = 1 : 2
                        
                        if (cell_type == 1)
                            cTxt = 'WT';
                        elseif (cell_type == 2)
                            cTxt = 'KO';
                        end
                        
                        
                        coreFileName = 'PDfewBurst_GPmVLmd1' ;
                        
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
                        
                        
                        Simulation_Code = [coreFileName '_rTC' num2str(rTC) '_wmTC' num2str(wmTC) '_' cTxt '_' 'GPmInput_Amp' num2str(GPmLight) '_Dur' num2str(GPmLightDur) '_GPmVLw_m' num2str(GPm_w_mn) '_sig' num2str(GPm_w_sg) ...
                            '_InGauss' num2str(InGauss_STDEV) '_IGmean' num2str(NoiseMEAN) '_IGmeanSig' num2str(IGmeanSig) '_W' num2str(W_Weight) '_SpecifiedPoisSpk_sig0.00Hz_T' num2str(TSTOP) '_trial' num2str(TRIAL_NO)];
                        % PDfewBurst_GPmVLmd1_0del_rTC120_wmTC4_KO_GPmInput_Amp0.3_Dur1000_GPmVLw_m0.04_sig0_InGauss0.2_IGmean-0.15_IGmeanSig0_W0.0012_SpecifiedPoisSpk_sig0.00Hz_T4000_trial1
                        disp('==================================================================================================')
                        disp(Simulation_Code)
                        disp('######  Download M1 ')
                        %M1
                        Name_postfix = [ 'M1_' Simulation_Code];
                        getCenterPart =1;
                        BaselineAct % Get all cell activity
                        
                        tmp.E = E;      tmp.I = I;  tmp.All = All; %tmp.VL = VL; tmp.M1 = M1;
                        if (cell_type == 1)
                            M1.WT = tmp;
                            cTxt = 'WT';
                        elseif (cell_type == 2)
                            M1.KO = tmp;
                            cTxt = 'KO';
                        end
                        
                        
                        figNameCode = [ saveTxt1 num2str(PARAM1(p1_ii))  '_' saveTxt2 num2str(PARAM2(p2_ii)) '_' cTxt];
                        
                        %                 CheckFileExist( dirLoc, Name_postfix  )
                        disp('######  Download VL ')
                        Name_postfix = [ Simulation_Code];
                        VL_PhotoactivationAll
                        
                        %                 if(0)
                        
                        tmp.E = E;      tmp.I = I;  tmp.All = All; %tmp.VL = VL; tmp.M1 = M1;
                        clear E I All somaVall somaVallE somaVallI
                        if (cell_type == 1)
                            VL.WT = tmp;
                            cTxt = 'WT'; Name_postfix_WT = Name_postfix;
                            clear tmp
                        elseif (cell_type == 2)
                            VL.KO = tmp;
                            cTxt = 'KO'; Name_postfix_KO = Name_postfix;
                            clear tmp
                        end
                        
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
                            if (Close_Fig_aftr_save)
                                close(fgSmpl);
                            end
                        end
                        
                        disp('Finish plot sample membrain voltage');
                        clear somaVall
                        toc(tt)
                    end
                    
                    
                    
                    Basal_Act.VL = VL;
                    Basal_Act.M1 = M1;
                    clear VL
                    clear M1
                    
                    if (DoTTest)  %for Baseline activity
                        % two-tailed t-test % during normal
                        
                        % disp('PoisSPk')
                        % [hP,pP] = ttest(WT.PoisSpk.fr_data, KO.PoisSpk.fr_data)
                        disp('--------------------------------------------------------------------------------------------------')
                        disp('Two-tailed t-test')
                        tmpWTfr = Basal_Act.VL.WT.All.fr_data;
                        tmpKOfr = Basal_Act.VL.KO.All.fr_data;
                        
                        [hE,pE] = ttest2(tmpWTfr, tmpKOfr, 'Tail','both');
                        disp([ ' p-value = ' num2str(pE)])
                        if(hE)
                            disp('H0 rejected : average normal firing rate of E cells in WT significantly different from KO (p < 0.05)')
                        else
                            disp('Do not rejected H0 : average normal firing rate of E cells in WT does not significantly different from KO (p > 0.05)')
                        end
                        Basal_Act.ttest2.hE = hE;
                        Basal_Act.ttest2.pE = pE;
                        
                        
                    end
                    %% Sample T-Test for Baseline activity
                    
                    if(DoSampleTtest)
                        tmpWTfr = Basal_Act.VL.WT.All.fr_data;
                        tmpKOfr = Basal_Act.VL.KO.All.fr_data;
                        
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
                            [ht,pt] = ttest2(tmpWTfr(sample), tmpKOfr(sample));
                            [h,p] = ranksum(tmpWTfr(sample), tmpKOfr(sample));
                            
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
                    
                    ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii } = Basal_Act;
                    clear Basal_Act
                end
            end
        end
    end
end

%%

if(SAVE_BASAL_ACT)
    save([dirLoc dirFig 'Saved_Activity_result.mat' ], 'ACT_Record','PARAMETERS','-v7.3');
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    GPM - VL Connnection
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rasterplot of VL  and GPm-VL chance of burst
SAVE_RASTER_FIG =1;
% GPm VL
CntBurstSpkWT = zeros(ACT_Rec_size);
CntBurstSpkKO = zeros(ACT_Rec_size);
ChanceOfBurstingWT = zeros(ACT_Rec_size);
ChanceOfBurstingKO = zeros(ACT_Rec_size);
TimeofFirstReboundSpk = zeros(ACT_Rec_size);
p1_ii = 1; p2_ii=1;
for p3_ii = 1 : length(PARAM3)
    fLA1 = figure; set(fLA1,'position',[302          49        1290         948]) ;  set(gcf,'PaperPositionMode','auto')
    fLA2 = figure; set(fLA2,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
    fB1 = figure; set(fB1,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
    fB2 = figure; set(fB2,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
    fRB = figure; set(fRB,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
    fBspk = figure; set(fBspk,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
    cnt = 0;
    nR = length(PARAM4); nC = length(PARAM5);
    for p4_ii = 1 : length(PARAM4)
        for p5_ii = 1 : length(PARAM5)
            la_ii = p3_ii; m_ii = p4_ii; s_ii = p5_ii;
            
            
            GPmLight = LightAmp_LST(la_ii);
            GPm_w_mn = GPmVLw_mean_LST(m_ii);
            GPm_w_sg = GPmVLw_sig_LST(s_ii);
            %               TRIAL_NO = p1_ii;
            %                 GPmLight = 0.3;
            %                 GPm_w_mn = 0.06; GPm_w_sg = 0.01;
            BasalAct = ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii}.VL;
            
            simTxt = ['VL ' titleTxt3 ' = ' num2str(PARAM3(p3_ii)) ', ' titleTxt4 ' = ' num2str(PARAM4(p4_ii)) ', ' titleTxt5 ' = ' num2str(PARAM5(p5_ii))];
            
            spkBin = BasalAct.WT.All.spktrain;
            baseline_fr_WT = get_avg_baseline(spkBin, cutTime, InjStartT);
            [freq_all_WT,spkTime_all_WT, fg_handle_WT ] = raster_from_spkbin_BurstRange( spkBin,PhotoStop, Tstop,PhotoInjT, PhotoStop, DelayT, BurstRange, [ simTxt ' : WT']);
            set(fg_handle_WT, 'position',[  449   450   791   528])
            
            spkBin = BasalAct.KO.All.spktrain;
            baseline_fr_KO = get_avg_baseline(spkBin, cutTime, InjStartT);
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
            
            WT_expBaselineSpk = baseline_fr_WT/1000*burstRange;
            WT_burstSpk = BasalAct.WT.All.BurstSpk; % BurstSpk  = Number of spike during burstRange
            CntBurstSpkWT(p1_ii,p2_ii,p3_ii,p1_ii,p2_ii,p3_ii) = mean(WT_burstSpk -WT_expBaselineSpk); % Number of spike during burstRange - expected number of spike from baseline activity
            
            KO_expBaselineSpk = baseline_fr_WT/1000*burstRange;
            KO_burstSpk = BasalAct.KO.All.BurstSpk; % BurstSpk  = Number of spike during burstRange
            CntBurstSpkKO(p1_ii,p2_ii,p3_ii,p4_ii,p5_ii) = mean(KO_burstSpk - KO_expBaselineSpk); % Number of spike during burstRange - expected number of spike from baseline activity
            
            ChanceOfBurstingWT(p1_ii,p2_ii,p3_ii,p4_ii,p5_ii) =  sum( round(WT_burstSpk -WT_expBaselineSpk) > 0)/(nE+nI);
            ChanceOfBurstingKO(p1_ii,p2_ii,p3_ii,p4_ii,p5_ii) =  sum( round(KO_burstSpk -KO_expBaselineSpk) > 0)/(nE+nI);
            
            cnt = cnt+1;
            figure(fLA1);
            subplot( nR,nC,cnt)
            [nb,xb]=hist(BasalAct.WT.All.fr_AfterLightOff,10); bh=bar(xb,nb); hold on;
            set(bh,'facecolor',[0 0 0]);
            xlim([0 max(xb)+1])
            title([ 'mean = ' num2str(GPm_w_mn) ', sig = ' num2str(GPm_w_sg)])
            figure(fLA2);
            subplot( nR,nC,cnt)
            [nb,xb]=hist(BasalAct.KO.All.fr_AfterLightOff,10); bh=bar(xb,nb); hold on;
            set(bh,'facecolor',[1 0 0]);
            xlim([0 max(xb)+1])
            title([ 'mean = ' num2str(GPm_w_mn) ', sig = ' num2str(GPm_w_sg)])
            
            figure(fB1);
            subplot( nR,nC,cnt)
            [nb,xb]=hist(BasalAct.WT.All.BurstSpk); bh=bar(xb,nb); hold on;
            set(bh,'facecolor',[0 0 0]);
            xlim([0 max(xb)+1])
            title([ 'mean = ' num2str(GPm_w_mn) ', sig = ' num2str(GPm_w_sg)])
            figure(fB2);
            subplot( nR,nC,cnt)
            [nb,xb]=hist(BasalAct.KO.All.BurstSpk,10); bh=bar(xb,nb); hold on;
            set(bh,'facecolor',[1 0 0]);
            xlim([0 max(xb)+1])
            title([ 'mean = ' num2str(GPm_w_mn) ', sig = ' num2str(GPm_w_sg)])
            
            
            % Distribution of first rebound spike timing 
            tmpRebound_WT = zeros(ncells,1); tmpRebound_KO = zeros(ncells,1);
            for id = 1 : ncells
                tmp = find(BasalAct.WT.All.BurstSpkTrain(id,:),1,'first');
                if ~isempty(tmp) 
                    tmpRebound_WT(id) = tmp;                  
                end
                
                tmp = find(BasalAct.KO.All.BurstSpkTrain(id,:),1,'first');
                if ~isempty(tmp) 
                    tmpRebound_KO(id) = tmp;                  
                end                
            end      
             tmpRebound_WT = tmpRebound_WT(tmpRebound_WT ~= 0);
             tmpRebound_KO = tmpRebound_KO(tmpRebound_KO ~= 0);
            figure(fRB); 
            subplot(nR,nC,cnt);
            [N,X] = hist(tmpRebound_WT,15);
            bar(X,N/ncells,'FaceColor','k','EdgeColor','w'); hold on; 
            [N,X] = hist(tmpRebound_KO,15);
            bar(X,N/ncells,'FaceColor','r','EdgeColor','w'); hold on; 
            title([ titleTxt4 ' = ' num2str(PARAM4(p4_ii)) ', ' titleTxt5 ' = ' num2str(PARAM5(p5_ii))])
            
            % Distribution of spikes in bursting range 
           
             tmpRebound_WT = sum(BasalAct.WT.All.BurstSpkTrain);
             tmpRebound_KO = sum(BasalAct.KO.All.BurstSpkTrain);
             xx = 1:1:size(BasalAct.WT.All.BurstSpkTrain,2);
            figure(fBspk); 
            subplot(length(PARAM2), length(PARAM3),cnt);
            bar(xx,tmpRebound_WT/ncells*100,'k'); hold on; bar(xx,tmpRebound_KO/ncells*100,'FaceColor','r','EdgeColor','r');
            ylabel('% Neurons'); xlabel('Time(ms)')
            title([ titleTxt4 ' = ' num2str(PARAM4(p4_ii)) ', ' titleTxt5 ' = ' num2str(PARAM5(p5_ii))])
            
            clear BasalAct spkBin
        end
    end
    figure(fLA1); suptitle(['Avg Fr after light off WT: ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]);
    figure(fB1); suptitle(['Number of Burst spike WT: ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]);
    figure(fLA2); suptitle(['Avg Fr after light off KO: ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]);
    figure(fB2);  suptitle(['Number of Burst spike KO: ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]);
    figure(fRB); suptitle(['Latency of first Rebound spike : ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]);  % Note equal to the latency of peak firing rate --> latency of peak firing rate calculate from instantaneous firing rate
    figure(fBspk); suptitle(['Percent of neuron that fire during bursting period : ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]);  % Note equal to the latency of peak firing rate --> latency of peak firing rate calculate from instantaneous firing rate
    
    if (SAVE_FIG)
        tmpTxt = get_Parameters_saveText(PARAMETERS, [3], [ p3_ii ]);
        ffig = [ dirLoc dirFig 'AvgFrLightOff' tmpTxt];
        saveas( fLA1, [ffig '_WT.jpg'], 'jpg')
        saveas( fLA1, [ffig '_WT.fig'], 'fig')
        saveas( fLA2, [ffig '_KO.jpg'], 'jpg')
        saveas( fLA2, [ffig '_KO.fig'], 'fig')        
        ffig = [ dirLoc dirFig 'NumBurstSpk'  tmpTxt];
        saveas( fB1, [ffig '_WT.jpg'], 'jpg')
        saveas( fB1, [ffig '_WT.fig'], 'fig')
        saveas( fB2, [ffig '_KO.jpg'], 'jpg')
        saveas( fB2, [ffig '_KO.fig'], 'fig')
        ffig = [ dirLoc dirFig 'FirstReboundSpike'  tmpTxt];
        saveas(fRB, [ffig '.fig'] , 'fig');
        saveas(fRB, [ffig '.jpg'] , 'jpg');
        ffig = [ dirLoc dirFig 'AllBurstSpike'  tmpTxt]; % Better minus by expected baseline activity
        saveas(fBspk, [ffig '.fig'] , 'fig');
        saveas(fBspk, [ffig '.jpg'] , 'jpg');
        
        if (Close_Fig_aftr_save)
            close(fLA1); close(fLA2); close(fB1); close(fB2); close(fRB); close(fBspk);
        end
        
    end
end

cnt =0;
p1_ii = 1; p2_ii=1;
for p3_ii = 1 : length(PARAM3)
    Bfg = figure;  set(Bfg,'position',[680   283   712   695]); set(gcf,'PaperPositionMode','auto');
    Bfg2 = figure;  set(Bfg2,'position',[680   283   712   695]);  set(gcf,'PaperPositionMode','auto');
     nR = length(PARAM3); nC = length(PARAM4);
    for p4_ii = 1 : length(PARAM4)
        cnt = cnt+1;
        tmpTxt = get_Parameters_titleText(PARAMETERS, [4], [ p4_ii ]);
        figure(Bfg)
        subplot(nR,nC, cnt); hold on;
        plot( PARAM5, squeeze(CntBurstSpkWT(p1_ii,p2_ii,p3_ii,p4_ii,:)),'*-k');
        plot( PARAM5, squeeze(CntBurstSpkKO(p1_ii,p2_ii,p3_ii,p4_ii,:)),'*-r');
        legend('WT','KO','location','Best')
        title(tmpTxt);
        xlabel(titleTxt5); ylabel('Avg Burst spike')
        
        figure(Bfg2)
        subplot(nR,nC, cnt); hold on;
        plot( PARAM5, squeeze(ChanceOfBurstingWT(p1_ii,p2_ii,p3_ii,p4_ii,:)),'*-k');
        plot( PARAM5, squeeze(ChanceOfBurstingKO(p1_ii,p2_ii,p3_ii,p4_ii,:)),'*-r');
         ylim([ 0 0.7])
        legend('WT','KO','location','Best')
        title(tmpTxt);
        xlabel(titleTxt5); ylabel('Chance of Burst spike')
        
    end
     tmpTxt = get_Parameters_titleText(PARAMETERS, [3], [ p3_ii ]);
    figure(Bfg)
    suptitle(tmpTxt);
    figure(Bfg2)
    suptitle([ 'Chance of Bursting, '  tmpTxt]);
    
    if (SAVE_FIG)
        tmpTxt = get_Parameters_saveText(PARAMETERS, [3], [ p3_ii ]);
        ffig = [ dirLoc dirFig 'NumBurstSpk' tmpTxt ];
        saveas(  Bfg, [ffig '.jpg'], 'jpg');        saveas(  Bfg, [ffig '.fig'], 'fig');
        ffig = [ dirLoc dirFig 'ChanceOfBurst' tmpTxt];
        saveas(  Bfg2, [ffig '.jpg'], 'jpg');        saveas(  Bfg2, [ffig '.fig'], 'fig');
        if (Close_Fig_aftr_save)
            close(Bfg); close(Bfg2);
        end
    end
end
