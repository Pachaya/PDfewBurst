
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

DoTTest = 1;
SAVE_BASAL_ACT = 0;
FNAME_SIM = 'GPmVL' ; %'PreLim_INOISE0.2_InSPK50_EE1_EI2_IE20_II12_mult125_w0.0006';
plotSampleV = 0;
SPECIFIED_RATIO = 0; %TestSim_WT_InGauss0.2_W0.0005_50.00Hz_rEE0_rEI1_rIE2_rII1_Wmult10
Basal_Act_set  = [];

NoiseMEAN = 0;

SAVE_FIG = 1;

dirLoc = 'Model1Min150403_TCconvergence/' ;
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
dirFig = 'Fig_testrTC/';
mkdir([dirLoc dirFig])
TSTOP = 4000;

rTC_LST = [50 75 100];
wmTC_LST = [20 30 40 50  ];


LightAmp_LST = [0.4; ];
GPmVLw_mean_LST = [0.08];
GPmVLw_sig_LST =[ 0.01 ];

% Candidates Model 1 :
% I amp = 0.3, mean 0.06, sigma = 0.01, 0.005 ?0.03
% I amp = 0.5, mean 0.04, sigma = 0.03, 0.01 ?0.03
%PDfewBurst_GPmVLmd1_rTC120_wmTC10_WT_GPmInput_Amp0.3_Dur1000_GPmVLw_m0.06_sig0.01_InGauss0.2_IGmean-0.15_IGmeanSig0_W0.0015_SpecifiedPoisSpk_sig0.00Hz_T4000_trial3
NUM_TRIAL = 5;
ncells = 1150;
TRIAL_LST = 1 : NUM_TRIAL;

CUTTIME = 500;
PhotoInjT = 1500;
PhotoStop = 2500;
DelayT = 0;
BurstRange = 100;


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
                        
                        
                        coreFileName = 'PDfewBurst_GPmVLmd1_0del' ;
                        
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
                        
                        
                
                        
                        %                 CheckFileExist( dirLoc, Name_postfix  )
                        disp('######  Download VL ')
                        Name_postfix = [ Simulation_Code];
                        %                         VL_PhotoactivationAll
                        
                        %                 if(0)
                        
                         tmp.E = E;      tmp.I = I;  tmp.All = All; %tmp.VL = VL; tmp.M1 = M1;
                        if (cell_type == 1)
                             VL.WT = tmp;
                            cTxt = 'WT'; Name_postfix_WT = Name_postfix;
                        elseif (cell_type == 2)
                            VL.KO = tmp;
                            cTxt = 'KO'; Name_postfix_KO = Name_postfix;
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
                        end
                        
                        disp('Finish plot sample membrain voltage');
                        toc(tt)
                    end
                    
                    
                    
                     Basal_Act.VL = VL;
                   Basal_Act.M1 = M1;
                    
                    
                end
                
                                   ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii } = Basal_Act;
                
            end
        end
    end
end
%% 
if(SAVE_BASAL_ACT)
    save([dirLoc dirFig 'Saved_Activity_result_' date '.mat' ], 'ACT_Record','PARAMETERS','-v7.3');
end

%% 
%Make only the neuron activity during burst

TimeofFirstReboundSpk = zeros(ACT_Rec_size);
p1_ii = 1; p2_ii=1;
for p3_ii = 1 : length(PARAM3)
    fRB = figure; set(fRB,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
    fBspk = figure; set(fBspk,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
    cumf = figure; set(fBspk,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
    cumfspk = figure; set(fBspk,'position',[302          49        1290         948]);  set(gcf,'PaperPositionMode','auto')
    cnt = 0;
    nR = length(PARAM4); nC = length(PARAM5);
    for p4_ii = 1 : length(PARAM4)
        for p5_ii = 1 : length(PARAM5)
            la_ii = p3_ii; m_ii = p4_ii; s_ii = p5_ii;
            
          
            BasalAct = ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii}.VL;
            hbin =1;
            histBin = 1:hbin:size(BasalAct.WT.All.BurstSpkTrain,2); 
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
            [NWT,X] = hist(tmpRebound_WT,hbin);
            pH= bar(X,NWT/ncells*100,'FaceColor','k','EdgeColor','w'); hold on; 
            cH = get(pH,'Children');
            set(cH,'FaceAlpha',0.7); % 0 = transparent, 1 = opaque.

            [NKO,X] = hist(tmpRebound_KO,hbin);
            pH=bar(X,NKO/ncells*100,'FaceColor','r','EdgeColor','w'); hold on; 
            cH = get(pH,'Children');
            set(cH,'FaceAlpha',0.7); % 0 = transparent, 1 = opaque.

            ylabel('% Neurons'); xlabel('Time(ms)')
            title([ titleTxt4 ' = ' num2str(PARAM4(p4_ii)) ', ' titleTxt5 ' = ' num2str(PARAM5(p5_ii))])
             legend('WT','KO')

             
            figure(cumf); 
            subplot(nR,nC,cnt);
            plot(cumsum(NWT),'k'); hold on;
             plot(cumsum(NKO),'r'); hold on;
          legend('WT','KO')


            % Distribution of spikes in bursting range 
           
             tmpRebound_WT = sum(BasalAct.WT.All.BurstSpkTrain);
             tmpRebound_KO = sum(BasalAct.KO.All.BurstSpkTrain);
             xx = 1:1:size(BasalAct.WT.All.BurstSpkTrain,2);
            figure(fBspk); 
            subplot(nR,nC,cnt);
            pH=bar(xx,tmpRebound_WT/ncells*100,'k'); hold on; 
            cH = get(pH,'Children');
            set(cH,'FaceAlpha',0.7); % 0 = transparent, 1 = opaque.

            pH=bar(xx,tmpRebound_KO/ncells*100,'FaceColor','r','EdgeColor','r');
                        cH = get(pH,'Children');
            set(cH,'FaceAlpha',0.7); % 0 = transparent, 1 = opaque.
            ylabel('% Neurons'); xlabel('Time(ms)')
            title([ titleTxt4 ' = ' num2str(PARAM4(p4_ii)) ', ' titleTxt5 ' = ' num2str(PARAM5(p5_ii))])
             legend('WT','KO')
            figure(cumfspk); 
            subplot(nR,nC,cnt);
            plot(cumsum(NWT),'k'); hold on;
            plot(cumsum(NKO),'r'); hold on;
            legend('WT','KO');
            clear BasalAct spkBin
        end
    end
    figure(fRB); suptitle(['Latency of first Rebound spike : ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]);  % Note equal to the latency of peak firing rate --> latency of peak firing rate calculate from instantaneous firing rate
    figure(cumf); suptitle(['Latency of first Rebound spike : ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]); 
    figure(fBspk); suptitle(['Percent of neuron that fire during bursting period : ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]);  % Note equal to the latency of peak firing rate --> latency of peak firing rate calculate from instantaneous firing rate
    figure(cumfspk); suptitle(['Percent of neuron that fire during bursting period : ' titleTxt3 '=' num2str(PARAM3(p3_ii)) ]);  %

end

if (SAVE_FIG)
        tmpTxt = get_Parameters_saveText(PARAMETERS, [3], [ p3_ii ]);
        ffig = [ dirLoc dirFig 'FirstReboundSpike'  tmpTxt];
        saveas(fRB, [ffig '.fig'] , 'fig');
        saveas(fRB, [ffig '.jpg'] , 'jpg');
        ffig = [ dirLoc dirFig 'AllBurstSpike'  tmpTxt]; % Better minus by expected baseline activity
        saveas(fBspk, [ffig '.fig'] , 'fig');
        saveas(fBspk, [ffig '.jpg'] , 'jpg');
        ffig = [ dirLoc dirFig 'FirstReboundSpikeCumSum'  tmpTxt];
        saveas(cumf, [ffig '.fig'] , 'fig');
        saveas(cumf, [ffig '.jpg'] , 'jpg');
        ffig = [ dirLoc dirFig 'AllBurstSpikeCumSum'  tmpTxt]; % Better minus by expected baseline activity
        saveas(cumfspk, [ffig '.fig'] , 'fig');
        saveas(cumfspk, [ffig '.jpg'] , 'jpg');
        
        
        if (Close_Fig_aftr_save)
            close(fRB); close(fBspk);
        end
        
end
