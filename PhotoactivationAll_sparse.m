%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%%% VL Basal activity : Center part  %%%%               %%%%%%%  Figure S4.A  %%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all
%%%%% NOTE :: This version did not include cell in the border 
%%% The range of connection in E cell is 200/sqrt(2)*3.5 =  494.9747 um
%%%                            I cell is 100/ sqrt(2)*3.5 = 247.4874 um
%%% Therefore, the cut border size is 500um
%%%% Get spiek train matrix from Membrane potential instead of raster plots
%%%% file that use automatic spike counts in NEURON
%% Parameter Settings
ssize = 1500; LowerBorder = 500; UpperBorder = ssize-LowerBorder;
THRESHOLD = 0;  % was -55 , in NEURON use 0
RES = 1; %temporal resolution of one bin
% CUTTIME = 500;
% PhotoInjT = 1000;
% PhotoStop = 1050;
% DelayT = 0;
% BurstRange = 100;
% getCenterPart = 0; 


% %% VL layer only [ give Poison Spike input to each cell in VL] 
%     dirLoc = 'VL_LocalConn/';
%     Name_postfix = 'TestSim_KO_InGauss0.2_W0.0005_50.00Hz_rEE0_rEI1_rIE2_rII0_Wmult1000';
% disp(Name_postfix)

%% Download Membrane potential of the network
[nE, nI, Tstop,somaVall] = get_save_vec_file_BigNet([dirLoc 'SomaVolt_' Name_postfix '.txt'] );
somaVallE = somaVall(1:nE,:);
somaVallI = somaVall(nE+1:nE+nI,:);
%% Download neurons location and get cell ID of the center neurons %note: It's not necessary (in this VL_Basal_ACt simulation) to find the Center ID everytime as the simulations use the same cells location

NNloc= importdata([dirLoc 'Neurons_location_' Name_postfix '.txt']); 
strE = 1; endE = nE;
strI = nE+1; endI= strI+nI-1;

Epos2D = NNloc(strE:endE,2:3); %2D -XY
E.xy = Epos2D; 
Ipos2D = NNloc(strI:endI,2:3); %2D -XY
I.xy = Ipos2D;
Epos3D = NNloc(strE:endE,2:4); %Not ceil , in 3D space
E.xyz= Epos3D;
Ipos3D = NNloc(strI:endI,2:4);
I.xyz = Ipos3D;    
if(getCenterPart)
%get Center Part
centerID_E = find((E.xy(:,1) > LowerBorder) & (E.xy(:,1) <UpperBorder)&(E.xy(:,2) > LowerBorder) & (E.xy(:,2) <UpperBorder));
centerID_I = find((I.xy(:,1) > LowerBorder) & (I.xy(:,1) <UpperBorder)&(I.xy(:,2) > LowerBorder) & (I.xy(:,2) <UpperBorder));
NcenterE = length(centerID_E); NcenterI = length(centerID_I); 
somaVall = somaVall ([centerID_E; nE+centerID_I],:);
somaVallE = somaVallE(centerID_E,:);
somaVallI = somaVallI(centerID_I,:);
N_E = nE; N_I = nI;
nE = size(somaVallE,1); nI = size(somaVallI,1);
end
%%
%  bb = importdata('TestSim/Input_SPKtrain_TestSim_KO_InGauss0.2_W0.0005_50.00Hz_rEE0_rEI2_rIE3_rII5_Wmult1000.txt');
% N_E = bb(1,1);
% N_I = bb(1,2);
% N_Trial = N_E + N_I;
% Tstop = bb(1,3);
% Trial_List = bb(2,:);
% bb = bb(3:end,:);
% allVec = zeros(N_Trial,size(bb,1));
% for i = 1:length(Trial_List)
%     allVec(i,:) = bb(:,i)';
% end
%% Making spiking train of cells 
PoisSpk = [];  E = []; I = [];
All.spktrain = zeros(nE+nI,Tstop+1); %t start at zero
E.spktrain = zeros(nE, Tstop+1); %t start at zero
I.spktrain = zeros(nI, Tstop+1); %t start at zero
fr_cells = zeros(nE+nI,1);
fr_cellsE = zeros(nE,1);
fr_cellsI = zeros(nI,1);
BurstSpk = zeros(nE+nI,1);
BurstSpkE = zeros(nE,1);
BurstSpkI = zeros(nI,1);
tmpT = PhotoStop + DelayT +1: PhotoStop + DelayT+BurstRange;
BurstSpkTrain = zeros(nE+nI,length(tmpT));
BurstSpkETrain = zeros(nE,length(tmpT));
BurstSpkITrain = zeros(nI,length(tmpT));
fr_AfterLightOff = zeros(nE+nI,1);
fr_AfterLightOffE = zeros(nE,1);
fr_AfterLightOffI = zeros(nI,1);

cuttime = CUTTIME;
for ii = 1:nE
    tmptrain = zeros(Tstop+1, 1);
    [pks,loc] = findpeaks(somaVallE(ii,:),'MINPEAKHEIGHT',THRESHOLD);
    tmptrain(round(loc*RES)) =1;
    E.spktrain(ii,:) = tmptrain;
    All.spktrain(ii,:) = tmptrain;
    fr_cellsE(ii) = sum(E.spktrain(ii,cuttime+1:Tstop))/ length(cuttime+1:Tstop)*1000; 
    fr_cells(ii) = fr_cellsE(ii);
    tmpT = PhotoStop + DelayT+1 : PhotoStop + DelayT+BurstRange;
    BurstSpkETrain(ii,:) = E.spktrain(ii,tmpT);  
    BurstSpkTrain(ii,:) = BurstSpkETrain(ii,:);
    BurstSpkE(ii) = sum(E.spktrain(ii,tmpT));  
    BurstSpk(ii) = BurstSpkE(ii);
    tmpT = PhotoStop + DelayT+BurstRange+ 1  : Tstop;
    fr_AfterLightOffE(ii) = sum(E.spktrain(ii,tmpT))/ length(tmpT)*1000;  % after photoactivation  + burstRange stop 
    fr_AfterLightOff(ii) = fr_AfterLightOffE(ii);
end

for ii = 1:nI
    tmptrain = zeros(Tstop+1, 1);
    [pks,loc] = findpeaks(somaVallI(ii,:),'MINPEAKHEIGHT',THRESHOLD);
    tmptrain(round(loc*RES)) =1;
    I.spktrain(ii,:) = tmptrain;
    All.spktrain(nE+ii,:) = tmptrain;
    fr_cellsI(ii) = sum(I.spktrain(ii,cuttime+1:Tstop))/ length(cuttime+1:Tstop)*1000; 
    fr_cells(nE+ii) = fr_cellsI(ii);
    tmpT = PhotoStop + DelayT+1 : PhotoStop + DelayT+BurstRange;
    BurstSpkITrain(ii,:) = I.spktrain(ii,tmpT);  
    BurstSpkTrain(nE+ii,:) = BurstSpkITrain(ii,:);
    BurstSpkI(ii) = sum(I.spktrain(ii,tmpT));
    BurstSpk(nE+ii) = BurstSpkI(ii);
    tmpT = PhotoStop + DelayT+ BurstRange + 1 : Tstop;
    fr_AfterLightOffI(ii) = sum(I.spktrain(ii,tmpT))/ length(tmpT)*1000;  
    fr_AfterLightOff(nE+ii) = fr_AfterLightOffI(ii);
end


E.spktrain = sparse(E.spktrain);
I.spktrain = sparse(I.spktrain);
All.spktrain = sparse(All.spktrain);

All.fr_data = fr_cells;
All.meanFR = mean(fr_cells);
All.var = var(fr_cells);
All.std = std(fr_cells);
All.sem = std(fr_cells)/sqrt(length(fr_cells));
All.BurstSpk = BurstSpk;
All.BurstSpkTrain = sparse(BurstSpkTrain);
All.fr_AfterLightOff = fr_AfterLightOff;
    
mE = mean(fr_cellsE);
varE = var(fr_cellsE);
stdE = std(fr_cellsE);
semE = stdE/sqrt(length(fr_cellsE));
E.fr_data = fr_cellsE; 
E.meanFR = mE;
E.var = varE;
E.std = stdE;
E.sem = semE;
E.BurstSpk = BurstSpkE;
E.BurstSpkTrain = sparse(BurstSpkETrain);
E.fr_AfterLightOff = fr_AfterLightOffE;


mI = mean(fr_cellsI);
varI = var(fr_cellsI);
stdI = std(fr_cellsI);
semI = stdI/sqrt(length(fr_cellsI));
I.fr_data = fr_cellsI; 
I.meanFR = mI;
I.var = varI;
I.std = stdI;
I.sem = semI;
I.BurstSpk = BurstSpkI;
I.BurstSpkTrain = sparse(BurstSpkITrain);
I.fr_AfterLightOff = fr_AfterLightOffI;


% %Test
% figure; 
% plot(somaVallE(ii,:)); hold on
% scatter(loc,pks);
% % xlim([0 1000]);

% % 
% InsertInputSPK = 0 ;
% if (InsertInputSPK)
%     aa = importdata([dirLoc 'raster_plots_paramOpt_' Name_postfix '.txt']);  %Just for now , Later move to get the spike train vectors from file
%     
%     nE = aa(1,1); nI=aa(1,2); Tstop =aa(1,3);
% tvec = round(aa(2:end,1)); 
% idvec = aa(2:end,2);
% strE = 1; endE = nE;
% strI = nE+1; endI= strI+nI-1;
% ssize = 1500;
% ncells = nE+nI;  
% totalCells = max(idvec(:))+1; 
% if (totalCells < nE + nI) 
%     totalCells = nE+nI;
% end
% nPoisSpk = totalCells - ncells; 
% strPoisSpk = endI+1; endPoisSpk = strPoisSpk + nPoisSpk -1;
% if(nPoisSpk > 0)   
%     PoisSpk.spktrain = zeros(nPoisSpk,Tstop+1); %t start at zero
%     spkMat = zeros(totalCells, Tstop+1);
%     for i = 1:length(tvec)
%         if(tvec(i) ~= 0)
%           spkMat(idvec(i)+1,tvec(i)+1) = 1;  %%%% spkmat( row = id, col = time(1-Tstop ms))
%         end
%     end   
%     %get average firing rate per cell
%     fr_PoisSpk = zeros(nPoisSpk,1);
% 
%     cuttime = CUTTIME; %500; %discard first 500 ms
%     PoisSpk.spktrain = zeros(nPoisSpk, Tstop+1); %t start at zero
%     for i = 1:nPoisSpk
%         PoisSpk.spktrain(i,:) = spkMat(ncells+i,1:end);
%         fr_PoisSpk(i) = sum(spkMat(ncells+i,cuttime+1:end)) / length(cuttime+1:Tstop)*1000; 
%     %     display(['Firing rate of input spike train#' num2str(i) ' is ' num2str(fr_PoisSpk(i)) ])
%     end
%     %Network
% mPoisSpk = mean(fr_PoisSpk);
% varPoisSpk = var(fr_PoisSpk);
% stdPoisSpk = std(fr_PoisSpk);
% semPoisSpk = std(fr_PoisSpk)/sqrt(length(fr_PoisSpk));
% PoisSpk.fr_data = fr_PoisSpk;
% PoisSpk.meanFR = mPoisSpk;
% PoisSpk.var = varPoisSpk;
% PoisSpk.std = stdPoisSpk;
% PoisSpk.sem = semPoisSpk;
% else
%     PoisSpk.fr_data = 0;
%     PoisSpk.meanFR = 0;
%     PoisSpk.var = 0;
%     PoisSpk.std = 0;
%     PoisSpk.sem = 0;
%     mPoisSpk =0;
%     stdPoisSpk=0;
%     semPoisSpk=0;
%     varPoisSpk=0;
%     
% end
% 
% end    
%     


disp('==================================================================================================')
disp([ 'Mean fr of E = ' num2str(mE)  ' , std = '  num2str(stdE) ' , SEM = ' num2str(semE)  ' , CV = ' num2str(stdE/mE) ] );
disp([ 'Mean fr of I = ' num2str(mI) ' , std = '  num2str(stdI) ' , SEM = ' num2str(semI)  ' , CV = ' num2str(stdI/mI) ] );
disp([ 'Mean fr of all = '  num2str(mean([fr_cellsE(:); fr_cellsI(:)]))]);
disp([ 'Burst Spike of all = ' num2str(mean(BurstSpk))]);
disp([ 'Fr after Light Off mean = ' num2str(mean(fr_AfterLightOff))]);

