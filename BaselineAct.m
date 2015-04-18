%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%
%%%% VL Basal activity : Center part  %%%%               %%%%%%%  Figure S4.A  %%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all
%%%%% NOTE :: This version did not include cell in the border if
%%%%% getCenterPart =1
%%% The range of connection in E cell is 200/sqrt(2)*3.5 =  494.9747 um
%%%                            I cell is 100/ sqrt(2)*3.5 = 247.4874 um
%%% Therefore, the cut border size is 500um
%%%% Get spike train matrix from Membrane potential instead of raster plots
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
if(getCenterPart)
    
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
%% Making spiking train of cells 
PoisSpk = [];  E = []; I = [];
All.spktrain = zeros(nE+nI,Tstop+1); %t start at zero
E.spktrain = zeros(nE, Tstop+1); %t start at zero
I.spktrain = zeros(nI, Tstop+1); %t start at zero
fr_cells = zeros(nE+nI,1);
fr_cellsE = zeros(nE,1);
fr_cellsI = zeros(nI,1);


cuttime = CUTTIME;
for ii = 1:nE
    tmptrain = zeros(Tstop+1, 1);
    [pks,loc] = findpeaks(somaVallE(ii,:),'MINPEAKHEIGHT',THRESHOLD);
    tmptrain(round(loc*RES)) =1;
    E.spktrain(ii,:) = tmptrain;
    All.spktrain(ii,:) = tmptrain;
    fr_cellsE(ii) = sum(E.spktrain(ii,cuttime+1:Tstop))/ length(cuttime+1:Tstop)*1000; 
    fr_cells(ii) = fr_cellsE(ii);
end

for ii = 1:nI
    tmptrain = zeros(Tstop+1, 1);
    [pks,loc] = findpeaks(somaVallI(ii,:),'MINPEAKHEIGHT',THRESHOLD);
    tmptrain(round(loc*RES)) =1;
    I.spktrain(ii,:) = tmptrain;
    All.spktrain(nE+ii,:) = tmptrain;
    fr_cellsI(ii) = sum(I.spktrain(ii,cuttime+1:Tstop))/ length(cuttime+1:Tstop)*1000; 
    fr_cells(nE+ii) = fr_cellsI(ii);
end

All.fr_data = fr_cells;
All.meanFR = mean(fr_cells);
All.var = var(fr_cells);
All.std = std(fr_cells);
All.sem = std(fr_cells)/sqrt(length(fr_cells));

    
mE = mean(fr_cellsE);
varE = var(fr_cellsE);
stdE = std(fr_cellsE);
semE = stdE/sqrt(length(fr_cellsE));
E.fr_data = fr_cellsE; 
E.meanFR = mE;
E.var = varE;
E.std = stdE;
E.sem = semE;


mI = mean(fr_cellsI);
varI = var(fr_cellsI);
stdI = std(fr_cellsI);
semI = stdI/sqrt(length(fr_cellsI));
I.fr_data = fr_cellsI; 
I.meanFR = mI;
I.var = varI;
I.std = stdI;
I.sem = semI;





disp('==================================================================================================')
disp([ 'Mean fr of E = ' num2str(mE)  ' , std = '  num2str(stdE) ' , SEM = ' num2str(semE)  ' , CV = ' num2str(stdE/mE) ] );
disp([ 'Mean fr of I = ' num2str(mI) ' , std = '  num2str(stdI) ' , SEM = ' num2str(semI)  ' , CV = ' num2str(stdI/mI) ] );
disp([ 'Mean fr of all = '  num2str(mean([fr_cellsE(:); fr_cellsI(:)]))]);


