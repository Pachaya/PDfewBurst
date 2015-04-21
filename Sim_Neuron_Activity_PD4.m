%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation activity  : 1 Pulse case 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time of event
% 0-500ms       : cut first 500ms out from initialized problem
% 500-1000ms    : normal activity
% 1000 - 1050ms : current injection period 
% 1050 -  burstRange 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setting 
ssize = 1500; LowerBorder = round(rTC/sqrt(2)*3 ); UpperBorder = ssize-LowerBorder;
 % Actually  round(rTC/sqrt(2)*3 )
THRESHOLD = -20;  % was -55 , in NEURON use 0
RES = 1; %temporal resolution of one bin
CUTTIME = 500;



%% Download spikes information  %[ spkTrainStruct, avgPoisSpk, avgE, avgI,avgAll ]=  get_avg_fr_from_raster_plots_cellsCN_PoisSpk(dirLoc, Name_postfix, CUTTIME, FLAG_disp)
aa = importdata([dirLoc 'raster_plots_paramOpt_' Name_postfix '.txt']);  
CUTTIME = avgFR_CUTTIME; %discard first 500 ms

M1 = [];  E = []; I = []; VL = []; 
nE = aa(1,1); nI=aa(1,2); tstop =aa(1,3);
spkMat = zeros(nE+nI, tstop);
if isempty(aa(2:end,2)) % no activity at all
    disp('Empty Rasterplot')
    E.spktrain = zeros(nE, tstop);
    E.fr_data =  zeros(nE, 1); 
    E.meanFR = 0;
    E.var = 0;
    E.std = 0;
    E.sem = 0;
     
    I.spktrain = zeros(nI, tstop);
    I.fr_data =  zeros(nI, 1);
    I.meanFR = 0;
    I.var = 0;
    I.std = 0;
    I.sem = 0;
    
    nE_M1  = 100;
    M1.normal.spkMatCut  =  zeros( nE_M1, tstop);
    M1.normal.fr_data = zeros( nE_M1,1);
    M1.inject.spkMatCut  =  zeros( nE_M1, tstop);
    M1.inject.fr_data = zeros( nE_M1,1);
    M1.burst.spkMatCut  =  zeros( nE_M1, tstop);
    M1.burst.fr_data = zeros( nE_M1,1);
    
    
    VL.normal.spkMatCut  =  zeros( nE, tstop);
    VL.normal.fr_data = zeros( nE,1);
    VL.inject.spkMatCut  =  zeros( nE, tstop);
    VL.inject.fr_data = zeros( nE,1);
    VL.burst.spkMatCut  =  zeros( nE, tstop);
    VL.burst.fr_data = zeros( nE,1);
    
    Condition = {['VL Fr : normal ' num2str(mean(VL.normal.fr_data)) ', injection ' num2str(mean(VL.inject.fr_data)) ', burst ' num2str(mean(VL.burst.fr_data))], ['M1 Fr : normal ' num2str(mean(M1.normal.fr_data)) ', injection ' num2str(mean(M1.inject.fr_data)) ', burst ' num2str(mean(M1.burst.fr_data))]};
    fg_rasterZoom = figure;
    set(fg_rasterZoom, 'position', [  947   108   598   873]);
    subplot(211);
    [freq_all_ALLnet,spkTime_all_Allnet] = raster_from_spkbin_noTitle( spkMat,InjStopT, Tstop, ['Full Raster plot' ]);
    rectangle('Position',[InjStartT,0,InjStopT-InjStartT,totalCells],'EdgeColor', 'k')
    rectangle('Position',[burstT1,0,burstT2-burstT1,totalCells],'EdgeColor', 'r')
    subplot(212);
    [freq_all_ALLnet,spkTime_all_Allnet] = raster_from_spkbin_noTitle( spkMat,InjStopT, Tstop, Condition);
    hold on;
    suptitle( [cTxt ', noise mean=' num2str(NoiseMEAN) ',noise sigma=' num2str(NoiseSIG) ', rTC' num2str(rTC) ', wmTC' num2str(wmTC) ])
    rectangle('Position',[InjStartT,0,InjStopT-InjStartT,totalCells],'EdgeColor', 'k')
    rectangle('Position',[burstT1,0,burstT2-burstT1,totalCells],'EdgeColor', 'r')
    xlim([InjStartT-50  burstT2+50]);
    ylim([ncells-50 ncells+50])
    if(SAVE_FIG)
        saveas(fg_rasterZoom, [dirLoc dirFig 'Raster' figNameCode '.fig'], 'fig')
        saveas(fg_rasterZoom, [dirLoc dirFig 'Raster' figNameCode '.jpg'], 'jpg')
    end        
    return
end

tvec = round(aa(2:end,1)); 
idvec = aa(2:end,2);
strE = 1; endE = nE;
strI = nE+1; endI= strI+nI-1;
ssize = 1500;
ncells = nE+nI;  
totalCells = max(idvec(:))+1; 
if (totalCells < nE + nI) 
    totalCells = nE+nI;
end

nM1  = totalCells - ncells;% no PoisSpk cell
strM1 = endI+1; endM1 = strM1 + nM1 -1;


%% Membrane potential of VL
[nE, nI, Tstop,somaVall] = get_save_vec_file_BigNet([dirLoc 'SomaVolt_' Name_postfix '.txt'] );
somaVallE = somaVall(1:nE,:);
somaVallI = somaVall(nE+1:nE+nI,:);


%% Membrane potential of M1 
[nE_M1, nI_M1, Tstop,somaVall_M1] = get_save_vec_file_BigNet([dirLoc 'SomaVolt_M1_' Name_postfix '.txt'] );
somaVallE_M1 = somaVall_M1(1:nE_M1,:);
somaVallI_M1 = somaVall_M1(nE_M1+1:nE_M1+nI_M1,:);
nM1 = nE_M1 + nI_M1; 
%% Download neurons location and get cell ID of the center neurons %note: It's not necessary (in this VL_Basal_ACt simulation) to find the Center ID everytime as the simulations use the same cells location
 
%getCenterPart_M1 = 1; 
NNloc= importdata([dirLoc 'Neurons_location_M1_' Name_postfix '.txt']); 
% strE = 1; endE = nE;
% strI = nE+1; endI= strI+nI-1;
strE_M1 = 1; endE_M1 = nE_M1;
strI_M1 = nE_M1 + 1; endI_M1 = strI_M1+nI_M1-1;

Epos2D = NNloc(strE_M1:endE_M1,2:3); %2D -XY
M1.E.xy = Epos2D; 
Ipos2D = NNloc(strI_M1:endI_M1,2:3); %2D -XY
M1.I.xy = Ipos2D;
Epos3D = NNloc(strE_M1:endE_M1,2:4); %Not ceil , in 3D space
M1.E.xyz= Epos3D;
Ipos3D = NNloc(strI_M1:endI_M1,2:4);
M1.I.xyz = Ipos3D;    
if(getCenterPart_M1)
%get Center Part of M1
centerID_E = find((M1.E.xy(:,1) > LowerBorder) & (M1.E.xy(:,1) <UpperBorder)&(M1.E.xy(:,2) > LowerBorder) & (M1.E.xy(:,2) <UpperBorder));
centerID_I = find((M1.I.xy(:,1) > LowerBorder) & (M1.I.xy(:,1) <UpperBorder)&(M1.I.xy(:,2) > LowerBorder) & (M1.I.xy(:,2) <UpperBorder));
NcenterE = length(centerID_E); NcenterI = length(centerID_I); 
somaVall_M1_backup = somaVall_M1;
somaVall_M1 = somaVall_M1([centerID_E; nE_M1+centerID_I],:);
somaVallE_M1 = somaVallE_M1(centerID_E,:);
somaVallI_M1 = somaVallI_M1(centerID_I,:);
N_E_M1 = nE_M1; N_I_M1 = nI_M1;
nE_M1 = size(somaVallE_M1,1); nI_M1 = size(somaVallI_M1,1);
nM1 = nE_M1 + nI_M1; 
end

%% Making spiking train of cells 
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
%     fr_cellsE(ii) = sum(E.spktrain(ii,:))/ length(cuttime+1:Tstop)*1000; 
%     fr_cells(ii) = fr_cellsE(ii);
end
for ii = 1:nI
    tmptrain = zeros(Tstop+1, 1);
    [pks,loc] = findpeaks(somaVallI(ii,:),'MINPEAKHEIGHT',THRESHOLD);
    tmptrain(round(loc*RES)) =1;
    I.spktrain(ii,:) = tmptrain;
    All.spktrain(nE+ii,:) = tmptrain;
%     fr_cellsI(ii) = sum(I.spktrain(ii,:))/ length(cuttime+1:Tstop)*1000; 
%     fr_cells(ii) = fr_cellsI(ii);
end
% 
% All.fr_data = fr_cells;
% All.meanFR = mean(fr_cells);
% All.var = var(fr_cells);
% All.std = std(fr_cells);
% All.sem = std(fr_cells)/sqrt(length(fr_cells));
% 
% mE = mean(fr_cellsE);
% varE = var(fr_cellsE);
% stdE = std(fr_cellsE);
% semE = stdE/sqrt(length(fr_cellsE));
% E.fr_data = fr_cellsE; 
% E.meanFR = mE;
% E.var = varE;
% E.std = stdE;
% E.sem = semE;
% 
% mI = mean(fr_cellsI);
% varI = var(fr_cellsI);
% stdI = std(fr_cellsI);
% semI = stdI/sqrt(length(fr_cellsI));
% I.fr_data = fr_cellsI; 
% I.meanFR = mI;
% I.var = varI;
% I.std = stdI;
% I.sem = semI;


%% M1 
M1.All.spktrain = zeros(nE_M1+nI_M1,Tstop+1); %t start at zero
M1.E.spktrain = zeros(nE_M1, Tstop+1); %t start at zero
M1.I.spktrain = zeros(nI_M1, Tstop+1); %t start at zero
fr_cells_M1 = zeros(nE_M1+nI_M1,1);
fr_cellsE_M1 = zeros(nE_M1,1);
fr_cellsI_M1 = zeros(nI_M1,1);
cuttime = CUTTIME;
for ii = 1:nE_M1
    tmptrain = zeros(Tstop+1, 1);
    [pks,loc] = findpeaks(somaVallE_M1(ii,:),'MINPEAKHEIGHT',THRESHOLD);
    tmptrain(round(loc*RES)) =1;
    M1.E.spktrain(ii,:) = tmptrain;
    M1.All.spktrain(ii,:) = tmptrain;
%     fr_cellsE_M1(ii) = sum(E.spktrain(ii,:))/ length(cuttime+1:Tstop)*1000; 
%     fr_cells_M1(ii) = fr_cellsE_M1(ii);
end
for ii = 1:nI
    tmptrain = zeros(Tstop+1, 1);
    [pks,loc] = findpeaks(somaVallI_M1(ii,:),'MINPEAKHEIGHT',THRESHOLD);
    tmptrain(round(loc*RES)) =1;
    M1.I.spktrain(ii,:) = tmptrain;
    M1.All.spktrain(nE+ii,:) = tmptrain;
%     fr_cellsI_M1(ii) = sum(M1.I.spktrain(ii,:))/ length(cuttime+1:Tstop)*1000; 
%     fr_cells_M1(ii) = fr_cellsI_M1(ii);
end

% M1.All.fr_data = fr_cells_M1;
% M1.All.meanFR = mean(fr_cells_M1);
% M1.All.var = var(fr_cells_M1);
% M1.All.std = std(fr_cells_M1);
% M1.All.sem = std(fr_cells_M1)/sqrt(length(fr_cells_M1));
% 
% mE_M1 = mean(fr_cellsE_M1);
% varE_M1 = var(fr_cellsE_M1);
% stdE_M1 = std(fr_cellsE_M1);
% semE_M1 = stdE/sqrt(length(fr_cellsE_M1));
% M1.E.fr_data = fr_cellsE_M1; 
% M1.E.meanFR = mE_M1;
% M1.E.var = varE_M1;
% M1.E.std = stdE_M1;
% M1.E.sem = semE_M1;
% 
% mI_M1 = mean(fr_cellsI_M1);
% varI_M1 = var(fr_cellsI_M1);
% stdI_M1 = std(fr_cellsI_M1);
% semI_M1 = stdI/sqrt(length(fr_cellsI_M1));
% M1.I.fr_data = fr_cellsI_M1; 
% M1.I.meanFR = mI_M1;
% M1.I.var = varI_M1;
% M1.I.std = stdI_M1;
% M1.I.sem = semI_M1;

% %Test
% figure; plot(somaVallE(1,:)); hold on 
% scatter(loc,pks);
% xlim([0 1000]);



%%
%% Average firing rate 
normalT1 = cutTime; normalT2 = InjStartT; 
injectT1 = InjStartT; injectT2 = InjStopT;
% reboundDelay = 70; %Input_I_dur;
burstT1 = InjStopT+reboundDelay;  burstT2 = burstT1+burstRange;
disp('==================================================================================================')
VL=[];
VL.normal = []; VL.inject =[]; VL.burst =[];
GET_ONLY_E = 1
if( GET_ONLY_E)
    spkMat = [E.spktrain; M1.E.spktrain;]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GetOnlyE
     ncells = nE; nM1 = nE_M1;
else
    spkMat = [All.spktrain; M1.All.spktrain;]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% combine VL and M1
end
disp(['VL: average Firing rate during normal activity : T = (' num2str(normalT1) ' , ' num2str(normalT2) ' )']) ;
VL.normal.spkMatCut = spkMat(1:ncells,normalT1 : normalT2);
VL.normal.fr_data = sum(VL.normal.spkMatCut,2)./ length(normalT1:normalT2)*1000;
disp(num2str(mean(VL.normal.fr_data)))

disp(['VL: average Firing rate during inhibition : T = (' num2str(injectT1) ' , ' num2str(injectT2) ' )']) ;
VL.inject.spkMatCut = spkMat(1:ncells, injectT1 : injectT2);
VL.inject.fr_data = sum(VL.inject.spkMatCut,2)./ length(injectT1:injectT2)*1000;
disp(num2str(mean(VL.inject.fr_data)))

disp(['VL: average Firing rate during bursting activity : T = (' num2str(burstT1) ' , ' num2str(burstT2) ' )']) ;
VL.burst.spkMatCut = spkMat(1:ncells, burstT1 : burstT2);
VL.burst.fr_data = sum(VL.burst.spkMatCut,2)./ length(burstT1:burstT2)*1000;
disp(num2str(mean(VL.burst.fr_data)))

disp('==================================================================================================')
M1.normal = []; M1.inject =[]; M1.burst =[];
disp(['M1: average Firing rate during normal activity : T = (' num2str(normalT1) ' , ' num2str(normalT2) ' )']) ;
M1.normal.spkMatCut = spkMat(ncells+1:ncells+nM1,normalT1 : normalT2);
M1.normal.fr_data = sum(M1.normal.spkMatCut,2)./ length(normalT1:normalT2)*1000;
disp(num2str(mean(M1.normal.fr_data)))

disp(['M1: average Firing rate during inhibition : T = (' num2str(injectT1) ' , ' num2str(injectT2) ' )']) ;
M1.inject.spkMatCut = spkMat(ncells+1:ncells+nM1, injectT1 : injectT2);
M1.inject.fr_data = sum(M1.inject.spkMatCut,2)./ length(injectT1:injectT2)*1000;
disp(num2str(mean(M1.inject.fr_data)))

disp(['M1: average Firing rate during bursting activity : T = (' num2str(burstT1) ' , ' num2str(burstT2) ' )']) ;
M1.burst.spkMatCut = spkMat(ncells+1:ncells+nM1, burstT1 : burstT2);
M1.burst.fr_data = sum(M1.burst.spkMatCut,2)./ length(burstT1:burstT2).*1000;
disp(num2str(mean(M1.burst.fr_data)))
disp('==================================================================================================')
%%
% [freq_all_ALLnet,spkTime_all_Allnet, fg_handle_Allnet] = raster_from_spkbin( spkMat,avgFR_CUTTIME, Tstop, [cTxt ', noise mean=' num2str(NoiseMEAN) ', rTC' num2str(rTC) ', wmTC' num2str(wmTC) ]);
% set(fg_handle_Allnet ,'position',[   1165         414         595         420]);
% hold on;
%%
Condition = {['VL Fr : normal ' num2str(mean(VL.normal.fr_data)) ', injection ' num2str(mean(VL.inject.fr_data)) ', burst ' num2str(mean(VL.burst.fr_data))], ['M1 Fr : normal ' num2str(mean(M1.normal.fr_data)) ', injection ' num2str(mean(M1.inject.fr_data)) ', burst ' num2str(mean(M1.burst.fr_data))]};
fg_rasterZoom = figure;
set(fg_rasterZoom, 'position', [  947   108   598   873]);
subplot(211);
[freq_all_ALLnet,spkTime_all_Allnet] = raster_from_spkbin_noTitle( spkMat,avgFR_CUTTIME, Tstop, ['Full Raster plot' ]);
rectangle('Position',[InjStartT,0,InjStopT-InjStartT,totalCells],'EdgeColor', 'k')
rectangle('Position',[burstT1,0,burstT2-burstT1,totalCells],'EdgeColor', 'r')
subplot(212);
[freq_all_ALLnet,spkTime_all_Allnet] = raster_from_spkbin_noTitle( spkMat,avgFR_CUTTIME, Tstop, Condition);
hold on;
suptitle( [cTxt ', noise mean=' num2str(NoiseMEAN) ',noise sigma=' num2str(NoiseSIG) ', rTC' num2str(rTC) ', wmTC' num2str(wmTC) ])
rectangle('Position',[InjStartT,0,InjStopT-InjStartT,totalCells],'EdgeColor', 'k')
rectangle('Position',[burstT1,0,burstT2-burstT1,totalCells],'EdgeColor', 'r')
xlim([InjStartT-100  burstT2+100]);
ylim([ncells-50 ncells+50])

if(SAVE_FIG)
        saveas(fg_rasterZoom, [dirLoc dirFig 'Raster' figNameCode '.fig'], 'fig')
        saveas(fg_rasterZoom, [dirLoc dirFig 'Raster' figNameCode '.jpg'], 'jpg')
end        

% line([InjStartT InjStartT], [0 totalCells],'Color', 'k');
% line([InjStopT InjStopT], [0 totalCells],'Color', 'k');
% line([burstT1 burstT1], [0 totalCells],'Color', 'r');
% line([burstT2 burstT2], [0 totalCells],'Color', 'r');

% xlim([1450 1700]);
% ylim([1100 1200]);

%   saveas(fg_rasterZoom,[dirLoc 'Fig_Rasterplot_' cTxt '_' Name_postfix '.fig'],'fig');
%   saveas(fg_rasterZoom,[dirLoc 'Fig_Rasterplot_' cTxt '_' Name_postfix '.jpg'],'jpg');

% %% from raster
% if(nM1 > 0)   
%     M1.spktrain = zeros(nM1,Tstop+1); %t start at zero
%     spkMat = zeros(totalCells, Tstop+1);
%     for i = 1:length(tvec)
%         if(tvec(i) ~= 0)
%           spkMat(idvec(i)+1,tvec(i)+1) = 1;  %%%% spkmat( row = id, col = time(1-Tstop ms))
%         end
%     end   
% %     get average firing rate per cell
%     fr_M1 = zeros(nM1,1);
% 
%     cuttime = avgFR_CUTTIME; %500; %discard first 500 ms
%     M1.spktrain = zeros(nM1, Tstop+1); %t start at zero
%     for i = 1:nM1
%         M1.spktrain(i,:) = spkMat(ncells+i,1:end);
%         fr_M1(i) = sum(spkMat(ncells+i,cuttime+1:end)) / length(cuttime+1:Tstop)*1000; 
% %         display(['Firing rate of input spike train#' num2str(i) ' is ' num2str(fr_M1(i)) ])
%     end
% %     Network
% mM1 = mean(fr_M1);
% varM1 = var(fr_M1);
% stdM1 = std(fr_M1);
% semM1 = std(fr_M1)/sqrt(length(fr_M1));
% M1.fr_data = fr_M1;
% M1.meanFR = mM1;
% M1.var = varM1;
% M1.std = stdM1;
% M1.sem = semM1;
% else
%     M1.fr_data = 0;
%     M1.meanFR = 0;
%     M1.var = 0;
%     M1.std = 0;
%     M1.sem = 0;
%     mM1 =0;
%     stdM1=0;
%     semM1=0;
%     varM1=0;
%     
% end















