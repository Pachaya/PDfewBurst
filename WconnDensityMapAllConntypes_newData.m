 
 
%% Load_Data

gauss = load('TC_Conn_INFO_gaussPgaussW_11-Jun-2015');
unif = load('TC_Conn_INFO_avgPuniformW_11-Jun-2015');
negExp = load('TC_Conn_INFO_avgPnegexpW_11-Jun-2015');
PARAMETERS = gauss.PARAMETERS;
AllConnData = []; 
AllConnData.gauss = gauss.TC_Conn_INFO;
AllConnData.unif = unif.TC_Conn_INFO;
AllConnData.negExp = negExp.TC_Conn_INFO;
tmp = AllConnData.gauss{1,1,1,1}; 
% length(tmp.TC_basedOnM1)

%%  Plot Sample Connection

% Load Cells Mosaics 

load('VL_M1_CellPosition') % Epos_VL  and Epos_M1
NNloc_VL = [zeros(length(Epos_VL),1) Epos_VL zeros(length(Epos_VL),1)];
NNloc_M1 = [zeros(length(Epos_M1),1) Epos_M1 50*ones(length(Epos_M1),1)];
% Plot NN loc 
figure; scatter(NNloc_M1(:,2), NNloc_M1(:,3)); axis square; 
rr = rand(166,2)*260;
figure; scatter(rr(:,1), rr(:,2)); axis square; xlim([-25 260+25]); ylim([-25 260+25]); box on;
rid = 2; wid = 3;
tmp = AllConnData.gauss{rid,wid,1,1};  Code = 'Convergent connection rule : Gaussian';
% tmp = AllConnData.unif{rid,wid,1,1};   Code = 'Convergent connection rule : Uniform';
% tmp = AllConnData.negExp{rid,wid,1,1};   Code = 'Convergent connection rule : Negative Exponential';

srclist = tmp.TCconn_raw.srcVLcell;
tarlist = tmp.TCconn_raw.tarM1Cell;
typelist = tmp.TCconn_raw.connType;
weightlist = tmp.TCconn_raw.connWeight;
delaylist = tmp.TCconn_raw.connDelay;

Build3DNetworkFigure_call

title(Code)
%% 
rid = 2; wid = 3; % Range = 100 , w = 200

ALLweightList = [ AllConnData.gauss{rid,wid,1,1}.TCconn_raw.connWeight;  AllConnData.unif{rid,wid,1,1}.TCconn_raw.connWeight;  AllConnData.negExp{rid,wid,1,1}.TCconn_raw.connWeight;  ];
maxC =prctile(ALLweightList,90); 
midC = prctile(ALLweightList,50); 
minC =prctile(ALLweightList,10);


AllsmoothW = [];
m3C =prctile(AllsmoothW(:),100); 
m2C = prctile(AllsmoothW(:),90); 
m1C =prctile(AllsmoothW(:),80);
%% Density plot for each connection types

rid = 2; wid = 3;
% tmp = AllConnData.gauss{rid,wid,1,1};  Code = 'Convergent connection rule : Gaussian';
% tmp = AllConnData.unif{rid,wid,1,1};   Code = 'Convergent connection rule : Uniform';
tmp = AllConnData.negExp{rid,wid,1,1};   Code = 'Convergent connection rule : Negative Exponential';

srclist = tmp.TCconn_raw.srcVLcell;
tarlist = tmp.TCconn_raw.tarM1Cell;
typelist = tmp.TCconn_raw.connType;
weightlist = tmp.TCconn_raw.connWeight;
delaylist = tmp.TCconn_raw.connDelay;
ALLweightList = [ ALLweightList; weightlist];


Range = PARAMETERS{1}.PARAM(rid);  
Rsize = ceil(Range / sqrt(2) * 4.5); 
[X,Y] = meshgrid(-Rsize : Rsize);
Wdense = zeros(size(X)); 


sampleM1ID = 87; 
M1_ID = 87;
nE_M1 = 166;
for M1_ID = 1 : nE_M1    
sampleNC = find( tarlist == tarID(sampleM1ID));
sumW = sum(weightlist(sampleNC));
for i = 1 :length(sampleNC)
    NC = sampleNC(i);
    srcID = srclist(NC);
    dff = Epos_VL(srcID,:) - Epos_M1(M1_ID,:) ;
%     dff = Epos_M1(M1_ID,:)  - Epos_VL(srcID,:) ;
    Wdense(mid + dff(1), mid + dff(2)) =  weightlist(NC) / sumW;    
end
end
figure; imagesc(-Rsize : Rsize, -Rsize : Rsize,Wdense); axis square xy;
caxis([minC maxC]); 

hcb=colorbar;
set(hcb,'YTick',[minC midC maxC])
set(hcb,'YTickLabel',{'10th', '50th', '90th'})

title(Code)

% Smooth image
fsig = 42;
smoothW = imgaussfilt(Wdense,fsig);
figure; imagesc(-Rsize : Rsize, -Rsize : Rsize, smoothW); axis square xy;
% colorbar(); %caxis([0.5 6]*10^(-5))

caxis([m1C m3C])
hcb=colorbar;
set(hcb,'YTick',[m1C m2C m3C])
set(hcb,'YTickLabel',{'80th', '90th', '100th'})

title(Code)
AllsmoothW = [ AllsmoothW; smoothW ];
%% 