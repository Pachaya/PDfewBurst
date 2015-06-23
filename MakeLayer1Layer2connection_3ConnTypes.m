clc
close all
clear all
% Make connection for layer 1 - layer 2 with three types of conenction
% rules
% MakeLayer1Layer2connection_3ConnTypes.m

PROP_SPEED = 1000 ; % [um/ms] speed of axonal propagation < 1 m/s >
layer1locTxt = 'Epos_E1150_I380_Regular_VL_Sweep5Itr5000_ceil_sorted';
layer2locTxt = 'ECenterpos_rTC250_E5625_I1892_ssize1500_Regular_sweep5Itr5000_ceil_sorted';

SimCode = 'Theory';
[L1pos, L1Distmat, L1nnDist]= GetCellsDistInfo(layer1locTxt, 1);
[L2pos, L2Distmat, L2nnDist]= GetCellsDistInfo(layer2locTxt, 1); %VL
DistMat = pdist2(L1pos,L2pos);
PLOT_FIG = 0;
%%
close all
W_SCALE = 0.00001; 
W_LST = [50:25:100]; % 50 : 25 : 100
% Range_LST = [10 25:25:300];
Range_LST = [50:50:200];

saveConn_g = cell(  length(Range_LST),  length(W_LST));     saveConn_u = cell(  length(Range_LST),  length(W_LST));     saveConn_e = cell(  length(Range_LST),  length(W_LST));
rw_NumCon_g = zeros(  length(Range_LST),  length(W_LST));     rw_NumCon_u = zeros(  length(Range_LST),  length(W_LST));     rw_NumCon_e = zeros(  length(Range_LST),  length(W_LST)); 
rw_sumW_g = zeros(  length(Range_LST),  length(W_LST));     rw_sumW_u = zeros(  length(Range_LST),  length(W_LST));     rw_sumW_e = zeros(  length(Range_LST),  length(W_LST));

TRIAL = 5;
rng(TRIAL)
for rr = 1 : length(Range_LST)
    for ww = 1 : length(W_LST)
tt = tic();
range = Range_LST(rr);
weight_factor = W_LST(ww);



%%
gConnMat = DistMat.*0;  uConnMat = DistMat.*0;  eConnMat = DistMat.*0;
gW_Mat = DistMat.*0;    uW_Mat = DistMat.*0;    eW_Mat = DistMat.*0;
gDelay_Mat = DistMat.*0;    uDelay_Mat = DistMat.*0;    eDelay_Mat = DistMat.*0;
Pmax = 0.85;

SumConnectivity_g =  zeros( length(L2pos), 1);
NumConnection_g = zeros( length(L2pos), 1);

SumConnectivity_s =  zeros( length(L2pos), 1);
NumConnection_s = zeros( length(L2pos), 1);

SumTotalConnStrength_g =  zeros( length(L2pos), 1);
SumTotalConnStrength_u =  zeros( length(L2pos), 1);
SumTotalConnStrength_e =  zeros( length(L2pos), 1);

SumConnStrength_g =  zeros( length(L2pos), 1);
SumConnStrength_u =  zeros( length(L2pos), 1);
SumConnStrength_e =  zeros( length(L2pos), 1);

% Note : For Volume comparison , take summation values of all the cells
% within range before make connection. This way, it's the closest to the
% real volumn

for l2_ID = 1 : length(L2pos) % for each cell in layer 2
    %%% only the cells in range
    potentialLayer1 = find(DistMat(:,l2_ID) <= range*3 ); % Cells that located within range
    
    %%% Convergent Connection Rule 1
%     disp('--------------------- G-G ---------------------');
    % ############## Connectivity : Gaussian Distribution
    pp = func_gauss3(Pmax, DistMat(potentialLayer1,l2_ID),range);
    xx = rand(size(pp));
    tmpVec_g =  xx < pp; sumProb = sum(pp);  
    tmpConnID = potentialLayer1(tmpVec_g);
    nCon_of_g = length(tmpConnID); 
       
    % Volume Comparison
    SumConnectivity_g(l2_ID) =  sumProb;
    NumConnection_g(l2_ID) = nCon_of_g;
    
    % Recorded Connection
    gConnMat(tmpConnID,l2_ID) = 1; % 1 = connection exist, 0 = no connection;    
    gDelay_Mat(tmpConnID,l2_ID) = DistMat(tmpConnID, l2_ID)./PROP_SPEED;
    
    % ############## Connection Strength : Gaussian Distribution
    tmpW = func_gauss3(weight_factor*W_SCALE, DistMat(potentialLayer1, l2_ID),range); %for all possible connection within range 
    WforConnectedCells =  tmpW(tmpVec_g);    
    
    %Volume Comparison
    SumTotalConnStrength_g(l2_ID) = sum(tmpW);
    SumConnStrength_g(l2_ID) =  sum(WforConnectedCells);
    
    % Recorded Connection Strength
    gW_Mat(tmpConnID,l2_ID) = WforConnectedCells;
    
    
%     disp('--------------------- Sparse Connectivity in S-U and S-E ---------------------');
    %%% Convergent Connection Rule 2 and 3 use the same kind of
    %%% connectivity = sparse connection from uniform connectivity  ----> Make Connectivity first
   
    % ############## Connectivity : Uniform Distribution
    pGauss = sumProb; % From Gaussian
    pp_s = ones(size(potentialLayer1)).*(pGauss / length(potentialLayer1)); % Average weight per one L1 cell
    xx = rand(size(pp_s));  sumProb_s = sum(pp_s); 
    tmpVec_s =  xx < pp_s;  
    tmpConnID_s = potentialLayer1(tmpVec_s);   
    nCon_of_s = length(tmpConnID_s); 
       
    % Volume Comparison
    SumConnectivity_s(l2_ID) =  sumProb_s;
    NumConnection_s(l2_ID) = nCon_of_s;
    
    % Recorded Connection
    % S-U
    uConnMat(tmpConnID_s,l2_ID) = 1; % 1 = connection exist, 0 = no connection;    
    uDelay_Mat(tmpConnID_s,l2_ID) = DistMat(tmpConnID_s, l2_ID)./PROP_SPEED;
    % S-E
    eConnMat(tmpConnID_s,l2_ID) = 1; % 1 = connection exist, 0 = no connection;    
    eDelay_Mat(tmpConnID_s,l2_ID) = DistMat(tmpConnID_s, l2_ID)./PROP_SPEED;
    
    
    %%% Convergent Connection Rule 2 : Connectivity - Uniform , Strength - Uniform
%     disp('-----------------------S-U ---------------------');

    % ############## Connection Strength : Uniform Distribution
    wGauss = SumConnStrength_g(l2_ID); % sum(func_gauss3(weight_factor*W_SCALE, DistMat(tmpConnID_s,l2_ID),sig)); % The possible weight in Gaussian 
    avgW = wGauss ./ length(tmpConnID_s); 
    tmpW = ones(size(tmpConnID_s)).*avgW;
    WforConnectedCells_u = tmpW; %  % For Uniform distribution
    
    %Volume Comparison
    SumTotalConnStrength_u(l2_ID) = sum(tmpW);
    SumConnStrength_u(l2_ID) =  sum(WforConnectedCells_u);
    
    % Recorded Connection Strength
    uW_Mat(tmpConnID_s,l2_ID) = WforConnectedCells_u;
    
    %%% Convergent Connection Rule 3 : Connectivity - Uniform , Strength -
    %%% Negative Exponential
%     disp('-----------------------S-E ---------------------');

    % ############## Connection Strength : Exponential Distribution
    wGauss = SumConnStrength_g(l2_ID); % sum(func_gauss3(weight_factor*W_SCALE, DistMat(tmpConnID_s,l2_ID),sig)); % The possible weight in Gaussian 
    avgW = wGauss./ length(tmpConnID_s); 
    tmpW = exprnd(avgW, size(tmpConnID_s));
   
    % control sum W to be same 
    tmpW = tmpW.* wGauss./ sum(tmpW); % already checked that distribution of tmpW is same before and after normalized
    WforConnectedCells_e = tmpW;%( tmpVec_s);  % For Exponential distribution
    
    %Volume Comparison
    SumTotalConnStrength_e(l2_ID) = sum(tmpW);
    SumConnStrength_e(l2_ID) =  sum(WforConnectedCells_e);
    
    % Recorded Connection Strength
    eW_Mat(tmpConnID_s,l2_ID) = WforConnectedCells_e;  
    
%     %%% Convergent Connection Rule 2 : Connectivity - Uniform , Strength - Uniform
%     disp('-----------------------S-U ---------------------');
% 
%     % ############## Connection Strength : Uniform Distribution
%     wGauss = sum(func_gauss3(weight_factor*W_SCALE, DistMat(potentialLayer1,l2_ID),sig)); % The possible weight in Gaussian 
%     avgW = wGauss ./ length(potentialLayer1); 
%     tmpW = ones(size(potentialLayer1)).*avgW;
%     WforConnectedCells_u = tmpW( tmpVec_s);  % For Uniform distribution
%     
%     %Volume Comparison
%     SumTotalConnStrength_u(l2_ID) = sum(tmpW);
%     SumConnStrength_u(l2_ID) =  sum(WforConnectedCells_u);
%     
%     % Recorded Connection Strength
%     gW_Mat(tmpConnID_s,l2_ID) = WforConnectedCells_u;
%     
%     
%     %%% Convergent Connection Rule 3 : Connectivity - Uniform , Strength -
%     %%% Negative Exponential
%     disp('-----------------------S-E ---------------------');
% 
%     % ############## Connection Strength : Exponential Distribution
%     wGauss = sum(func_gauss3(weight_factor*W_SCALE, DistMat(potentialLayer1,l2_ID),sig)); % The possible weight in Gaussian 
%     avgW = wGauss./ length(potentialLayer1); 
%     tmpW = exprnd(avgW, size(potentialLayer1));
%     WforConnectedCells_e = tmpW( tmpVec_s);  % For Exponential distribution
%     
%     %Volume Comparison
%     SumTotalConnStrength_e(l2_ID) = sum(tmpW);
%     SumConnStrength_e(l2_ID) =  sum(WforConnectedCells_e);
%     
%     % Recorded Connection Strength
%     gW_Mat(tmpConnID_s,l2_ID) = WforConnectedCells_e;  
if(l2_ID == 1)&&(PLOT_FIG)
    figure; set(gcf,'position',[  350         380        1369         420]);
    subplot(131); hist(WforConnectedCells); title(['Gaussian, sum W =' num2str(sum(WforConnectedCells)) ]);
    subplot(132); hist(WforConnectedCells_u); title(['Uniform, sum W =' num2str(sum(WforConnectedCells_u)) ]);
    subplot(133); hist(WforConnectedCells_e); title(['Exponential, sum W =' num2str(sum(WforConnectedCells_e)) ]);
    suptitle('Distribution of weight in one layer2 cell (one cell example)');
end

end

if(PLOT_FIG)
% Connectivity
figure; set(gcf,'position',[   657   384   712   420]);
subplot(121); hist(SumConnectivity_g); title(['Gaussian, mean =' num2str(mean(SumConnectivity_g)) ]); 
subplot(122); hist(SumConnectivity_s); title(['Sparse, mean =' num2str(mean(SumConnectivity_s)) ]); 
suptitle('Average sum of connectivity of all cells within range (= Volume of Connectivity)');


figure; set(gcf,'position',[  657   384   712   420]);
subplot(121); hist(NumConnection_g); title(['Gaussian, mean =' num2str(mean(NumConnection_g)) ]); 
subplot(122); hist(NumConnection_s); title(['Sparse, mean =' num2str(mean(NumConnection_s)) ]); 
suptitle('Average Number of Connection');

%Connection Strength
% figure; set(gcf,'position',[  350         380        1369         420]);
% subplot(131); hist(SumTotalConnStrength_g); title(['Gaussian, mean =' num2str(mean(SumTotalConnStrength_g)) ]); 
% subplot(132); hist(SumTotalConnStrength_u); title(['Uniform, mean =' num2str(mean(SumTotalConnStrength_u)) ]); 
% subplot(133); hist(SumTotalConnStrength_e); title(['Exponential, mean =' num2str(mean(SumTotalConnStrength_e)) ]); 
% suptitle('Average Total W of all possible connection within range (= Volume of Connection Strength)');

figure; set(gcf,'position',[  350         380        1369         420]);
subplot(131); hist(SumConnStrength_g); title(['Gaussian, mean =' num2str(mean(SumConnStrength_g)) ]); 
subplot(132); hist(SumConnStrength_u); title(['Uniform, mean =' num2str(mean(SumConnStrength_u)) ]); 
subplot(133); hist(SumConnStrength_e); title(['Exponential, mean =' num2str(mean(SumConnStrength_e)) ]); 
suptitle('Average summation of effective strength ( strength of connected cells)');
end
%%
     
%save
tmpStruct.Conn_Mat = sparse(gConnMat);   tmpStruct.W_Mat = sparse(gW_Mat);  tmpStruct.Delay_Mat = sparse(gDelay_Mat); 
saveConn_g{rr,ww} = tmpStruct;

tmpStruct.Conn_Mat = sparse(uConnMat);   tmpStruct.W_Mat = sparse(uW_Mat);  tmpStruct.Delay_Mat = sparse(uDelay_Mat); 
saveConn_u{rr,ww} = tmpStruct;

tmpStruct.Conn_Mat = sparse(eConnMat);   tmpStruct.W_Mat = sparse(eW_Mat);  tmpStruct.Delay_Mat = sparse(eDelay_Mat); 
saveConn_e{rr,ww} = tmpStruct;


% G-G
[i,j,val] = find(gConnMat);
conn_list = zeros(length(val), 4);
cnt = 0;
for  m = 1: length(i)
        cnt = cnt + 1; 
        conn_list(cnt,:) =[i(m)-1 j(m)-1 gW_Mat(i(m),j(m)) gDelay_Mat(i(m),j(m))];
end
%Save to file
fname = sprintf('ConvergentInput_GG_Wscale%g_W%g_range%g_Trial%g.txt', W_SCALE, weight_factor, range, TRIAL);
fid = fopen(fname,'w');
fprintf( fid,'%d %d \n', size(conn_list,1),  size(gW_Mat,1)); % Total number of Conn , Total number of cells in Layer1 
fprintf( fid,'%d %d %1.9f %g\n', conn_list' ); % SourceID TargetID Weight Delay
fclose(fid);

% S-U
[i,j,val] = find(uConnMat);
conn_list = zeros(length(val), 4);
cnt = 0;
for  m = 1: length(i)
        cnt = cnt + 1; 
        conn_list(cnt,:) =[i(m)-1 j(m)-1 uW_Mat(i(m),j(m)) uDelay_Mat(i(m),j(m))];
end
%Save to file
fname = sprintf('ConvergentInput_SU_Wscale%g_W%g_range%g_Trial%g.txt', W_SCALE, weight_factor, range, TRIAL);
fid = fopen(fname,'w');
fprintf( fid,'%d %d \n', size(conn_list,1),  size(uW_Mat,1)); % Total number of Conn , Total number of cells in Layer1 
fprintf( fid,'%d %d %1.9f %g\n', conn_list' ); % SourceID TargetID Weight Delay
fclose(fid);


% S-E
[i,j,val] = find(eConnMat);
conn_list = zeros(length(val), 4);
cnt = 0;
for  m = 1: length(i)
        cnt = cnt + 1; 
        conn_list(cnt,:) =[i(m)-1 j(m)-1 eW_Mat(i(m),j(m)) eDelay_Mat(i(m),j(m))];
end
%Save to file
fname = sprintf('ConvergentInput_SE_Wscale%g_W%g_range%g_Trial%g.txt', W_SCALE, weight_factor, range, TRIAL);
fid = fopen(fname,'w');
fprintf( fid,'%d %d \n', size(conn_list,1),  size(eW_Mat,1)); % Total number of Conn , Total number of cells in Layer1 
fprintf( fid,'%d %d %1.9f %g\n', conn_list' ); % SourceID TargetID Weight Delay
fclose(fid);

rw_NumCon_g(rr,ww) = mean(NumConnection_g);  rw_NumCon_s(rr,ww) = mean(NumConnection_s); 
rw_sumW_g(rr,ww) = mean(SumConnStrength_g); rw_sumW_u(rr,ww) = mean(SumConnStrength_u);  rw_sumW_e(rr,ww) = mean(SumConnStrength_e); 

% Report
disp('==================================================================================================');
disp(['Done :: range = ' num2str(range) ', weight factor = ' num2str(weight_factor) ])
disp(['Average Number of Connection::: Gaussian ' num2str(mean(NumConnection_g)) ', Sparse' num2str(mean(NumConnection_s)) ])
disp(['Average summation of effective strength::: G-G = ' num2str(mean(SumConnStrength_g)) ', S-U = ' num2str(mean(SumConnStrength_u)) ', S-E = ' num2str(mean(SumConnStrength_e))])
toc(tt)
disp('==================================================================================================');


    end
end

%% SAVE sim data
% save(['NewGennConn' date '.mat'],'saveConn_g','saveConn_u','saveConn_e','W_LST','Range_LST','-v7.3')
