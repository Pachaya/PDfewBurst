% Make connection for layer 1 - layer 2 with three types of conenction
% rules
% MakeLayer1Layer2connection_3ConnTypes.m

PROP_SPEED = 1000 ; % [um/ms] speed of axonal propagation < 1 m/s >
layer1locTxt = 'Epos_E1150_I380_Regular_VL_Sweep5Itr5000_ceil_sorted';
layer2locTxt = 'ECenterpos_rTC250_E5625_I1892_ssize1500_Regular_sweep5Itr5000_ceil_sorted';

SimCode = 'Theory';
[L1pos, L1Distmat, L1nnDist]= GetCellsDistInfo(layer1locTxt, 1);
[L2pos, L2Distmat, L2nnDist]= GetCellsDistInfo(layer2locTxt, 1); %VL

%%
rng('shuffle')
% Range_LST = [10 25:25:300];
Range_LST = [150];
saveConn = cell(size(Range_LST));
ctrlW = 0.001; 
DistMat = pdist2(L1pos,L2pos);

W_SCALE = 0.00001; 
W_LST = [100]; % 50 : 25 : 100
% Three Kinds of Connecting Rules 
% First Gaussian
% Then Uniform  - same posibility of getting W
%      and NegExp - w_sampling = exprnd(mu of W);
ww = 1;
for rr = 1 : length(Range_LST)
tt = tic();
range = Range_LST(rr);
range = 100;
weight_factor = W_LST(ww);
sig = range / sqrt(2);


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
    potentialLayer1 = find(DistMat(:,l2_ID) <= sig*3 ); % Cells that located within range
    
    %%% Convergent Connection Rule 1
    disp('--------------------- G-G ---------------------');
    % ############## Connectivity : Gaussian Distribution
    pp = func_gauss3(Pmax, DistMat(potentialLayer1,l2_ID),sig);
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
    tmpW = func_gauss3(weight_factor*W_SCALE, DistMat(potentialLayer1, l2_ID),sig); %for all possible connection within range 
    WforConnectedCells =  tmpW(tmpVec_g);    
    
    %Volume Comparison
    SumTotalConnStrength_g(l2_ID) = sum(tmpW);
    SumConnStrength_g(l2_ID) =  sum(WforConnectedCells);
    
    % Recorded Connection Strength
    gW_Mat(tmpConnID,l2_ID) = WforConnectedCells;
    
    
    disp('--------------------- Sparse Connectivity in S-U and S-E ---------------------');
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
    disp('-----------------------S-U ---------------------');

    % ############## Connection Strength : Uniform Distribution
    wGauss = sum(func_gauss3(weight_factor*W_SCALE, DistMat(potentialLayer1,l2_ID),sig)); % The possible weight in Gaussian 
    avgW = wGauss ./ length(potentialLayer1); 
    tmpW = ones(size(potentialLayer1)).*avgW;
    WforConnectedCells_u = tmpW( tmpVec_s);  % For Uniform distribution
    
    %Volume Comparison
    SumTotalConnStrength_u(l2_ID) = sum(tmpW);
    SumConnStrength_u(l2_ID) =  sum(WforConnectedCells_u);
    
    % Recorded Connection Strength
    gW_Mat(tmpConnID_s,l2_ID) = WforConnectedCells_u;
    
    
    %%% Convergent Connection Rule 3 : Connectivity - Uniform , Strength -
    %%% Negative Exponential
    disp('-----------------------S-E ---------------------');

    % ############## Connection Strength : Exponential Distribution
    wGauss = sum(func_gauss3(weight_factor*W_SCALE, DistMat(potentialLayer1,l2_ID),sig)); % The possible weight in Gaussian 
    avgW = wGauss./ length(potentialLayer1); 
    tmpW = exprnd(avgW, size(potentialLayer1));
    WforConnectedCells_e = tmpW( tmpVec_s);  % For Uniform distribution
    
    %Volume Comparison
    SumTotalConnStrength_e(l2_ID) = sum(tmpW);
    SumConnStrength_e(l2_ID) =  sum(WforConnectedCells_e);
    
    % Recorded Connection Strength
    gW_Mat(tmpConnID_s,l2_ID) = WforConnectedCells_e;  
    
end



% Connectivity
figure; set(gcf,'position',[  350         380        1369         420]);
subplot(121); hist(SumConnectivity_g); title(['Gaussian, mean =' num2str(mean(SumConnectivity_g)) ]); 
subplot(122); hist(SumConnectivity_s); title(['Sparse, mean =' num2str(mean(SumConnectivity_s)) ]); 
suptitle('Average Total possible connectivity (= Volume of Connectivity)');


figure; set(gcf,'position',[  350         380        1369         420]);
subplot(121); hist(NumConnection_g); title(['Gaussian, mean =' num2str(mean(NumConnection_g)) ]); 
subplot(122); hist(NumConnection_s); title(['Sparse, mean =' num2str(mean(NumConnection_s)) ]); 
suptitle('Average Number of Connection');

%Connection Strength
figure; set(gcf,'position',[  350         380        1369         420]);
subplot(131); hist(SumTotalConnStrength_g); title(['Gaussian, mean =' num2str(mean(SumTotalConnStrength_g)) ]); 
subplot(132); hist(SumTotalConnStrength_u); title(['Uniform, mean =' num2str(mean(SumTotalConnStrength_u)) ]); 
subplot(133); hist(SumTotalConnStrength_e); title(['Exponential, mean =' num2str(mean(SumTotalConnStrength_e)) ]); 
suptitle('Average Total possible connectivity (= Volume of Connectivity)');


figure; set(gcf,'position',[  350         380        1369         420]);
subplot(131); hist(SumConnStrength_g); title(['Gaussian, mean =' num2str(mean(SumConnStrength_g)) ]); 
subplot(132); hist(SumConnStrength_u); title(['Uniform, mean =' num2str(mean(SumConnStrength_u)) ]); 
subplot(133); hist(SumConnStrength_e); title(['Exponential, mean =' num2str(mean(SumConnStrength_e)) ]); 
suptitle('Average summation of effective strength ( strength of connected cells)');
%%
%      
% %save 
% saveConn{rr} = sparse(gW_Mat);
% 
% tmptmpW  = sum(gW_Mat,2);
% IDID = find(tmptmpW  >= 10^(-10));
% cutW_mat =  gW_Mat(IDID, :); %minimum spike generator required for this type of connection
% 
% 
% [i,j,val] = find(cutW_mat);
% conn_list = [i,j,val]; % source (layer 0) , target (layer1),  weight    
% %/// Can convert to sparse matrix via data = spconvert(conn_list );  or
% %    full matrix via W_Mat = full(spconvert(conn_list );
% %Save to file
% fname = sprintf('CorrelatedPoissonInput_ctrlW%g_range%g.txt', ctrlW,range);
% fid = fopen(fname,'w');
% fprintf( fid,'%d %d \n', size(conn_list,1),  size(cutW_mat,1)); % Total number of Conn , Total number of cell to b e ui
% fprintf( fid,'%d %d %1.9f\n', conn_list' );
% fclose(fid);
% 
% 
% % Report
% disp('==================================================================================================');
% fprintf('Range = %g\t\tSigma = %g\n', range, sig);
% disp(['avg #layer0 / layer1  = ' num2str(mean(sum(gConnMat)))]);
% disp(['avg sumW = ' num2str(mean(sum(gW_Mat)))]);
% toc(tt)
% disp('==================================================================================================');

end
