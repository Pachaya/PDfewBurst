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
% 
% ConnMat = DistMat.*0;
% W_Mat = DistMat.*0;
% Pmax = 0.85; 
% 
% % connectivity
% pp = func_gauss3(Pmax, DistMat,sig);
% xx = rand(size(pp));
% tmpVec =  xx < pp;
% ConnMat(tmpVec) = 1; % 1 = exist connection, 0 = no connection;
% % Connection Strength -->  sum W = 0.001
% tmpW = func_gauss3(Pmax, DistMat,sig);
% tmpW = tmpW.*tmpVec;
% W_Mat = tmpW.*ctrlW./(repmat(sum(tmpW), length(L0pos),1)); % normalized, so that the sum W = control W = 0.001


%%
gConnMat = DistMat.*0;  uConnMat = DistMat.*0;  eConnMat = DistMat.*0;
gW_Mat = DistMat.*0;    uW_Mat = DistMat.*0;    eW_Mat = DistMat.*0;
gDelay_Mat = DistMat.*0;    uDelay_Mat = DistMat.*0;    eDelay_Mat = DistMat.*0;
Pmax = 0.85;

SumConnectivity_g =  zeros( length(L2pos), 1);
NumConnection_g = zeros( length(L2pos), 1);

SumConnectivity_u =  zeros( length(L2pos), 1);
NumConnection_u = zeros( length(L2pos), 1);

for l2_ID = 1 : length(L2pos) % for each cell in layer 2
    %%% Convergent Connection Rule 1
    disp('-------------------------- G-G ---------------------');
    % connectivity --- Determine
    pp = func_gauss3(Pmax, DistMat(:,l2_ID),sig);
    xx = rand(size(pp));
    tmpVec =  xx < pp; % sumProb = sum(pp(tmpVec));   disp(['Summation of connectivity : ' num2str(sumProb)]);
     sumProb = sum(pp);   disp(['Summation of connectivity : ' num2str(sumProb)]);
    gConnMat(tmpVec,l2_ID) = 1; % 1 = exist connection, 0 = no connection;
    nCon_of_g = sum(tmpVec); disp(['Number of Connection in G : ' num2str(nCon_of_g)])
    gDelay_Mat(tmpVec,l2_ID) = DistMat(tmpVec, l2_ID)./PROP_SPEED;
    
    SumConnectivity_g(l2_ID) =  sumProb;
    NumConnection_g(l2_ID) = nCon_of_g;
    % Connection Strength
    tmpW = func_gauss3(weight_factor*W_SCALE, DistMat(tmpVec,l2_ID),sig); % Only the cells that are connected
    sumW_of_g = sum(tmpW); disp(['Summation of connection strength : ' num2str(sumW_of_g)]);
    gW_Mat(tmpVec,l2_ID) = tmpW;
    
    %%% Convergent Connection Rule 2 : Connectivity - Uniform , Strength - Uniform
    disp('-------------------------- S-U ---------------------');
    % only the cells in range
    potentialLayer1 = find(DistMat(:,l2_ID) <= sig*3 );
    %connectivity
    pGauss  = sum(func_gauss3(Pmax, DistMat(potentialLayer1,l2_ID),sig));
    pp2 = ones(size(potentialLayer1)).*(pGauss / length(potentialLayer1)); % Average weight per one L1 cell
    xx = rand(size(pp2));
    tmpVec2 =  xx < pp2;  %sumProb2 = sum(pp2(tmpVec2));    disp(['Summation of connectivity : ' num2str(sumProb2)]);
    sumProb2 = sum(pp2);    disp(['Summation of connectivity : ' num2str(sumProb2)]);
    nCon_of_u = sum(tmpVec2);
    disp(['Number of Connection in U : ' num2str(nCon_of_u)])
    SumConnectivity_u(l2_ID) =  sumProb2;
    NumConnection_u(l2_ID) = nCon_of_u;
    %  disp(['Summation of connection strength : ' num2str(sumW_of_g)]);
    
    
end

figure; 
subplot(121); hist(SumConnectivity_g); title(['G-G, mean =' num2str(mean(SumConnectivity_g)) ]); 
subplot(122); hist(SumConnectivity_u); title(['S-U, mean =' num2str(mean(SumConnectivity_u)) ]); 
suptitle('Summation of connected connectivity');


figure; 
subplot(121); hist(NumConnection_g); title(['G-G, mean =' num2str(mean(NumConnection_g)) ]); 
subplot(122); hist(NumConnection_u); title(['S-U, mean =' num2str(mean(NumConnection_u)) ]); 
suptitle('Number of Connection');
%%
     
%save 
saveConn{rr} = sparse(gW_Mat);

tmptmpW  = sum(gW_Mat,2);
IDID = find(tmptmpW  >= 10^(-10));
cutW_mat =  gW_Mat(IDID, :); %minimum spike generator required for this type of connection


[i,j,val] = find(cutW_mat);
conn_list = [i,j,val]; % source (layer 0) , target (layer1),  weight    
%/// Can convert to sparse matrix via data = spconvert(conn_list );  or
%    full matrix via W_Mat = full(spconvert(conn_list );
%Save to file
fname = sprintf('CorrelatedPoissonInput_ctrlW%g_range%g.txt', ctrlW,range);
fid = fopen(fname,'w');
fprintf( fid,'%d %d \n', size(conn_list,1),  size(cutW_mat,1)); % Total number of Conn , Total number of cell to b e ui
fprintf( fid,'%d %d %1.9f\n', conn_list' );
fclose(fid);


% Report
disp('==================================================================================================');
fprintf('Range = %g\t\tSigma = %g\n', range, sig);
disp(['avg #layer0 / layer1  = ' num2str(mean(sum(gConnMat)))]);
disp(['avg sumW = ' num2str(mean(sum(gW_Mat)))]);
toc(tt)
disp('==================================================================================================');

end
