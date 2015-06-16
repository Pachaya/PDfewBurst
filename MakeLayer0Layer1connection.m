% Make connection for layer 0 - layer 1 


layer0locTxt = 'Epos4600_noI_ssize3000_repmatVLpos';
layer1locTxt = 'Epos_E1150_I380_Regular_VL_Sweep5Itr5000_ceil_sorted';


[L0pos, L0Distmat, L0nnDist]= GetCellsDistInfo(layer0locTxt, 1);
[L1pos, L1Distmat, L1nnDist]= GetCellsDistInfo(layer1locTxt, 1); %VL

%%
rng('shuffle')
Range_LST = [10 25:25:300];
saveConn = cell(size(Range_LST));
ctrlW = 0.001; 
DistMat = pdist2(L0pos,L1pos);

for rr = 1 : length(Range_LST)
tt = tic();
range = Range_LST(rr);
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
ConnMat = DistMat.*0;
W_Mat = DistMat.*0;

for l1_ID = 1 : length(L1pos)


% connectivity
pp = func_gauss3(Pmax, DistMat(:,l1_ID),sig);
xx = rand(size(pp));
tmpVec =  xx < pp;
ConnMat(tmpVec,l1_ID) = 1; % 1 = exist connection, 0 = no connection;
% Connection Strength -->  sum W = 0.001
tmpW = func_gauss3(1, DistMat(tmpVec,l1_ID),sig);
normW = tmpW./sum(tmpW).*ctrlW; % normalized, so that the sum W = control W = 0.001
W_Mat(tmpVec,l1_ID) = normW;
end

%%
     
%save 
saveConn{rr} = sparse(W_Mat);

tmptmpW  = sum(W_Mat,2);
IDID = find(tmptmpW  >= 10^(-10));
cutW_mat =  W_Mat(IDID, :); %minimum spike generator required for this type of connection


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
disp(['avg #layer0 / layer1  = ' num2str(mean(sum(ConnMat)))]);
disp(['avg sumW = ' num2str(mean(sum(W_Mat)))]);
toc(tt)
disp('==================================================================================================');

end
