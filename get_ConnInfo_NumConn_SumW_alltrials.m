

dirMat = 'AllTrials/';
load([ dirMat 'FRmatrix_L2_rowSynchLvl_colCnvrgntTypes_Trial1_10_30-Jun-2015']);
 CnvrgntTypesINFO{2}.CodeName = 'UU';  CnvrgntTypesINFO{3}.CodeName = 'UE'; 
 
%%  Get Number of Cells 
load([dirMat 'NewGennConn_r50_100_w_50_100_trial1']);

NumSrc_Matrix = cell{3,1}; 

for ct_ii = 1 : 3
    switch ct_ii 
        case 1 
            tmp = saveConn_g;
        case 2
            tmp = saveConn_u;
        case 3 
            tmp = saveConn_e;
    end
    tmpMat = zeros(6,6);
    tmpMat_w = zeros(6,6);
    for rr = 1 : 6
        for ww = 1:6 
            ff = full(tmp{rr,ww}.Conn_Mat); 
            tmpMat(rr,ww) = mean(sum(ff));
            ff = full(tmp{rr,ww}.Conn_Mat); 
            tmpMat(rr,ww) = mean(sum(ff));
            