% Get only average FR in 10 Trials 

%% Load
N_TRIAL = 10; 
nr = 3; nc = 3; 
saveMat = cell(N_TRIAL, nr,nc);
NumCnvrgntTypes = 3;
tstop =5000;


for TRIAL_NO = 1 : N_TRIAL
% TRIAL_NO = 1 ; 
dirMat = 'AllTrials/';
load([ dirMat 'CnvrgntTypes_20Hz__rTC_50_100_wmTC_50_100_trial_' num2str(TRIAL_NO)])

matsize = [length(PARAMETERS{1}.PARAM) length(PARAMETERS{2}.PARAM)];
cnt = 0; 
p3=1;p4=1;
for ct_ii = 1 : NumCnvrgntTypes
   ACT_Record = CnvrgntTypes{ct_ii}.ACT_Record;
   for p5 = 1 : length(PARAMETERS{5}.PARAM)  %Synchronization Level
       
   tmpMat  = zeros(matsize);
   for ii = 1 : length(PARAMETERS{1}.PARAM)
       for jj = 2 : length(PARAMETERS{2}.PARAM)
           tmpMat(ii,jj) = mean(sum (ACT_Record{ii,jj,p3,p4,p5}.L2.spktrain,2)./tstop*1000);
       end
   end
   cnt = cnt +1;                
   saveMat{TRIAL_NO, p5, ct_ii} = tmpMat;
   end
end
disp(['End TRIAL#' num2str(TRIAL_NO)])

end
%% Change to Matrix 

FRmat_allTrials = cell(nr,nc);
for ct_ii = 1: 3
    for p5 = 1 : 3        
        tmptmpMat = zeros(length(PARAMETERS{1}.PARAM), length(PARAMETERS{2}.PARAM), N_TRIAL);
        for ii  = 1 : N_TRIAL
            tmptmpMat(:,:, ii) = saveMat{ii, p5, ct_ii};
            FRmat_allTrials{p5,ct_ii} = tmptmpMat;
        end
    end
end
PARAMETERS{1}.titleTxt = 'Range_c_o_n';
PARAMETERS{2}.titleTxt = 'W_c_o_n';
dirMat = 'AllTrials/';

CnvrgntTypesINFO = CnvrgntTypes;
for  ii = 1 : 3   
   field = 'ACT_Record';
   CnvrgntTypesINFO{ii} = rmfield(CnvrgntTypesINFO{ii},field);
end
save([dirMat 'FRmatrix_L2_rowSynchLvl_colCnvrgntTypes_Trial1_10_' date '.mat' ],'FRmat_allTrials','PARAMETERS','CnvrgntTypesINFO','-v7.3');



