% Load from Mat files 

% Connection Types
NumCnvrgntTypes = 3; CnvrgntTypes = cell(NumCnvrgntTypes,1);
CnvrgntTypes{1}.FileName = 'GPmVLmd1_0del_KO2'; 
CnvrgntTypes{2}.FileName = 'GPmVLmd1_0del_KO2_avgPWuniformTC'; 
CnvrgntTypes{3}.FileName = 'GPmVLmd1_0del_KO2_avgPnegExpWTC' ; 
CnvrgntTypes{1}.TitleName = 'Connectivity: Gaussian, Strength: Gaussian'; 
CnvrgntTypes{2}.TitleName = 'Connectivity: Uniform, Strength: Uniform'; 
CnvrgntTypes{3}.TitleName = 'Connectivity: Uniform, Strength: Negative exponential'; 
CnvrgntTypes{1}.CodeName = 'gaussPgaussW'; 
CnvrgntTypes{2}.CodeName = 'avgPuniformW'; 
CnvrgntTypes{3}.CodeName = 'avgPnegecpW'; 
CnvrgntTypes{1}.leg = 'P:Gauss, W:Gauss'; 
CnvrgntTypes{2}.leg = 'P:uniform, W:uniform';
CnvrgntTypes{3}.leg = 'P:uniform, W:Neg Exp'; 

    matFile = {'sparseActivity_gaussPgaussW05-Jun-2015','sparseActivity_avgPuniformW05-Jun-2015', 'sparseActivity_avgPnegexpW05-Jun-2015'};

    for ct_ii = 1 : NumCnvrgntTypes
        tload = tic();
        CnvrgntTypes{ct_ii}.matFile = matFile{ct_ii};
        load([matFile{ct_ii} '.mat' ]);
        CnvrgntTypes{ct_ii}.ACT_Record = ACT_Record;
        disp('============================================================================================');
        disp(CnvrgntTypes{ct_ii}.TitleName)
        toc(tload);
        disp('============================================================================================');
        
    end
%%  Analysis 

avgFR_Matrix = cell(NumCnvrgntTypes,2,3); 


    for ct_ii = 1 : NumCnvrgntTypes 
    ACT_Record = CnvrgntTypes{ct_ii}.ACT_Record;
    tmpdirFig =  ['RestartAnalysisOnConnTypes/'];
     if ct_ii == 1 
        dirLoc = [PATH 'OscInput_Sim/'];
        dirFig = ['../OscInput_varyTCtype_Sim/'    tmpdirFig ];    
     else     
        dirLoc = [PATH 'OscInput_varyTCtype_Sim/']; % [PATH 'OscInput_Sim/']; % /OscInput_varyTCtype_Sim/
        dirFig = [   tmpdirFig ];
     end
%      mkdir([dirLoc dirFig])
    
    CtypeTxt =  CnvrgntTypes{ct_ii}.TitleName;
    codeTxt =   CnvrgntTypes{ct_ii}.CodeName;
    
    mat = get_M1FR_matrix(ACT_Record);
    
    
    M1_BaselineActivity_Matrix_noIndvPlot  % ---> plot raw M1 activity 
% M1_BaselineActivity_Matrix_noIndvPlot_normVL % ---> plot M1 activity  - VL activity (normalized with VL activity)
% Diff_M1act_OSC_F_AMP_oneType_normVL

    end