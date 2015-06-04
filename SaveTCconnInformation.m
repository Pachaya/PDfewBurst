% Get Interlayer Connection Information


p3_ii = 1; p6_ii = 1;

for ct_ii = 1 : NumCnvrgntTypes 
 ACT_Record = CnvrgntTypes{ct_ii}.ACT_Record;
     tmpdirFig =  ['CompareConTypes_' num2str( PoisInputFr) 'Hz_' get_Parameters_RangeTxt( PARAMETERS,[1,2,4,5,6]) '/'];
%   tmpdirFig =  ['RestartAnalysisOnConnTypes/'];
 if ct_ii == 1 
    dirLoc = [PATH 'OscInput_Sim/'];
    dirFig = ['../OscInput_varyTCtype_Sim/'    tmpdirFig ];    
 else     
    dirLoc = [PATH 'OscInput_varyTCtype_Sim/']; % [PATH 'OscInput_Sim/']; % /OscInput_varyTCtype_Sim/
    dirFig = [   tmpdirFig ];
 end
 TC_Conn_INFO = cell(ACT_Rec_size);
 
 for p1_ii = 1 : length(PARAM1)
     for p2_ii = 1 : length(PARAM2)
         for p4_ii = 1 : length(PARAM4)
             for p5_ii = 1 : length(PARAM5)
      
   
     
                 Name_postfix_WT  = ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii }.VL.WT.Simulation_Code;  % The connection is same for WT and KO
                 Name_postfix_KO  = ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii }.VL.KO.Simulation_Code;
        
        
                 [TC_basedOnM1, TCconn_raw ]  = ExtractTC_info( dirLoc, Name_postfix_WT, 0  );
                 
                 TC_Conn_INFO{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii,p6_ii}.TC_basedOnM1  = TC_basedOnM1;
                 TC_Conn_INFO{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii,p6_ii}.TCconn_raw = TCconn_raw;
             end
         end
     end
 end
 CnvrgntType{ct_ii}.TC_Conn_INFO = TC_Conn_INFO; 
 
 % save 
 
    contype =   CnvrgntTypes{ct_ii}.CodeName;

    save([dirLoc dirFig 'TC_Conn_INFO_' contype '_' date '.mat' ], 'TC_Conn_INFO','PARAMETERS','-v7.3');

    disp('============================================================================================');
    disp('                      Saved interlayer connection information;');
    disp('============================================================================================');
    disp(CnvrgntTypes{ct_ii}.TitleName)
    disp('============================================================================================');
    disp('============================================================================================');
%     clear TC_Conn_INFO
end

 save([dirLoc dirFig 'ActivityResult_allCnvrgntTypes_' date '.mat' ], 'CnvrgntTypes','PARAMETERS','-v7.3');