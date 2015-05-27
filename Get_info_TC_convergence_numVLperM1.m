%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get info for TC convergenc --> number of VL per M1 , average maxW , average summation of weight
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p3_ii = 1; p4_ii = 1; p5_ii =1;
ssize = 1500; dist = ssize/2; pmax =0.85; W_scale = 1E-05;
fprintf('RangeTC\t\t wmTC\t\t numVL/M1\t\t sumW\t\t maxW\t\n');
numVL_M1_LST = zeros(length(PARAM1),length(PARAM2));
sumW_LST = 		 zeros(length(PARAM1),length(PARAM2));
maxW_LST  = zeros(length(PARAM1),length(PARAM2));

for p1_ii = 1 : length(PARAM1)
    for p2_ii = 1 : length(PARAM2)
%         disp('##################################################################################################')
%         disp( [lblTxt1 ' = ' num2str(PARAM1(p1_ii)) ]);
%         disp( [lblTxt2 ' = ' num2str(PARAM2(p2_ii)) ]);
        rangeTC = PARAM1(p1_ii);      sigTC = rangeTC/sqrt(2);
        wmTC = PARAM2(p2_ii);
        
        Name_postfix_WT  = ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii }.VL.WT.Simulation_Code;  % The connection is same for WT and KO
        Name_postfix_KO  = ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii }.VL.KO.Simulation_Code;
        display = 0;
        [TC_basedOnM1_WT, TC_sumW_WT, TC_maxW_WT, TC_numVL_WT ]               = ExtractTC_info( dirLoc, Name_postfix_WT,  display);
%         disp([ num2str(rangeTC) '        ' num2str(mean(TC_numVL_WT)) '        ' num2str(specified_wmTC) '        ' num2str(mean(TC_sumW_WT)) '        ' num2str(mean(TC_maxW_WT))  ])
        fprintf(' %3.0f\t\t%7.4f\t\t%7.4X\t\t%7.4X\t\t%7.4X\n', rangeTC, wmTC, mean(TC_numVL_WT), mean(TC_sumW_WT), mean(TC_maxW_WT));
        numVL_M1_LST(p1_ii,p2_ii) = mean(TC_numVL_WT);
        sumW_LST(p1_ii,p2_ii) = 	mean(TC_sumW_WT);
        maxW_LST(p1_ii,p2_ii)  = mean(TC_maxW_WT);
    end
end
VLperM1 = numVL_M1_LST(:,1);
