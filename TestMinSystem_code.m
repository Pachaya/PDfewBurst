   





dirLoc = 'TestMinSystem/';
dirFig = 'Fig/';
mkdir([dirLoc dirFig])




InputFR_LST = [10:10:90];
         
            
  ACT_Record = cell(length(InputFR_LST),1);
  WT = []; KO = []; ACT =[]; ACT.WT =[]; ACT.KO =[];
  
  
for inF_ii = 1 : length(InputFR_LST)
            %SomaVolt_WT_M1_Nsample20_SynchLvl1_TSTOP1200_InputFR245_wSPK0.032_wVLM1_0.0003
            
            
            W_SPK = 0.032;
            W_VL_M1 = 0.0003;
            SynchLvl = 0; 
            InputFR = InputFR_LST(inF_ii);
            CUTTIME = 500;
            Nsample = 20; 
            TSTOP =  200500; 
            coreName ='';
            simCode =[ coreName 'Nsample' num2str(Nsample) '_SynchLvl' num2str(SynchLvl) '_TSTOP' num2str(TSTOP) '_InputFR' num2str(InputFR) ...
                '_wSPK' num2str(W_SPK) '_wVLM1_' num2str(W_VL_M1)];
            disp('=========================================================================');
            disp(['Input FR = ' num2str(InputFR)]);
            
            %WT_VL
            Name_postfix  = ['WT_VL_' simCode ];
            VL_BaselineAct
            tmp.E = E;      tmp.I = I;  tmp.All = All; tmp.All.somaVall = somaVall;
            WT.VL = tmp;
            
            %KO_VL
            Name_postfix  = ['KO_VL_' simCode ];
            VL_BaselineAct
            tmp.E = E;      tmp.I = I;  tmp.All = All; tmp.All.somaVall = somaVall;
            KO.VL = tmp;           
            
            ACT.WT =  WT; ACT.KO = KO;
            
            ACT_Record{inF_ii} = ACT;

% %             disp( spkBin = ACT_Record{p1_ii,p2_ii,p3_ii}.WT.VL.All.spktrain)
            
end

  MeanOutWT  = zeros( length(InputFR_LST),1);
  MeanOutKO  = zeros( length(InputFR_LST),1);
  
for inF_ii = 1 : length(InputFR_LST)
    
MeanOutWT(inF_ii) = mean(ACT_Record{inF_ii}.WT.VL.All.fr_data);
MeanOutKO(inF_ii) = mean(ACT_Record{inF_ii}.KO.VL.All.fr_data);

end

 
WTWT = [26.8600   17.2000    3.6650    5.9900    5.6200    1.1600   19.3700   11.8000   28.4000    5.2250  22.9950   14.7750    4.7800   12.1750    6.3200   14.9650   18.1800   21.6300    7.8700    6.9200   4.1300    9.1700    4.4150    5.7400    7.0950   10.4250    3.4350    9.6450    8.7300    6.9300   13.2963   11.5820   32.7672   11.9365    4.8350    0.9550];
KOKO = [ 7.0750   10.6050    3.5100    1.1800    8.9450   12.4250   12.6800   12.0200   16.2300   25.3400   17.1500    6.4350   19.0950    4.3050   21.9150   26.8100   15.0150    8.8100];

WT_Expected_input = 11.43.*exp(0.06084.*WTWT);   %%%% Fitting Equation
KO_Expected_input = 11.43.*exp(0.06084.*KOKO);

xx = 0:1:35;
yfit = 11.43.*exp(0.06084.*xx);
figure; plot(xx,yfit,'.-k'); hold on;
scatter(WTWT,WT_Expected_input,'xb');
scatter(KOKO,KO_Expected_input,'xr');
scatter(MeanOutWT,InputFR_LST,'b');
scatter(MeanOutKO,InputFR_LST,'r');
xlabel('Output'); ylabel('Input');
title([ coreName ' : Output - Input Function'])

%% save
save([dirLoc dirFig 'Save_run_result_' coreName '200s' '.m'], 'ACT_Record', 'InputFR_LST')
