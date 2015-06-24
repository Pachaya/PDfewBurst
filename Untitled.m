% Soma Volt to SPKtrain 
%% Download Membrane potential of the network
[nL2, nI, Tstop,somaVall] = get_save_vec_file_BigNet([dirLoc 'SomaVolt_' Name_postfix '.txt'] );
somaVallE = somaVall(1:nE,:);
somaVallI = somaVall(nE+1:nE+nI,:);




for ii = 1:nL2
    tmptrain = zeros(Tstop+1, 1);
    [pks,loc] = findpeaks(somaVallE(ii,:),'MINPEAKHEIGHT',THRESHOLD);
    tmptrain(round(loc*RES)) =1;
    E.spktrain(ii,:) = tmptrain;
    All.spktrain(ii,:) = tmptrain;
    fr_cellsE(ii) = sum(E.spktrain(ii,cuttime+1:Tstop))/ length(cuttime+1:Tstop)*1000; 
    fr_cells(ii) = fr_cellsE(ii);
    tmpT = PhotoStop + DelayT+1 : PhotoStop + DelayT+BurstRange;
    BurstSpkETrain(ii,:) = E.spktrain(ii,tmpT);  
    BurstSpkTrain(ii,:) = BurstSpkETrain(ii,:);
    BurstSpkE(ii) = sum(E.spktrain(ii,tmpT));  
    BurstSpk(ii) = BurstSpkE(ii);
    tmpT = PhotoStop + DelayT+BurstRange+ 1  : Tstop;
    fr_AfterLightOffE(ii) = sum(E.spktrain(ii,tmpT))/ length(tmpT)*1000;  % after photoactivation  + burstRange stop 
    fr_AfterLightOff(ii) = fr_AfterLightOffE(ii);
end