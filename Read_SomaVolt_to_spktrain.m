% Soma Volt to SPKtrain 
%% Download Membrane potential of the network
THRESHOLD = 0; RES = 1;
[nL2, nI, Tstop,somaVall] = get_save_vec_file_BigNet([dirLoc 'SomaVolt_' SimCode '.txt'] );
somaVallE = somaVall(1:nL2,:);

spktrain = zeros(nL2, tstop);

for ii = 1:nL2
    tmptrain = zeros(Tstop, 1);
    [pks,loc] = findpeaks(somaVallE(ii,:),'MINPEAKHEIGHT',THRESHOLD);
    tmptrain(round(loc*RES)) = 1; 
    spktrain(ii,:) = tmptrain;    
end
spktrain = sparse(spktrain);