function [box, sumAC] = CrosscorrFromSpktrain_VL_M1_Return (spktrainVecVL, spktrainVecM1 , Trange,cutTime )
% AutoCorrelation from Spktrain

% spktrainVec = WT.All.spktrain;  Trange = 50;
% spktrainVec = WT.VL.burst.spkMatCut; Trange = 40; ttl = 'WT VL';
% size(spktrainVec)
box = -Trange:1:Trange; %ConsiderResolutionLater
sumAC = zeros(1,length(box));
 
Upperborder = size(spktrainVecVL,2) - length(sumAC);
for id = 1 : size(spktrainVecVL,1)
    spktime = find(spktrainVecVL(id,:) ==1);
    for spk = 1 : length(spktime)
        for id_m1 = 1 : size(spktrainVecM1,1)
            if(spktime(spk)-Trange > cutTime) && (spktime(spk)+Trange <Upperborder)
                sumAC = sumAC + spktrainVecM1(id_m1,box+spktime(spk));
            end
        end
    end
end

% figure; bar(box,sumAC); title(ttl)



end