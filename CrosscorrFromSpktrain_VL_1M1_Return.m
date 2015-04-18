function [box, sumAC] = CrosscorrFromSpktrain_VL_1M1_Return (spktrainVec1, spktrainM1 , Trange,cutTime )
% AutoCorrelation from Spktrain

% spktrainVec = WT.All.spktrain;  Trange = 50;
% spktrainVec = WT.VL.burst.spkMatCut; Trange = 40; ttl = 'WT VL';
% size(spktrainVec)
box = -Trange:1:Trange; %ConsiderResolutionLater
sumAC = zeros(1,length(box));
 
Upperborder = size(spktrainVec1,2) - length(sumAC);
for id = 1 : size(spktrainVec1,1)
    spktime = find(spktrainVec1(id,:) ==1);
    for spk = 1 : length(spktime)
%         if(spktime(spk)> cutTime) && (spktime(spk) <Upperborder)
        if(spktime(spk)-Trange > cutTime) && (spktime(spk)+Trange <Upperborder)
            sumAC = sumAC + spktrainM1(box+spktime(spk));
        end
    end
end

% figure; bar(box,sumAC); title(ttl)



end