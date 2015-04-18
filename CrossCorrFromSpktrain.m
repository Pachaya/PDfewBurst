function [box,sumAC, cntSample] = CrossCorrFromSpktrain(spktrainVec1,spktrainVec2, Trange, cutTime )
% AutoCorrelation from Spktrain
%   - Cut edges data box out 
%   - Don't count reference spike
%   - don't count outher spikes in the ref ID
%   - Data =  all other cells ID ( not include the reference ID)
%   - Change the ref cell to other cells

bin = 1; % 1ms
box = -Trange:bin:Trange; 
sumAC = zeros(1,length(box));
Upperborder = size(spktrainVec1,2) - length(sumAC);
cntSample =0;
for ref_id = 1 : size(spktrainVec1,1)
    spktime = find(spktrainVec1(ref_id,:) ==1);
    for spk = 1 : length(spktime)
        for data_id = 1 : size(spktrainVec2,1)
                if(spktime(spk)-Trange > cutTime) && (spktime(spk)+Trange <Upperborder)
                    sumAC = sumAC + spktrainVec2(data_id,box+spktime(spk));
                    cntSample = cntSample +1;
                end
       
        end
    end
end

% figure; bar(box,sumAC); title(ttl)

end
