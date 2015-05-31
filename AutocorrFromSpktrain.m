function [box,sumAC, cntSample] = AutocorrFromSpktrain(spktrainVec, Trange, cutTime )
% AutoCorrelation from Spktrain
%   - Cut edges data box out 
%   - Don't count reference spike
%   - don't count other spikes in the ref ID
%   - Data =  all other cells ID ( not include the reference ID)
%   - Change the ref cell to other cells

bin = 1; % 1ms
box = -Trange:bin:Trange; 
sumAC = zeros(1,length(box));
Upperborder = size(spktrainVec,2) - length(sumAC);
cntSample =0;
for ref_id = 1 : size(spktrainVec,1)
    spktime = find(spktrainVec(ref_id,:) ==1);
    for spk = 1 : length(spktime)
        for data_id = 1 : size(spktrainVec,1)
            if(data_id ~= ref_id)
                if(spktime(spk)-Trange > cutTime) && (spktime(spk)+Trange <Upperborder)
                    sumAC = sumAC + spktrainVec(data_id,box+spktime(spk));
                    cntSample = cntSample +1;
                end
            end
        end
    end
end

% figure; bar(box,sumAC); title(ttl)

end
