function [baseline_fr] = get_avg_baseline(spkBin, cutTime, InjStartT)
% Return average firing rate from cutTime - InjStartT

    cut_SpkBin = spkBin(:,cutTime+1 : InjStartT);
    baseline_fr = sum(cut_SpkBin,2)./length(cutTime+1 : InjStartT) * 1000;
    
end

