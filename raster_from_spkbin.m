function [freq_all,spkTime_all, fg_handle ] = raster_from_spkbin(  spkBinAll,cutTime, Tstop, condition)
%Make a raster plot from input spike bin
%   Detailed explanation goes here



freq_all = zeros(size( spkBinAll,1),1);
spkTime_all =  cell(size( spkBinAll,1),1);
fg_handle = figure;
hold on
for tr = 1 : size( spkBinAll,1)
    spkBin_1tr = spkBinAll(tr,:);
    spkT = find(spkBin_1tr ==1);
    freq = sum(spkT > cutTime) / (Tstop - cutTime)*1000;
    freq_all(tr) = freq;
    spkTime_all{tr}=spkT;
    for jj = 1 : length(spkT)
        line([spkT(jj) spkT(jj)],[tr-1 tr])
    end
end
avgFreq = mean(freq_all);
avgSTD = std(freq_all);
ylim([0  size( spkBinAll,1)+1])
xlim ([0 Tstop])
xlabel('Time (ms)')
title(['Raster Plot ( # trial = ' num2str( size( spkBinAll,1)) ', ' condition  ')' ', Avg Fr: ' num2str(avgFreq) ', S.D.:' num2str(avgSTD) ]) %%% or make another plot for PoisSpk



end

