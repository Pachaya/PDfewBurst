function [freq_all,spkTime_all, fg_handle ] = raster_from_spkbin_BurstRange(  spkBinAll,cutTime, Tstop,PhotoInjT, PhotoStop, DelayT, BurstRange, condition)
%Make a raster plot from input spike bin
%   

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
title({['Raster Plot ( ' condition  ' )'], ['Avg Fr: ' num2str(avgFreq) ', S.D.:' num2str(avgSTD) ]}) %%% or make another plot for PoisSpk

burstT1 =PhotoStop +DelayT;  burstT2 =PhotoStop +DelayT+BurstRange; 
totalCells =  size( spkBinAll,1);
rectangle('Position',[PhotoInjT,0,PhotoStop-PhotoInjT,totalCells],'EdgeColor', 'k')
rectangle('Position',[burstT1,0,burstT2-burstT1,totalCells],'EdgeColor', 'r')

end

