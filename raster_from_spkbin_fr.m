function [freq_all,spkTime_all, fg_handle ] = raster_from_spkbin_fr(  spkBinAll,cutTime, Tstop, condition, SIGMA_fr_All)
%Make a raster plot from input spike bin
%   Detailed explanation goes here
RES = 1; 
freq_all = zeros(size( spkBinAll,1),1); 
spkTime_all =  cell(size( spkBinAll,1),1);
fg_handle = figure; subplot(6,1,1:3)
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

%Firing Rate 
%All
sigma =  SIGMA_fr_All; % ms
edges =[-3*sigma:1:3*sigma];
kernel = normpdf(edges,0,sigma);
kernel = kernel*RES; % Normalized with bin width
s_All = conv(sum(spkBinAll(1:end,:)),kernel);
center = ceil(length(edges)/2);
s_All = s_All(center:Tstop+center-1);
s_All = s_All/(size(spkBinAll,1))*1000;

subplot(6,1,5:6)
plot(s_All,'k');  hold on;
% legend('All','Location', 'NorthWest')
title('Average firing rate Fr(t)')
ylabel('Average firing rate (Hz)')
xlabel('Time(ms)')


end

