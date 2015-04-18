function [ tmp_MS,segment_fr_data_all ] = get_mean_std_fr_of_segment_fromSpkTrain( tmp_spk, binT, binsize )
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here

tmp_MS = zeros(size(tmp_spk,1),2); %mean std
segment_fr_data_all = zeros(size(tmp_spk,1), length(binT)-1);
for id = 1 : size(tmp_spk,1)
% binsize = 200; 
% binT = cuttime:binsize:Tstop;
segment_fr_data = zeros(length(binT)-1,1);
for ii = 1:length(binT)-1
    segment_fr_data(ii) = sum(tmp_spk(id,binT(ii):binT(ii+1)))./binsize*1000;
end
    tmp_MS(id,:) =[mean(segment_fr_data) std(segment_fr_data)];
    segment_fr_data_all(id,:) =segment_fr_data;
end


end

