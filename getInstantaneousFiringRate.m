function [ sAll ] = getInstantaneousFiringRate(spkTrain, sigma, RES)
%Return the instantenous firing rate with gaussian windowing of the specify SpkBin and sigma 
Tstop = size(spkTrain,2); %ms
ncells = size(spkTrain,1);
        edges =[-3*sigma:1:3*sigma];
        kernel = normpdf(edges,0,sigma);
        kernel = kernel*RES; % Normalized with bin width
        sAll = conv(sum(spkTrain(1:end,:)),kernel);
        center = ceil(length(edges)/2);
        sAll = sAll(center:Tstop+center-1);
        sAll = sAll/ncells*1000;

end

