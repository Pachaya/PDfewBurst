function [spktrain, spktime ] = PoissonSpikeTrain( Fr,T )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
 time = 0:1:T; % input time resolution is 1ms
 dfr = Fr/1000;
 rr = rand(T,1);
 time(rr < dfr) = 1;
 spktrain = time;
 spktime = find(time == 1); 
end

