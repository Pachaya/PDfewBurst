

%% Making population of poisson spike

N = 20; 
SynLevel = 1; %% 0 - 1
Fr = 40; %Hz
T = 2000; %2s
shiftDel = 100;
[spktrain, spktime ] = PoissonSpikeTrain( Fr,T );
spkMat = zeros(T,N);


% Test Full Synchronization
spkMat = repmat(spktrain,N,1);

[freq_all,spkTime_all] = raster_from_spkbin_noTitle(  spkMat,500, T, ['Test Full']);

