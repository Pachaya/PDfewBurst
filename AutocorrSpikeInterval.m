function [Npeak, width, sumISI] = AutocorrSpikeInterval(spktrainVec, Trange, ttl,cutTime )
% AutoCorrelation from inter-spike-interval (across cells of same type)
% when spiketrain is given

% spktrainVec = WT.E.spktrain(:,T1:T2); Trange = 10;
% spktrainVec = WT.VL.normal.spkMatCut; Trange = 40; ttl = 'WT VL';
% cutTime = 0; 

box = -Trange:1:Trange; %ConsiderResolutionLater
% sumAC = zeros(1,length(box));
sumISI = [];

Upperborder = size(spktrainVec,2) - length(box);

for iii = 1 : size(spktrainVec,1)
    spktime = find(spktrainVec(iii,:) ==1);
    for spk = 1 : length(spktime) % For each spike in reference cell #iii
        if(spktime(spk)-Trange > cutTime) && (spktime(spk)+Trange <Upperborder)
            for id = 1 : size(spktrainVec,1)
                tmpISI = box((spktrainVec(id,box+spktime(spk))==1));
%                 tmpISI = box(tmpISI);
                if(id ~= iii)  && ~isempty(tmpISI)
                    %                 sumAC = sumAC + spktrainVec(id,box+spktime(spk));  
                    for cntISI = 1 : length(tmpISI)
                        sumISI(end + 1) = tmpISI(cntISI);%-spktime(spk);
                    end
                end
            end
        end
    end
end

% [m,s] = normfit(sumISI);
Npeak=hist(sumISI,box); 
normISIpeak = Npeak./Npeak(ceil(length(box)/2)); 
figure; hist(sumISI,box); title(ttl);
xlabel('time (ms)'); ylabel('Number of spikes');
hold on;
% figure;
% dd= -3.5*s:s*3.5;
% plot(dd,max(Npeak)*exp( -dd.^2/2/s^2 ),'r','LineWidth',2);
width = 2.355*s; %width at half of maximum, normal fit
% title([ttl ' , [mu,sig] = [' num2str(m) ',' num2str(s) '], width = ' num2str(width) ]); 
xlabel('time (ms)'); ylabel('Number of spikes');
figure; bar(box,normISIpeak); title([ttl ' Normalized'])
xlabel('time (ms)'); ylabel('Spikes Ratio');

end