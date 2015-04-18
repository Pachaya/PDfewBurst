function [N_E,N_I, Tstop,allVec] = get_save_vec_file_BigNet(fname )
%Read saveed Vector files and return number of trials(N_Trial), run
%time(Tstop), and array of the saved data (allVec)
%   .
% fname = [dirLoc  'SomaVolt_' 'TestSim_WT_InGauss0.2_W0.0006_50.00Hz_rEE0.0625_rEI0.5_rIE5_rII5_Wmult500.txt'];
aa = importdata(fname);
N_E = aa(1,1);
N_I = aa(1,2);
N_Trial = N_E + N_I;
Tstop = aa(1,3);
Trial_List = aa(2,:);
aa = aa(3:end,:);
allVec = zeros(N_Trial,size(aa,1));
for i = 1:length(Trial_List)
    allVec(i,:) = aa(:,i)';
end
if(N_Trial < 3)
    tmp = allVec;
    allVec = zeros(N_Trial,size(aa,1));
    for i = 1:N_Trial
        allVec(i,:) = tmp(i,:);
    end
end
        
    
end
