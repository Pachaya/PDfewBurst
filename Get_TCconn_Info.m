function [ TC_Conn ] = Get_TCconn_Info(  Name_postfix, dirLoc )
%Get information about TC connection 
%   


% get Name_postfix , dirLoc

tmpdat = importdata([dirLoc 'TC_ConnectionWDParam_' Name_postfix '.txt']); 
N_Conn = tmpdat(1); tmpdat = tmpdat(2:end);
srclist = tmpdat(1:5:end); tarlist = tmpdat(2:5:end); 
typelist = tmpdat(3:5:end); wlist = tmpdat(4:5:end); 
delaylist = tmpdat(4:5:end); 

% Count Number of VL that connect to M1  === same target cell
minTar = min(tarlist); maxTar = max(tarlist); 
idid = minTar : maxTar;
NumVL_of_M1 = zeros(length(minTar :maxTar),1);
for tarID = 1 : length(idid)
    NumVL_of_M1(tarID) = sum(tarlist == idid(tarID));
end
disp(['Average Number of VL per M1 is ' num2str(mean(NumVL_of_M1))]);
TC_Conn.N = N_Conn; TC_Conn.srclist = srclist; TC_Conn.tarlist = tarlist; 
TC_Conn.typelist = typelist; TC_Conn.wlist = wlist; TC_Conn.delaylist = delaylist; 
TC_Conn.NumVL_of_M1 = NumVL_of_M1; 



end

