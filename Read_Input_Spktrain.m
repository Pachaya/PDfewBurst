% load_input_spktrain(SimCode)
% SimCode = 'SE_InputFR20_OscF40_OscrltAmp1_Wscale1e-05_W50_range200_Trial1_T5000'; 

fname = [dirLoc 'InputSpkTrain_' SimCode '.txt'];

fid = fopen(fname);

tline = fgets(fid);
% first Line is nE nI tstop , can be ignore 
nL1 = 1150; 
tstop = 5000;

SPKTRAIN = zeros(nL1,tstop);

tline = fgets(fid);
while ischar(tline)
    aaa = strsplit(tline,'\t'); % first number is ID, 2nd - end = spike timing
    L1ID = str2num(aaa{1}) + 1; % NEURON numbering start from 0
    for ii = 2 : length(aaa)
        if(str2num(aaa{ii}) ~= 0)
         SPKTRAIN(L1ID, str2num(aaa{ii})) = 1;
        end
    end
         
%     disp(tline)
    tline = fgets(fid);
end

fclose(fid);

SPKTRAIN = sparse(SPKTRAIN); 

% figure; imagesc(SPKTRAIN)
% mean(sum(SPKTRAIN,2)) / 5