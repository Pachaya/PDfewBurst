% Make connection for layer 0 - layer 1 


layer0locTxt = 'Epos4600_noI_ssize3000_repmatVLpos';
layer1locTxt = 'Epos_E1150_I380_Regular_VL_Sweep5Itr5000_ceil_sorted';


[L0pos, L0Distmat, L0nnDist]= GetCellsDistInfo(layer0locTxt, 1);
[L1pos, L1Distmat, L1nnDist]= GetCellsDistInfo(layer1locTxt, 1); %VL