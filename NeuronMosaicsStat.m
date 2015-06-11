


% 
VLrawTxt = 'Epos_E1150_I380_Regular_VL_Sweep5Itr5000_ceil_sorted';
centerM1rawTxt =  'ECenterpos_rTC250_E5625_I1892_ssize1500_Regular_sweep5Itr5000_ceil_sorted';
fullM1rawTxt = 'Epos_E5625_I1892_ssize1500_Regular_sweep5Itr5000_ceil_sorted';


% coor = VLpos(:,1) + 1i*VLpos(:,2);
% dd = abs(coor);

%VL
[VLCellpos, VLDistmat, VLnnDist]= GetCellsDistInfo(VLrawTxt, 1);
% M1 full & center
[Cellpos, Distmat, nnDist]= GetCellsDistInfo(centerM1rawTxt, 1);
[Cellpos, Distmat, nnDist]= GetCellsDistInfo(fullM1rawTxt, 1);

ExpandCellpos = Cellpos *2;

Distmat = squareform (pdist(ExpandCellpos));
tmpnn = sort(Distmat,2);
nnDist = tmpnn(:,2);

figure; hist(nnDist,25); 
xlabel('Nearest Neighbor Distance (um)');
ylabel('# of Cells'); 
title({'Distribution of nearest neighbor distance', ['Mean = ' num2str(mean(nnDist)) ', std = ' num2str(std(nnDist))]});



%% get VL pos and stitch it for 4 times to make 3000 by 3000


VLpos_j_copy = VLpos +  repmat([0 1500],[1150 1]);
VLpos_i_copy = VLpos + repmat([1500 0],[1150 1]);
VLpos_ij_copy = VLpos +  repmat([1500 1500],[1150 1]);

repVLpos = [VLpos; VLpos_i_copy; VLpos_j_copy; VLpos_ij_copy; ];


Cellpos = repVLpos;
figure; scatter(Cellpos(:,1),Cellpos(:,2),20); axis square;
xlabel('um');  ylabel('um'); title('Distribution of cell position');


Distmat = squareform (pdist(Cellpos));
tmpnn = sort(Distmat,2);
nnDist = tmpnn(:,2);


figure; hist(nnDist,25); 
xlabel('Nearest Neighbor Distance (um)');
ylabel('# of Cells'); 
title({'Distribution of nearest neighbor distance', ['Mean = ' num2str(mean(nnDist)) ', std = ' num2str(std(nnDist))]});

%sort cell position
layer0_Pos = sortrows(repVLpos,1);