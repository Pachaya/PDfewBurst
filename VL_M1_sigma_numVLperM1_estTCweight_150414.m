% finding the relation between sigma and number of VL cells per M1 in
% thalamocortical connection
%% Parameters Setting 
 RangeE = 200 ; 
 SigmaTC = RangeE/sqrt(2);
%  Sim_Postfix = 

%% Download cell positions
dirLoc = '../';
% VL 
VL_E_fname = 'Epos_E1150_I380_Regular_VL_Sweep5Itr5000_ceil_sorted.txt';
VL_I_fname = 'Ipos_E1150_I380_Regular_VL_Sweep5Itr5000_ceil_sorted.txt';
% VL_Epos = importdata([dirLoc 'Neurons_location_' Name_postfix '.txt']);
tmpData = importdata([dirLoc VL_E_fname]);
VL_nE=tmpData(1);  
VL_Epos = zeros(VL_nE,2);
VL_Epos(:,1) = tmpData(2:2:end);  
VL_Epos(:,2) = tmpData(3:2:end);  

tmpData = importdata([dirLoc VL_I_fname ]);
VL_nI=tmpData(1);  
VL_Ipos = zeros(VL_nI,2);
VL_Ipos(:,1) = tmpData(2:2:end);  
VL_Ipos(:,2) = tmpData(3:2:end);  
% M1
M1_E_fname = 'Epos_E5625_I1892_ssize1500_Regular_sweep5Itr5000_ceil_sorted.txt';
M1_I_fname = 'Ipos_E5625_I1892_ssize1500_Regular_sweep5Itr5000_ceil_sorted.txt';

tmpData = importdata([dirLoc M1_E_fname ]); 
M1_nE=tmpData(1);  
M1_Epos = zeros(M1_nE,2);
M1_Epos(:,1) = tmpData(2:2:end);  
M1_Epos(:,2) = tmpData(3:2:end);  

tmpData = importdata([dirLoc M1_I_fname ]); 
M1_nI=tmpData(1);  
M1_Ipos = zeros(M1_nI,2);
M1_Ipos(:,1) = tmpData(2:2:end);  
M1_Ipos(:,2) = tmpData(3:2:end);  


%% Plot  the M1's E and VL's E together
figure; 
scatter(M1_Epos(:,1),M1_Epos(:,2),20, 'ok'); hold on;
scatter(M1_Epos(:,1),M1_Epos(:,2), '.k'); hold on;
scatter(VL_Epos(:,1),VL_Epos(:,2),20, 'or'); hold on;
scatter(VL_Epos(:,1),VL_Epos(:,2), '.r'); hold on;
axis equal xy

%% Finding distance between M1 and VL when projecting on the same plane 
D = pdist2(M1_Epos, VL_Epos);  % D(1,:) = distance from M1 cell#1 to all VL cell,  (D(:,1) = distance from M1 cell#1 to all VL cell
% Nearest neighbor 
nearestM1ofVL = sort(D,1);  nearestM1ofVL  =  nearestM1ofVL(1,:);
nearestVLofM1 = sort(D,2); nearestVLofM1 = nearestVLofM1(:,1);
figure; 
subplot(211)
hist(nearestM1ofVL,50); title(['Distribution of nearest M1 of each VL cell (E cell only), average = ' num2str(mean(nearestM1ofVL))]);
subplot(212)
hist(nearestVLofM1,100); title(['Distribution of nearest VL of each M1 cell (E cell only), average = ' num2str(mean(nearestVLofM1))]);

%% Get Center Cell ID 
RangeTC_List = 10 :10:320; %max about 300
avgVL_M1 = zeros(length(RangeTC_List),1);
for isig = 1 : length(RangeTC_List)
    
Sigma_Val = RangeTC_List(isig)/sqrt(2); %200/sqrt(2); 
disp([ 'TC range = ' num2str(RangeTC_List(isig)) ' after sqrt2 =' num2str(Sigma_Val)])
bordersize = ceil(Sigma_Val*3.5);
ssize = 1500; LowerBorder = bordersize ; UpperBorder = ssize-bordersize ;
THRESHOLD = -55; RES = 1; %temporal resolution of one bin
CUTTIME = 500;
% VL to     M1  : 60 ----> 210
% VL to     M1  : 50 ----> 190
% VL to     M1  : 40 ----> 170
if(LowerBorder > ssize/2)
    disp(['Error:: The sigma is too big for this neural patch'] )
end
getCenterPart = 1;

if(getCenterPart)
%get Center Part
centerID_VL = find((VL_Epos(:,1) > LowerBorder) & (VL_Epos(:,1) <UpperBorder)&(VL_Epos(:,2) > LowerBorder) & (VL_Epos(:,2) <UpperBorder));
centerID_M1 = find((M1_Epos(:,1) > LowerBorder) & (M1_Epos(:,1) <UpperBorder)&(M1_Epos(:,2) > LowerBorder) & (M1_Epos(:,2) <UpperBorder));
NcenterE_VL = length(centerID_VL); NcenterE_M1 = length(centerID_M1); 
% somaVall = somaVall ([centerID_E; nE+centerID_I],:);
% somaVallE = somaVallE(centerID_E,:);
% somaVallI = somaVallI(centerID_I,:);
% N_E = nE; N_I = nI;
% nE = size(somaVallE,1); nI = size(somaVallI,1);
end

%%  Counting the number of VL cell that connect to a M1 cell (only cells in M1 center)
Pconn = 0.85; 
NumVLperM1 = 1:length(centerID_M1);
VLtoM1 = cell(length(centerID_M1),1);

for m_id = 1 : length(centerID_M1)
    dist = D(centerID_M1(m_id),:);
    yy=Pconn*func_gauss(dist, Sigma_Val);
    ppp = rand(1,VL_nE) ; 
    VLtoM1{m_id} = find(ppp<yy);
    NumVLperM1(m_id) = length(VLtoM1{m_id});
end    
disp(['Average number of VL connection to M1 is approximately = ' num2str(mean(NumVLperM1))]);
avgVL_M1(isig) = mean(NumVLperM1);
% m_id = 1;
% figure; 
% scatter(M1_Epos(:,1),M1_Epos(:,2),20, 'ok'); hold on;
% scatter(M1_Epos(centerID_M1(m_id),1),M1_Epos(centerID_M1(m_id),2), '.k'); hold on;
% scatter(VL_Epos(:,1),VL_Epos(:,2),20, 'or'); hold on;
% scatter(VL_Epos( VLtoM1{m_id},1),VL_Epos( VLtoM1{m_id},2), '.r'); hold on;
% axis equal xy

end
%% 

% TC range = 10 after sqrt2 =7.0711
% Average number of VL connection to M1 is approximately = 0.13368
% TC range = 20 after sqrt2 =14.1421
% Average number of VL connection to M1 is approximately = 0.55583
% TC range = 30 after sqrt2 =21.2132
% Average number of VL connection to M1 is approximately = 1.2062
% TC range = 40 after sqrt2 =28.2843
% Average number of VL connection to M1 is approximately = 2.1674
% TC range = 50 after sqrt2 =35.3553
% Average number of VL connection to M1 is approximately = 3.389
% TC range = 60 after sqrt2 =42.4264
% Average number of VL connection to M1 is approximately = 4.9442
% TC range = 70 after sqrt2 =49.4975
% Average number of VL connection to M1 is approximately = 6.6889
% TC range = 80 after sqrt2 =56.5685
% Average number of VL connection to M1 is approximately = 8.7506
% TC range = 90 after sqrt2 =63.6396
% Average number of VL connection to M1 is approximately = 11.0602
% TC range = 100 after sqrt2 =70.7107
% Average number of VL connection to M1 is approximately = 13.6679
% TC range = 110 after sqrt2 =77.7817
% Average number of VL connection to M1 is approximately = 16.5564
% TC range = 120 after sqrt2 =84.8528
% Average number of VL connection to M1 is approximately = 19.6114
% TC range = 130 after sqrt2 =91.9239
% Average number of VL connection to M1 is approximately = 23.0566
% TC range = 140 after sqrt2 =98.9949
% Average number of VL connection to M1 is approximately = 26.6759
% TC range = 150 after sqrt2 =106.066
% Average number of VL connection to M1 is approximately = 30.7037
% TC range = 160 after sqrt2 =113.1371
% Average number of VL connection to M1 is approximately = 35.0288
% TC range = 170 after sqrt2 =120.2082
% Average number of VL connection to M1 is approximately = 39.5193
% TC range = 180 after sqrt2 =127.2792
% Average number of VL connection to M1 is approximately = 44.307
% TC range = 190 after sqrt2 =134.3503
% Average number of VL connection to M1 is approximately = 49.1454
% TC range = 200 after sqrt2 =141.4214
% Average number of VL connection to M1 is approximately = 54.7388
% TC range = 210 after sqrt2 =148.4924
% Average number of VL connection to M1 is approximately = 60.0914
% TC range = 220 after sqrt2 =155.5635
% Average number of VL connection to M1 is approximately = 65.8969
% TC range = 230 after sqrt2 =162.6346
% Average number of VL connection to M1 is approximately = 72.4656
% TC range = 240 after sqrt2 =169.7056
% Average number of VL connection to M1 is approximately = 79.0992
% TC range = 250 after sqrt2 =176.7767
% Average number of VL connection to M1 is approximately = 85.7229
% TC range = 260 after sqrt2 =183.8478
% Average number of VL connection to M1 is approximately = 92.9123
% TC range = 270 after sqrt2 =190.9188
% Average number of VL connection to M1 is approximately = 99.4154
% TC range = 280 after sqrt2 =197.9899
% Average number of VL connection to M1 is approximately = 107.7742
% TC range = 290 after sqrt2 =205.061
% Average number of VL connection to M1 is approximately = 112.375
% TC range = 300 after sqrt2 =212.132
% Average number of VL connection to M1 is approximately = 121
% TC range = 310 after sqrt2 =219.2031
% Error:: The sigma is too big for this neural patch
% Average number of VL connection to M1 is approximately = NaN
% TC range = 320 after sqrt2 =226.2742
% Error:: The sigma is too big for this neural patch
% Average number of VL connection to M1 is approximately = NaN

% rTC maximum at 300


%% Plot 
figure;
plot(RangeTC_List,avgVL_M1,'o')
xlabel('Thalamocortical connection range(um)')
ylabel('average #VL/M1')
title('TC range VS. average # of VL connection per M1')

%% get center ID VL
figure; 
scatter(VL_Epos(:,1), VL_Epos(:,2),20,'r'); axis square
hold on;
   scatter(VL_Epos(centerID_VL,1), VL_Epos(centerID_VL,2),'.r'); axis square
for ii = 1:length(centerID_VL)
    scatter(VL_Epos(centerID_VL(ii),1), VL_Epos(centerID_VL(ii),2),'.k'); axis square
    title(['#' num2str(ii) ': ' num2str(centerID_VL(ii))])
%     k = waitforbuttonpress;
end

%%% Center ID = centerID_VL(19) => 576
%% get center ID M1
figure; 
scatter(M1_Epos(:,1), M1_Epos(:,2),20,'r'); axis square
hold on;
scatter(M1_Epos(centerID_M1,1), M1_Epos(centerID_M1,2),'.r'); axis square
for ii = 1:length(centerID_M1)
    scatter(M1_Epos(centerID_M1(ii),1), M1_Epos(centerID_M1(ii),2),'.k'); axis square
    title(['#' num2str(ii) ': ' num2str(centerID_M1(ii))])
%     k = waitforbuttonpress;
end
%%% Center ID = centerID_M1(87) => 2816
% scatter(VL_Epos(:,1), VL_Epos(:,2),50,'b'); axis square
% scatter(VL_Epos(centerID_VL(19),1), VL_Epos(centerID_VL(19),2),'.g'); axis square


%% Test  AP occur when w = 0.001454
% Area under Gaussian = N???sqrt(2?)  --> connectivity and weighting factor
% Volume under 2D Gaussian = peak * sigma^2 *2*pi
% 0.001454 = (pmax * sigma^2 *2*pi )*( wmax*W_scale * sigma^2 *2*pi)
% 
% wmax = 0.001454 /((pmax * sigma^2 *2*pi )*( W_scale * sigma^2 *2*pi))
%

ssize = 1500; dist = ssize/2; pmax =0.85; W_scale = 1E-05;
sigmaLST = RangeTC_List./sqrt(2);
est_wMax = zeros(size(sigmaLST));
disp('rangeTC   sigmaTC   #VL/M1   Weight    Area')
for ii = 1:length(sigmaLST)
    sigma = sigmaLST(ii);
    est_wMax(ii) = 0.001454 /((pmax * sigma^2 *2*pi )*( W_scale * sigma^2 *2*pi));
    AA = (pmax * sigma^2 *2*pi )*( est_wMax(ii) *W_scale * sigma^2 *2*pi);
    disp([num2str(RangeTC_List(ii)) '   ' num2str(sigma) '   ' num2str(avgVL_M1(ii)) '   ' num2str(est_wMax(ii)) '   ' num2str(AA)]);
end
% 
% rangeTC   sigmaTC   #VL/M1   Weight    Area
% 10   7.0711   0.13711   9.651   0.001454
% 20   14.1421   0.55031   4.8255   0.001454
% 30   21.2132   1.2221   3.217   0.001454
% 40   28.2843   2.1802   2.4127   0.001454
% 50   35.3553   3.4284   1.9302   0.001454
% 60   42.4264   4.9015   1.6085   0.001454
% 70   49.4975   6.7161   1.3787   0.001454
% 80   56.5685   8.8197   1.2064   0.001454
% 90   63.6396   11.1055   1.0723   0.001454
% 100   70.7107   13.6325   0.9651   0.001454
% 110   77.7817   16.5113   0.87736   0.001454
% 120   84.8528   19.5508   0.80425   0.001454
% 130   91.9239   22.9488   0.74238   0.001454
% 140   98.9949   26.6988   0.68935   0.001454
% 150   106.066   30.8536   0.6434   0.001454
% 160   113.1371   34.8985   0.60319   0.001454
% 170   120.2082   39.3853   0.5677   0.001454
% 180   127.2792   44.1146   0.53616   0.001454
% 190   134.3503   48.9601   0.50795   0.001454
% 200   141.4214   54.3524   0.48255   0.001454
% 210   148.4924   60.0095   0.45957   0.001454
% 220   155.5635   66.0576   0.43868   0.001454
% 230   162.6346   71.6469   0.41961   0.001454
% 240   169.7056   78.2975   0.40212   0.001454
% 250   176.7767   84.759   0.38604   0.001454
