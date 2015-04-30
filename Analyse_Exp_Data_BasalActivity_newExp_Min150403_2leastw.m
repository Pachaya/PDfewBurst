%%%%%%%% Experimental Data Analysis
%%% Experiment data of base line activity (200 s during light off)
close all
WTWT_MS = [25.74821853	27.39890124
16.52706588	25.9163684
3.50543818	11.93101112
5.845730716	15.40425353
5.400675084	21.72911395
1.125140643	8.423298331
18.57732217	27.19027526
11.3064133	24.31725375
27.24840605	31.59219356
5.020627578	14.47925626
22.2127766	28.32763396
14.05675709	21.38098608
4.610576322	16.06268177
11.67145893	22.91561379
6.06575822	19.1599309
14.34179272	22.39785946
17.48218527	23.8527491
20.7575947	31.86267298
7.525940743	17.91221672
6.675834479	15.67482791
3.985498187	18.66560717
8.776097012	18.28861304
4.235529441	17.6278251
5.440680085	14.99529643
6.810851356	18.83897102
10.03625453	29.08088049
3.435429429	13.22282158
9.641205151	19.10890028
8.731091386	19.15905372
6.920865108	16.29030254
12.09151144	24.52680203
10.47130891	21.84684684
29.74871859	38.16396845
10.87635954	22.11881327
4.645580698	14.43167856
0.915114389	8.257100899
];

KOKO_MS = [6.755844481	15.90683901
10.16627078	20.22281133
3.340417552	12.33211359
1.135141893	7.571389928
8.501062633	20.41188695
11.92149019	25.84691186
12.17652207	21.86482018
11.53144143	23.12522809
15.49193649	24.50673476
24.41305163	32.46445587
16.42205276	26.14040405
6.165770721	16.70502284
18.35229404	30.58544215
4.170521315	13.69144952
21.0176272	30.56056366
25.8832354	34.68778861
14.4168021	24.95070003
8.521065133	19.42922795
];
 

WTWT =[25.74821853
16.52706588
3.50543818
5.845730716
5.400675084
1.125140643
18.57732217
11.3064133
27.24840605
5.020627578
22.2127766
14.05675709
4.610576322
11.67145893
6.06575822
14.34179272
17.48218527
20.7575947
7.525940743
6.675834479
3.985498187
8.776097012
4.235529441
5.440680085
6.810851356
10.03625453
3.435429429
9.641205151
8.731091386
6.920865108
12.09151144
10.47130891
29.74871859
10.87635954
4.645580698
0.915114389];

KOKO = [6.755844481
10.16627078
3.340417552
1.135141893
8.501062633
11.92149019
12.17652207
11.53144143
15.49193649
24.41305163
16.42205276
6.165770721
18.35229404
4.170521315
21.0176272
25.8832354
14.4168021
8.521065133];


%% Distribution of mean firing rate
fWT = figure; 
hist(WTWT, 20);
ylabel('Number of cells'); xlabel('Mean firing rate (Hz)'); title('WT');
xlim([0 35]);
fKO = figure; 
hist(KOKO,20); 
ylabel('Number of cells'); xlabel('Mean firing rate (Hz)'); title('KO');
xlim([0 35]);

%% 

sct_fg = figure; scatter(WTWT_MS(:,1),WTWT_MS(:,2),'b');
hold on; scatter(KOKO_MS(:,1),KOKO_MS(:,2),'r');
legend('WT','KO','Location','NorthWest');
title('Distribution of Mean and sigma of real experiment data')
xlabel('Mean average firing rate'); ylabel('sigma of mean firing rate'); 
% figure; hist(WTWT_MS(:,2) ./ WTWT_MS(:,1), 30 ); title('WT sig/mean'); xlabel('Average firing rate (Hz)'); ylabel('Number of cell');
% figure; hist(KOKO_MS(:,2) ./ KOKO_MS(:,1), 30); title('KO sig/mean'); xlabel('Average firing rate (Hz)'); ylabel('Number of cell');
% 

%% fit with polyval 
p = polyfit(WTWT_MS(:,1),WTWT_MS(:,2),3)
xx = 0:1:35;
yy = polyval(p,xx);
figure(sct_fg); hold on;
plot(xx,yy,'b')

p2 = polyfit(KOKO_MS(:,1),KOKO_MS(:,2),3) % 2 for power 2 , 3- for power3

xx = 0:1:35;
yy = polyval(p2,xx);
figure(sct_fg); hold on;
plot(xx,yy,'r')


%all
sct_fg = figure; scatter(WTWT_MS(:,1),WTWT_MS(:,2),'b');
hold on; scatter(KOKO_MS(:,1),KOKO_MS(:,2),'r');
legend('WT','KO','Location','NorthWest');
title('Distribution of Mean and sigma of real experiment data')
xlabel('Mean average firing rate'); ylabel('sigma of mean firing rate'); 
%fit all
all_MU = [KOKO_MS(:,1); WTWT_MS(:,1)];
all_SIG = [KOKO_MS(:,2); WTWT_MS(:,2)];
pp = polyfit(all_MU,all_SIG , 2)
% 0.8500   11.7696 
xx = 0:1:35;
yy = polyval(pp,xx);
figure(sct_fg); hold on;
plot(xx,yy,'k')


%% all data analysis
figure; 
hist(all_MU, 20);
xlabel('Average output firing rate')
ylabel('Number of Cells')
title({'Distribution of average firing rate measured in experiment',['min = ' num2str(min(all_MU)) ', med = ' num2str(median(all_MU)) ', max = ' num2str(max(all_MU))]})
%%%    21.0406
figure(sct_fg); scatter(median(all_MU),polyval(pp,median(all_MU)),80,'*k');

figure; 
hist(all_MU, 20);
xlabel('Average output firing rate')
ylabel('Number of Cells')
title('Distribution of average firing rate measured in experiment')
%%

p = polyfit(WTWT_MS(:,1),WTWT_MS(:,2),1)
% p =
% 0.7866   12.3619 % p(1) is the slope and p(2) is the intercept of the linear predictor.
a = p(1); b = p(2);
xx = 0:1:35;
yy = a.*xx+b;
figure(sct_fg); hold on;
plot(xx,yy,'b')

p2 = polyfit(KOKO_MS(:,1),KOKO_MS(:,2),1)
% p =
% 0.9928   10.1768 % p(1) is the slope and p(2) is the intercept of the linear predictor.
a2 = p2(1); b2 = p2(2);
xx = 0:1:35;
yy = a2.*xx+b2;
figure(sct_fg); hold on;
plot(xx,yy,'r')

%fit all
pp = polyfit([KOKO_MS(:,1); WTWT_MS(:,1)], [KOKO_MS(:,2); WTWT_MS(:,2)], 1)
% 0.8500   11.7696 
a3 = pp(1); b3 = pp(2);
xx = 0:1:35;
yy = a3.*xx+b3;
figure(sct_fg); hold on;
plot(xx,yy,'k')

figure;
scatter(WTWT_MS(:,1),WTWT_MS(:,2),'b');
hold on; scatter(KOKO_MS(:,1),KOKO_MS(:,2),'r');
plot(xx,yy,'c');
xlabel('Mean firing rate'); ylabel('sigma of firing rate)')


%% Expected input for simulation
X = WTWT;
WT_Expected_input = 64.78 * exp(0.0147*X) - 64.67; 

X= KOKO;
KO_Expected_input = 64.78 * exp(0.0147*X) - 64.67 ; 


 Codename = 'deltThisFold'; %'PoisInputFr150403_leastW_newExp';
 dirLoc = [Codename '\'];
 mkdir(dirLoc)
% distribution of expected input
myl = 50 ; %myl = max(max(Input_FR_WT(:)), max(Input_FR_KO(:)));
fg_exp_in = figure; set(fg_exp_in, 'position',[ 420         319        1163         420]);
subplot(121); hist(WT_Expected_input,10); xlim([0 myl]);
ylabel('Number of cells'); xlabel('Mean firing rate (Hz)'); title('Distribution of expected input for WT');
subplot(122); hist(KO_Expected_input,10); xlim([0 myl]);
ylabel('Number of cells'); xlabel('Mean firing rate (Hz)'); title('Distribution of expected input for KO');
saveas(fg_exp_in,[dirLoc Codename '_ExpInDist.fig'],'fig');
saveas(fg_exp_in,[dirLoc Codename '_ExpInDist.jpg'],'jpg');
%% 
PoisInputFr_list_WT = WT_Expected_input;

PoisInputFr_list_KO = KO_Expected_input;

% y = datasample(data,k) returns k observations sampled uniformly at random, with replacement, from the data in data.
rng(1,'twister');

N_TRIAL = 30;  
N_CELL = 1150; 
Input_FR_WT = zeros(N_CELL, N_TRIAL);
Input_FR_KO = zeros(N_CELL, N_TRIAL);

for nn = 1 : N_TRIAL
    Input_FR_WT(:,nn) = datasample(PoisInputFr_list_WT,N_CELL);
    Input_FR_KO(:,nn) = datasample(PoisInputFr_list_KO,N_CELL);
    

TRIAL = nn;

% dlmwrite([dirLoc Codename '_WT_' num2str(TRIAL) '.txt'],[N_CELL; Input_FR_WT(:,nn)],'delimiter','\n')
% dlmwrite([dirLoc Codename '_KO_' num2str(TRIAL) '.txt'],[N_CELL; Input_FR_KO(:,nn)],'delimiter','\n')
disp('WRITE')

myl = 50;
fg_sample = figure; set(fg_sample, 'position',[ 420         319        1163         420]);
subplot(121); hist(Input_FR_WT(:,TRIAL),15); xlim([0 myl]);
ylabel('Number of cells'); xlabel('Mean firing rate (Hz)'); title('Distribution of sample expected input for WT');
subplot(122); hist(Input_FR_KO(:,TRIAL),15); xlim([0 myl]);
ylabel('Number of cells'); xlabel('Mean firing rate (Hz)'); title('Distribution of sample expected input for KO');
suptitle(['Trial' num2str(TRIAL)])
saveas(fg_sample,[dirLoc Codename '_' num2str(TRIAL) '.fig'],'fig')
saveas(fg_sample,[dirLoc Codename '_' num2str(TRIAL) '.jpg'],'jpg')
end


for nn = 1 : N_TRIAL 
    Sample_WT(:,nn) = datasample(WTWT,N_CELL);
    Sample_KO(:,nn) = datasample(KOKO,N_CELL);
    disp('-------------------------------------------------------------------')
    disp(['Trial # = ' num2str(nn)])
    disp(['WT: Sample mean = ' num2str(mean(Sample_WT(:,nn))) ', std = ' num2str(std(Sample_WT(:,nn)))]);
    disp(['KO: Sample mean = ' num2str(mean(Sample_KO(:,nn))) ', std = ' num2str(std(Sample_KO(:,nn)))]);
end
    