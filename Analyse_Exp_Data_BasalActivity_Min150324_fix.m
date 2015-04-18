%%%%%%%% Experimental Data Analysis
%%% Experiment data of base line activity (200 ms during light off)
%%
WTWT_MS = [26.86	28.03281345
17.2	26.53013345
3.665	12.27952048
5.99	15.57402385
5.62	22.21885032
1.16	8.547720921
19.37	27.93744195
11.8	25.00476264
28.4	32.55107443
5.225	14.83667385
22.995	28.85713199
14.775	21.87142359
4.78	16.33355283
12.175	23.5040088
6.32	19.61901407
14.965	22.83668024
18.18	24.4532686
21.63	32.89191899
7.87	18.27853558
6.92	15.97953435
4.13	19.17782766
9.17	18.7284528
4.415	18.07618927
5.74	15.39746847
7.095	19.35737104
10.425	29.78641098
3.435	13.16139916
9.645	19.07929447
8.73	19.21023761
6.93	16.2483261
13.2962963	25.50678635
11.58201058	22.93539692
32.76719577	39.36207525
11.93650794	23.02569219
4.835	14.73940064
0.955	8.514519195
];

KOKO_MS = [7.075	16.29041274
10.605	20.72167099
3.51	12.66885798
1.18	7.836789962
8.945	21.08654916
12.425	26.52371229
12.68	22.36247102
12.02	23.79055167
16.23	25.20846154
25.34	33.11225532
17.15	26.69768896
6.435	17.14723093
19.095	31.36405487
4.305	13.9172981
21.915	31.42699846
26.81	35.46521003
15.015	25.60530666
8.81	19.84018912
];
 
WTWT = [26.8600   17.2000    3.6650    5.9900    5.6200    1.1600   19.3700   11.8000   28.4000    5.2250  22.9950   14.7750    4.7800   12.1750    6.3200   14.9650   18.1800   21.6300    7.8700    6.9200   4.1300    9.1700    4.4150    5.7400    7.0950   10.4250    3.4350    9.6450    8.7300    6.9300   13.2963   11.5820   32.7672   11.9365    4.8350    0.9550];
KOKO = [ 7.0750   10.6050    3.5100    1.1800    8.9450   12.4250   12.6800   12.0200   16.2300   25.3400   17.1500    6.4350   19.0950    4.3050   21.9150   26.8100   15.0150    8.8100];
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
% f = 11.43*exp(0.06084*x);
Codename = 'PoisInputFr150324_fix';
 
WT_Expected_input = 0.1881.*WTWT.*WTWT + 1.5709.*WTWT; 
KO_Expected_input = 0.1881.*WTWT.*WTWT + 1.5709.*WTWT; 

% distribution of expected input
myl = 85; %myl = max(max(Input_FR_WT(:)), max(Input_FR_KO(:)));
fg_exp_in = figure; set(fg_exp_in, 'position',[ 420         319        1163         420]);
subplot(121); hist(WT_Expected_input,10); xlim([0 myl]);
ylabel('Number of cells'); xlabel('Mean firing rate (Hz)'); title('Distribution of expected input for WT');
subplot(122); hist(KO_Expected_input,10); xlim([0 myl]);
ylabel('Number of cells'); xlabel('Mean firing rate (Hz)'); title('Distribution of expected input for KO');
saveas(fg_exp_in,[Codename '_ExpInDist.fig'],'fig');
saveas(fg_exp_in,[Codename '_ExpInDist.jpg'],'jpg');
%% 
PoisInputFr_list_WT = WT_Expected_input;

PoisInputFr_list_KO = KO_Expected_input;

% y = datasample(data,k) returns k observations sampled uniformly at random, with replacement, from the data in data.
rng(1,'twister');

N_TRIAL = 5;  
N_CELL = 1150; 
Input_FR_WT = zeros(N_CELL, N_TRIAL);
Input_FR_KO = zeros(N_CELL, N_TRIAL);

for nn = 1 : N_TRIAL
    Input_FR_WT(:,nn) = datasample(PoisInputFr_list_WT,N_CELL);
    Input_FR_KO(:,nn) = datasample(PoisInputFr_list_KO,N_CELL);
    

TRIAL = nn;

dlmwrite([Codename '_WT_' num2str(TRIAL) '.txt'],[N_CELL; Input_FR_WT(:,nn)],'delimiter','\n')
dlmwrite([Codename '_KO_' num2str(TRIAL) '.txt'],[N_CELL; Input_FR_KO(:,nn)],'delimiter','\n')

myl = 85;
fg_sample = figure; set(fg_sample, 'position',[ 420         319        1163         420]);
subplot(121); hist(Input_FR_WT(:,TRIAL),15); xlim([0 myl]);
ylabel('Number of cells'); xlabel('Mean firing rate (Hz)'); title('Distribution of sample expected input for WT');
subplot(122); hist(Input_FR_KO(:,TRIAL),15); xlim([0 myl]);
ylabel('Number of cells'); xlabel('Mean firing rate (Hz)'); title('Distribution of sample expected input for KO');
suptitle(['Trial' num2str(TRIAL)])
saveas(fg_sample,[Codename '_' num2str(TRIAL) '.fig'],'fig')
saveas(fg_sample,[Codename '_' num2str(TRIAL) '.jpg'],'jpg')
end


for nn = 1 : N_TRIAL
    Sample_WT(:,nn) = datasample(WTWT,N_CELL);
    Sample_KO(:,nn) = datasample(KOKO,N_CELL);
    disp('-------------------------------------------------------------------')
    disp(['Trial # = ' num2str(nn)])
    disp(['WT: Sample mean = ' num2str(mean(Sample_WT(:,nn))) ', std = ' num2str(std(Sample_WT(:,nn)))]);
    disp(['KO: Sample mean = ' num2str(mean(Sample_KO(:,nn))) ', std = ' num2str(std(Sample_KO(:,nn)))]);
end
    