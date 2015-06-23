% Fig 2
dirMat = 'MatFiles/';
load([ dirMat 'CnvrgntTypes_20Hz__rTC_50_100_wmTC_50_100_trial_1_1_17-Jun-2015'])
%%
NumCnvrgntTypes = length(CnvrgntTypes);

% Average FR matrix  : L2
tstop =5000;
nr = length(PARAMETERS{5}.PARAM); nc=NumCnvrgntTypes;
saveMat = cell(nr,nc);
matsize = [length(PARAMETERS{1}.PARAM) length(PARAMETERS{2}.PARAM)];
cnt = 0; 
p3=1;p4=1;
f1 = figure; set(gcf,'position',[   286         101        1371         872]); 
f2 = figure; set(gcf,'position',[   286         101        1371         872]); 
Ncontour =3;
for ct_ii = 1 : NumCnvrgntTypes
    ACT_Record = CnvrgntTypes{ct_ii}.ACT_Record;
   for p5 = 1 : length(PARAMETERS{5}.PARAM)  %Synchronization Level
       
   tmpMat  = zeros(matsize);
   for ii = 1 : length(PARAMETERS{1}.PARAM)
       for jj = 2 : length(PARAMETERS{2}.PARAM)
           tmpMat(ii,jj) = mean(sum (ACT_Record{ii,jj,p3,p4,p5}.L2.spktrain,2)./tstop*1000);
       end
   end
   cnt = cnt +1; 
   
            xt = PARAMETERS{2}.PARAM; yt = PARAMETERS{1}.PARAM;
            xl = PARAMETERS{2}.titleTxt; yl = PARAMETERS{1}.titleTxt;
            titleTxt = ['Type : ' CnvrgntTypes{ct_ii}.CodeName ', Synch lvl = ' num2str(PARAMETERS{5}.PARAM(p5)) ];            
%    subplot(nr,nc,cnt); imagesc(tmpMat); axis xy; colormap('gray'); colorbar();
            figure(f1); subplot(nr,nc,cnt);  plot_paramMat( tmpMat, titleTxt ,xl,yl, xt,yt);
            figure(f2); subplot(nr,nc,cnt);  plot_contourparamMat( tmpMat, titleTxt ,xl,yl, xt,yt, Ncontour)            
   saveMat{p5, ct_ii} = tmpMat;
   end
end
%%
 Average FR matrix  : L1
tstop =5000;
nr = length(PARAMETERS{5}.PARAM); nc=NumCnvrgntTypes;
saveMat = cell(nr,nc);
matsize = [length(PARAMETERS{1}.PARAM) length(PARAMETERS{2}.PARAM)];
cnt = 0; 
p3=1;p4=1;
f1 = figure; set(gcf,'position',[   286         101        1371         872]); 
f2 = figure; set(gcf,'position',[   286         101        1371         872]); 
Ncontour =3;
for ct_ii = 1 : NumCnvrgntTypes
    ACT_Record = CnvrgntTypes{ct_ii}.ACT_Record;
   for p5 = 1 : length(PARAMETERS{5}.PARAM)  %Synchronization Level
       
   tmpMat  = zeros(matsize);
   for ii = 1 : length(PARAMETERS{1}.PARAM)
       for jj = 2 : length(PARAMETERS{2}.PARAM)
           tmpMat(ii,jj) = mean(sum (ACT_Record{ii,jj,p3,p4,p5}.L1.spktrain,2)./tstop*1000);
       end
   end
   cnt = cnt +1; 
   
            xt = PARAMETERS{2}.PARAM; yt = PARAMETERS{1}.PARAM;
            xl = PARAMETERS{2}.titleTxt; yl = PARAMETERS{1}.titleTxt;
            titleTxt = ['Type : ' CnvrgntTypes{ct_ii}.CodeName ', Synch lvl = ' num2str(PARAMETERS{5}.PARAM(p5)) ];            
%    subplot(nr,nc,cnt); imagesc(tmpMat); axis xy; colormap('gray'); colorbar();
            figure(f1); subplot(nr,nc,cnt);  plot_paramMat( tmpMat, titleTxt ,xl,yl, xt,yt);
            figure(f2); subplot(nr,nc,cnt);  plot_contourparamMat( tmpMat, titleTxt ,xl,yl, xt,yt, Ncontour)            
   saveMat{p5, ct_ii} = tmpMat;
   end
end
%% % Instant Average FR matrix  : L1
tstop =5000;
nr = length(PARAMETERS{5}.PARAM); nc=NumCnvrgntTypes;
saveMat = cell(nr,nc);
cnt = 0; 
p3=1;p4=1;
f1 = figure; set(gcf,'position',[   286         101        1371         872]); 
f2 = figure; set(gcf,'position',[   286         101        1371         872]); 
Ncontour =3;
for ct_ii = 1 : NumCnvrgntTypes
    ACT_Record = CnvrgntTypes{ct_ii}.ACT_Record;
   for p5 = 1 : length(PARAMETERS{5}.PARAM)  %Synchronization Level
       
   cnt2 = 0; L1_instFR = zeros(ii*jj,tstop);
   for ii = 1 : length(PARAMETERS{1}.PARAM)
       for jj = 2 : length(PARAMETERS{2}.PARAM)
           cnt2 = cnt2+1;
               L1_instFR(cnt2,:)  = (sum (ACT_Record{p1,p2,p3,p4,p5}.L2.spktrain)./tstop*1000);
       end
   end
    p1 = 4; p2 = 3;    cnt = cnt +1; 
    titleTxt = ['Type : ' CnvrgntTypes{ct_ii}.CodeName ', Synch lvl = ' num2str(PARAMETERS{5}.PARAM(p5)) ];  
%     L1_instFR  = (sum (ACT_Record{p1,p2,p3,p4,p5}.L1.spktrain)./tstop*1000);
    figure(f1); subplot(nr,nc,cnt); hold on; 
    for i = 1 : ii*jj
        plot(L1_instFR(i,:)); 
    end
        xlim([1000 1200]);
    title(titleTxt)

%    
%             xt = PARAMETERS{2}.PARAM; yt = PARAMETERS{1}.PARAM;
%             xl = PARAMETERS{2}.titleTxt; yl = PARAMETERS{1}.titleTxt;
%             titleTxt = ['Type : ' CnvrgntTypes{ct_ii}.CodeName ', Synch lvl = ' num2str(PARAMETERS{5}.PARAM(p5)) ];            
%    subplot(nr,nc,cnt); imagesc(tmpMat); axis xy; colormap('gray'); colorbar();
%             figure(f1); subplot(nr,nc,cnt);  plot_paramMat( tmpMat, titleTxt ,xl,yl, xt,yt);
%             figure(f2); subplot(nr,nc,cnt);  plot_contourparamMat( tmpMat, titleTxt ,xl,yl, xt,yt, Ncontour)            
%    saveMat{p5, ct_ii} = tmpMat;

   end
end

%%
figure; 
p1 = 1; p2=1; p3 = 1; p4 =1; p5 =1; % static
data1= sum(CnvrgntTypes{1}.ACT_Record{p1,p2,p3,p4,p5}.L2.spktrain)./tstop*1000;   xlim([1000 1200]);
data2= sum(CnvrgntTypes{2}.ACT_Record{p1,p2,p3,p4,p5}.L2.spktrain)./tstop*1000;   xlim([1000 1200]);
data3= sum(CnvrgntTypes{3}.ACT_Record{p1,p2,p3,p4,p5}.L2.spktrain)./tstop*1000;   xlim([1000 1200]);
subplot(131);  imagesc(data1); 
subplot(132);  imagesc(data2); 
subplot(133);  imagesc(data3); 

%%
figure; 
p1 = 1; p2=1; p3 = 1; p4 =1; p5 =1; % static
data1= CnvrgntTypes{1}.ACT_Record{p1,p2,p3,p4,p5}.L1.spktrain;
data2= CnvrgntTypes{2}.ACT_Record{p1,p2,p3,p4,p5}.L1.spktrain;
data3= CnvrgntTypes{3}.ACT_Record{p1,p2,p3,p4,p5}.L1.spktrain;
subplot(131);  imagesc(data1); 
subplot(132);  imagesc(data2); 
subplot(133);  imagesc(data3); 

%% figure; 
p1 = 1; p2=1; p3 = 1; p4 =1; p5 =1; % static
data1= CnvrgntTypes{1}.ACT_Record{p1,p2,p3,p4,p5}.L1.spktrain; 
data2= CnvrgntTypes{2}.ACT_Record{p1,p2,p3,p4,p5}.L1.spktrain;
data3= CnvrgntTypes{3}.ACT_Record{p1,p2,p3,p4,p5}.L1.spktrain;
subplot(131);  imagesc(data1); 
subplot(132);  imagesc(data2); 
subplot(133);  imagesc(data3); 