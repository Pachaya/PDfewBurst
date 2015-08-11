% Fig 2
dirMat = 'AllTrials/';
load([ dirMat 'CnvrgntTypes_20Hz__rTC_50_100_wmTC_50_100_trial_1'])
% NewGennConn_r50_100_w_50_100_trial1 
% CnvrgntTypes_20Hz__rTC_50_100_wmTC_50_100_trial_1
PARAMETERS{1}.titleTxt = 'Range_c_o_n';
PARAMETERS{2}.titleTxt = 'W_c_o_n';
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
            figure(f1); subplot(nr,nc,cnt);  plot_paramMat( tmpMat, titleTxt ,xl,yl, xt,yt); caxis([0 35]);
            figure(f2); subplot(nr,nc,cnt);  plot_contourparamMat( tmpMat, titleTxt ,xl,yl, xt,yt, Ncontour); caxis([0 35]);
            
   saveMat{p5, ct_ii} = tmpMat;
   end
end

%%  Fix w = 80 -> p2 = 4
p2 = 4; 
%varies input
ff1 = figure;
RangeLabel = cell(length(PARAMETERS{1}.PARAM),1);
for ii = 1: length(PARAMETERS{1}.PARAM)
    RangeLabel{ii} = num2str(PARAMETERS{1}.PARAM(ii));
end
LEG = {'Static','Week Oscillation', 'Strong Oscillation'};
for ct_ii = 1 : NumCnvrgntTypes
    subplot(1,3,ct_ii);  C = 'kbr';
    for p5 = 1 : length(PARAMETERS{5}.PARAM)  %Synchronization Level
        tmpMat = saveMat{p5, ct_ii};
        plot(tmpMat(:,p2),'*-','Color',C(p5)); hold on;
    end      
    set(gca,'XTickLabel',RangeLabel);   
    xlabel('Range'); 
    ylabel('Average FR');
    title(CnvrgntTypes{ct_ii}.CodeName);
    legend(LEG,'location','best');
end
suptitle(['Fixed Weight = ' num2str(PARAMETERS{2}.PARAM(p2))]);

%varies con types
ff2 = figure;
SynLabel = {'Static','Week Oscillation', 'Strong Oscillation'};
LEG = {'G-G','S-U', 'S-E'};
for p5 = 1 : length(PARAMETERS{5}.PARAM)  %Synchronization Level  
    subplot(1,3,p5);  C = 'kbr';
    for ct_ii = 1 : NumCnvrgntTypes
        tmpMat = saveMat{p5, ct_ii};
        plot(tmpMat(:,p2),'*-','Color',C(ct_ii)); hold on;
    end      
    set(gca,'XTickLabel',RangeLabel);   
    xlabel('Range'); 
    ylabel('Average FR');
    title(SynLabel{p5});
    legend(LEG,'location','best');
end
suptitle(['Fixed Weight = ' num2str(PARAMETERS{2}.PARAM(p2))]);


%%  Fix r = 80 -> p1 = 4
p1 = 4; 
%varies input
ff1 = figure;
RangeLabel = cell(length(PARAMETERS{2}.PARAM),1);
for ii = 1: length(PARAMETERS{2}.PARAM)
    RangeLabel{ii} = num2str(PARAMETERS{2}.PARAM(ii));
end
LEG = {'Static','Week Oscillation', 'Strong Oscillation'};
for ct_ii = 1 : NumCnvrgntTypes
    subplot(1,3,ct_ii);  C = 'kbr';
    for p5 = 1 : length(PARAMETERS{5}.PARAM)  %Synchronization Level
        tmpMat = saveMat{p5, ct_ii};
        plot(tmpMat(p1,:),'*-','Color',C(p5)); hold on;
    end      
    set(gca,'XTickLabel',RangeLabel);   
    xlabel('Range'); 
    ylabel('Average FR');
    title(CnvrgntTypes{ct_ii}.CodeName);
    legend(LEG,'location','best');
end
suptitle(['Fixed Range = ' num2str(PARAMETERS{1}.PARAM(p1))]);

%varies con types
ff2 = figure;
SynLabel = {'Static','Week Oscillation', 'Strong Oscillation'};
LEG = {'G-G','S-U', 'S-E'};
for p5 = 1 : length(PARAMETERS{5}.PARAM)  %Synchronization Level  
    subplot(1,3,p5);  C = 'kbr';
    for ct_ii = 1 : NumCnvrgntTypes
        tmpMat = saveMat{p5, ct_ii};
        plot(tmpMat(p1,:),'*-','Color',C(ct_ii)); hold on;
    end      
    set(gca,'XTickLabel',RangeLabel);   
    xlabel('Range'); 
    ylabel('Average FR');
    title(SynLabel{p5});
    legend(LEG,'location','best');
end
suptitle(['Fixed Range = ' num2str(PARAMETERS{1}.PARAM(p1))]);

%%
% Average FR matrix  : L1
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
figure(f1); suptitle( 'Average FR matrix of L1')
figure(f2); suptitle( 'Average FR matrix of L1')
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
 p1 = 4; p2 = 3;   
   cnt2 = 0; L1_instFR = zeros(ii*jj,tstop);
   for ii = 1 : length(PARAMETERS{1}.PARAM)
       for jj = 2 : length(PARAMETERS{2}.PARAM)
           cnt2 = cnt2+1;
               L1_instFR(cnt2,:)  = (sum (ACT_Record{p1,p2,p3,p4,p5}.L2.spktrain)./tstop*1000);
       end
   end
   cnt = cnt +1; 
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