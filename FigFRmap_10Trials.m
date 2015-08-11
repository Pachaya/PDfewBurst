% Fig 2
dirMat = 'AllTrials/';
load([ dirMat 'FRmatrix_L2_rowSynchLvl_colCnvrgntTypes_Trial1_10_30-Jun-2015']);
 CnvrgntTypesINFO{2}.CodeName = 'UU';  CnvrgntTypesINFO{3}.CodeName = 'UE'; 
 
%%  Get Number of Cells 
load([dirMat 'NewGennConn_r50_100_w_50_100_trial1']);
%% num Cell
numMat_G = zeros(6,6); numMat_U = zeros(6,6); numMat_E = zeros(6,6);
for rr = 1 :6 
    for ww = 1: 6
        numMat_G(rr,ww) = mean(sum(saveConn_g{rr,ww}.Conn_Mat));
        numMat_U(rr,ww) = mean(sum(saveConn_u{rr,ww}.Conn_Mat));
        numMat_E(rr,ww) = mean(sum(saveConn_e{rr,ww}.Conn_Mat));
    end
end

figure; 
subplot(131); imagesc( numMat_G); title('GG'); c=colorbar; ylabel(c,'# cells');
subplot(132); imagesc( numMat_U); title('UU'); c=colorbar; ylabel(c,'# cells');
subplot(133); imagesc( numMat_E); title('UE'); c=colorbar; ylabel(c,'# cells');
m1 = mean(numMat_G,2); m2 = mean(numMat_U,2); m3 = mean(numMat_E,2);
numCell_LST = mean([m1 m2 m3],2);
%% num W per num
numMat_G = zeros(6,6); numMat_U = zeros(6,6); numMat_E = zeros(6,6);
for rr = 1 :6 
    for ww = 1: 6
        numMat_G(rr,ww) = mean(sum(saveConn_g{rr,ww}.Conn_Mat));
        numMat_U(rr,ww) = mean(sum(saveConn_u{rr,ww}.Conn_Mat));
        numMat_E(rr,ww) = mean(sum(saveConn_e{rr,ww}.Conn_Mat));
    end
end

figure; 
subplot(131); imagesc( numMat_G); title('GG'); c=colorbar; ylabel(c,'# cells');
subplot(132); imagesc( numMat_U); title('UU'); c=colorbar; ylabel(c,'# cells');
subplot(133); imagesc( numMat_E); title('UE'); c=colorbar; ylabel(c,'# cells');
m1 = mean(numMat_G,2); m2 = mean(numMat_U,2); m3 = mean(numMat_E,2);
numCell_LST = mean([m1 m2 m3],2);
%%
NumCnvrgntTypes = length(CnvrgntTypesINFO);
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
Synchlvl_LST = {'Static','Weak Oscillation', 'Strong Oscillation'};
for ct_ii = 1 : NumCnvrgntTypes
   for p5 = 1 : length(PARAMETERS{5}.PARAM)  %Synchronization Level
       
   tmpMat  = mean(FRmat_allTrials{p5,ct_ii}, 3);
   cnt = cnt +1;    
            xt = PARAMETERS{2}.PARAM; yt = PARAMETERS{1}.PARAM;
%             xl = PARAMETERS{2}.titleTxt; yl = PARAMETERS{1}.titleTxt;
xl = 'W [a.u.]'; yl = 'Range [um]';
            titleTxt = ['Type : ' CnvrgntTypesINFO{ct_ii}.CodeName ', ' Synchlvl_LST{p5} ];            
%    subplot(nr,nc,cnt); imagesc(tmpMat); axis xy; colormap('gray'); colorbar();
            figure(f1); subplot(nr,nc,cnt);  plot_paramMat( tmpMat, titleTxt ,xl,yl, xt,yt); caxis([0 35]); 
            c=colorbar; ylabel(c,'Fr(Hz)')
            figure(f2); subplot(nr,nc,cnt);  plot_contourparamMat( tmpMat, titleTxt ,xl,yl, xt,yt, Ncontour); caxis([0 35]);
              c=colorbar; ylabel(c,'Fr(Hz)')
            
   saveMat{p5, ct_ii} = tmpMat;
   end
end

%%  Fix w = 80 -> p2 = 4

p2 = 4; 
%varies input
ff1 = figure; set(ff1, 'position',[ 407         381        1382         420]);
RangeLabel = cell(length(PARAMETERS{1}.PARAM),1);
for ii = 1: length(PARAMETERS{1}.PARAM)
    RangeLabel{ii} = num2str(PARAMETERS{1}.PARAM(ii));
end
xlen = length(PARAMETERS{1}.PARAM);
LEG = {'S','WO', 'SO'};
for ct_ii = 1 : NumCnvrgntTypes
    subplot(1,3,ct_ii);  C = 'kbr';
    LEG2 = cell(length(PARAMETERS{5}.PARAM),1);
    for p5 = 1 : length(PARAMETERS{5}.PARAM)  %Synchronization Level
        tmpMat = saveMat{p5, ct_ii};
        plot(tmpMat(:,p2),'*-','Color',C(p5)); hold on;
        
        x = 1: length(PARAMETERS{1}.PARAM); y = tmpMat(:,p2)';
        p = polyfit(x,y,1)'; 
        yfit = polyval(p,x);
        LEG2{p5} = [LEG{p5}]; % ' y = ' num2str(p(1)) 'x + ' num2str(p(2))];
    end      
    set(gca,'XTick',1:length(PARAMETERS{1}.PARAM));  
    set(gca,'XTickLabel',RangeLabel);     set(gca,'FontSize',14);
    xlabel('Range'); 
    ylabel('Average FR');
    title(CnvrgntTypesINFO{ct_ii}.CodeName);
    legend(LEG2,'location','best');
    box off
end
suptitle(['Fixed Weight = ' num2str(PARAMETERS{2}.PARAM(p2))]);
%%
%varies con types
ff2 = figure; set(ff2, 'position',[ 407         381        1382         420]);
SynLabel = {'Static','Week Oscillation', 'Strong Oscillation'};
LEG = {'GG','UU', 'UE'};
for p5 = 1 : length(PARAMETERS{5}.PARAM)  %Synchronization Level  
    subplot(1,3,p5);  C = 'kbr';
      LEG2 = cell(length(PARAMETERS{5}.PARAM),1);
    for ct_ii = 1 : NumCnvrgntTypes
        tmpMat = saveMat{p5, ct_ii};
        plot(tmpMat(:,p2),'*-','Color',C(ct_ii)); hold on;
        
        x = 1: length(PARAMETERS{1}.PARAM); y = tmpMat(:,p2)';
        p = polyfit(x,y,1)'; 
        yfit = polyval(p,x);
        LEG2{ct_ii} = [LEG{ct_ii} ' y = ' num2str(p(1)) 'x + ' num2str(p(2))];
    end      
    set(gca,'XTickLabel',RangeLabel);     set(gca,'FontSize',14);
    xlabel('Range'); 
    ylabel('Average FR');
    title(SynLabel{p5});
    box off
    legend(LEG2,'location','best');
end
suptitle(['Fixed Weight = ' num2str(PARAMETERS{2}.PARAM(p2))]);


%%  Fix r = 80 -> p1 = 4
p1 = 4; 
%varies input
ff3 = figure; set(ff3, 'position',[ 220         247        1435         554]);
RangeLabel = cell(length(PARAMETERS{2}.PARAM),1);
for ii = 1: length(PARAMETERS{2}.PARAM)
    RangeLabel{ii} = num2str(PARAMETERS{2}.PARAM(ii));
end

LEG = {'Static','Weak oscilation', 'Strong oscilation'};
for ct_ii = 1 : NumCnvrgntTypes

    subplot(1,3,ct_ii);  C = 'kbr';
    LEG2 = cell(length(PARAMETERS{5}.PARAM),1);
    for p5 = 1 : length(PARAMETERS{5}.PARAM)  %Synchronization Level
        tmpMat = saveMat{p5, ct_ii};
        plot(tmpMat(p1,:),'*-','Color',C(p5)); hold on;
        
        x = 1: length(PARAMETERS{1}.PARAM); y = tmpMat(p1,:);
        p = polyfit(x,y,1)'; 
        yfit = polyval(p,x);
        LEG2{p5} = [LEG{p5}];% ' y = ' num2str(p(1)) 'x + ' num2str(p(2))];
    end      
    set(gca,'XTick',1:length(PARAMETERS{1}.PARAM));  
    set(gca,'XTickLabel',RangeLabel);     set(gca,'FontSize',14);
    xlabel('W[a.u]'); 
    ylabel('Average FR');
    title(CnvrgntTypesINFO{ct_ii}.CodeName);
    legend(LEG2,'location','best'); 
    box off
    
end
suptitle(['Fixed Range = ' num2str(PARAMETERS{1}.PARAM(p1))]);
%%
%varies con types

ff4 = figure; set(ff4, 'position',[ 407         381        1382         420]);
SynLabel = {'Static','Week Oscillation', 'Strong Oscillation'};
LEG = {'GG','UU', 'UE'};
for p5 = 1 : length(PARAMETERS{5}.PARAM)  %Synchronization Level  
    subplot(1,3,p5);  C = 'kbr';
     LEG2 = cell(length(PARAMETERS{5}.PARAM),1);
    for ct_ii = 1 : NumCnvrgntTypes
        tmpMat = saveMat{p5, ct_ii};
        plot(tmpMat(p1,:),'*-','Color',C(ct_ii)); hold on;
        
        x = 1: length(PARAMETERS{1}.PARAM); y = tmpMat(p1,:);
        p = polyfit(x,y,1)'; 
        yfit = polyval(p,x);
        LEG2{ct_ii} = [LEG{ct_ii} ' y = ' num2str(p(1)) 'x + ' num2str(p(2))];
    end      
    set(gca,'XTickLabel',RangeLabel);    set(gca,'FontSize',14);
    xlabel('W'); 
    ylabel('Average FR');
    title(SynLabel{p5});
    legend(LEG2,'location','best');
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
    ACT_Record = CnvrgntTypesINFO{ct_ii}.ACT_Record;
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
            titleTxt = ['Type : ' CnvrgntTypesINFO{ct_ii}.CodeName ', Synch lvl = ' num2str(PARAMETERS{5}.PARAM(p5)) ];            
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
    ACT_Record = CnvrgntTypesINFO{ct_ii}.ACT_Record;
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
    titleTxt = ['Type : ' CnvrgntTypesINFO{ct_ii}.CodeName ', Synch lvl = ' num2str(PARAMETERS{5}.PARAM(p5)) ];  
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
%             titleTxt = ['Type : ' CnvrgntTypesINFO{ct_ii}.CodeName ', Synch lvl = ' num2str(PARAMETERS{5}.PARAM(p5)) ];            
%    subplot(nr,nc,cnt); imagesc(tmpMat); axis xy; colormap('gray'); colorbar();
%             figure(f1); subplot(nr,nc,cnt);  plot_paramMat( tmpMat, titleTxt ,xl,yl, xt,yt);
%             figure(f2); subplot(nr,nc,cnt);  plot_contourparamMat( tmpMat, titleTxt ,xl,yl, xt,yt, Ncontour)            
%    saveMat{p5, ct_ii} = tmpMat;

   end
end

%%
figure; 
p1 = 1; p2=1; p3 = 1; p4 =1; p5 =1; % static
data1= sum(CnvrgntTypesINFO{1}.ACT_Record{p1,p2,p3,p4,p5}.L2.spktrain)./tstop*1000;   xlim([1000 1200]);
data2= sum(CnvrgntTypesINFO{2}.ACT_Record{p1,p2,p3,p4,p5}.L2.spktrain)./tstop*1000;   xlim([1000 1200]);
data3= sum(CnvrgntTypesINFO{3}.ACT_Record{p1,p2,p3,p4,p5}.L2.spktrain)./tstop*1000;   xlim([1000 1200]);
subplot(131);  imagesc(data1); 
subplot(132);  imagesc(data2); 
subplot(133);  imagesc(data3); 

%%
figure; 
p1 = 1; p2=1; p3 = 1; p4 =1; p5 =1; % static
data1= CnvrgntTypesINFO{1}.ACT_Record{p1,p2,p3,p4,p5}.L1.spktrain;
data2= CnvrgntTypesINFO{2}.ACT_Record{p1,p2,p3,p4,p5}.L1.spktrain;
data3= CnvrgntTypesINFO{3}.ACT_Record{p1,p2,p3,p4,p5}.L1.spktrain;
subplot(131);  imagesc(data1); 
subplot(132);  imagesc(data2); 
subplot(133);  imagesc(data3); 

%% figure; 
p1 = 1; p2=1; p3 = 1; p4 =1; p5 =1; % static
data1= CnvrgntTypesINFO{1}.ACT_Record{p1,p2,p3,p4,p5}.L1.spktrain; 
data2= CnvrgntTypesINFO{2}.ACT_Record{p1,p2,p3,p4,p5}.L1.spktrain;
data3= CnvrgntTypesINFO{3}.ACT_Record{p1,p2,p3,p4,p5}.L1.spktrain;
subplot(131);  imagesc(data1); 
subplot(132);  imagesc(data2); 
subplot(133);  imagesc(data3); 