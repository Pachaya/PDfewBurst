%Plot Neuron Network 3D
clc
% close all
clear all
strVL = 1; nVL =1150; nM1 = 166;
ssize = 1500;
halfsize = ssize/2;


Efile_VL = 'Neurons_location_VL_PDfewBurst_sample.txt';
Efile_M1 = 'Neurons_location_M1_PDfewBurst_sample.txt'; % only at center 
TCfile = 'TC_ConnectionWDParam_PDfewBurst_GPmVLmd1_0del_rTC250_specified_wmTC0.002_KO_GPmInput_Amp0.5_Dur1000_GPmVLw_m0.5_sig0_InGauss0.2_IGmean0_IGmeanSig0_W0.029_SpecifiedPoisSpk_sig0.00Hz_T4000_trial1.txt';
NNloc_VL= importdata(Efile_VL);  
NNloc_M1= importdata(Efile_M1); 

%for visualization purpose 
NNloc_M1(:,4) = 300+NNloc_M1(:,4);
NNloc_M1(:,2) = NNloc_M1(:,2)*2-halfsize ;
NNloc_M1(:,3) = NNloc_M1(:,3)*2-halfsize ;

figure;
hold on
scatter3(NNloc_VL(:,2),NNloc_VL(:,3),NNloc_VL(:,4),'r','MarkerFaceColor',[204/255 0 0]) %plot cells
scatter3(NNloc_M1(:,2),NNloc_M1(:,3),NNloc_M1(:,4),'b','MarkerFaceColor',[0 0 204/255]) %plot cells
%view(65,25)
view(-30,21)
grid on
axes_handle = gca;
axis equal
box on
%xlim([-10 110]);
%ylim([-10 110]);

NetConn = importdata(TCfile); 

totalNC = NetConn(1,1); NetConn = NetConn(2:end);
NetConn = reshape(NetConn,[5,totalNC]);
srclist = NetConn(1,:);
tarlist = NetConn(2,:); 
typelist = NetConn(3,:); % EE:0 EI:1 IE:2 II:3
weightlist = NetConn(4,:); % EE:0 EI:1 IE:2 II:3
delaylist = NetConn(5,:); % EE:0 EI:1 IE:2 II:3

srclist = srclist+1; %Because the neuron location variables (NNloc)use MATLAB indexin which start index at 1
tarlist= tarlist-nVL*3+1; % M1 cells' ID

% N_EE = length(find(typelist == 0))
% EE_NC = [srclist(find(typelist == 0)) tarlist(find(typelist == 0))];
% N_EI = length(find(typelist == 1))
% EI_NC = [srclist(find(typelist == 1)) tarlist(find(typelist == 1))];
% N_IE = length(find(typelist == 2))
% IE_NC = [srclist(find(typelist == 2)) tarlist(find(typelist == 2))];
% N_II = length(find(typelist == 3))
% II_NC = [srclist(find(typelist == 3)) tarlist(find(typelist == 3))];

%plot connection
hold on
if 0 %Plot all connection without differentiate connection type
    for i = 1:totalNC
        plot3([NNloc_VL(srclist(i),2) NNloc_M1(tarlist(i),2)],[ NNloc_VL(srclist(i),3) NNloc_M1(tarlist(i),3)],[ NNloc_VL(srclist(i),4) NNloc_M1(tarlist(i),4)],'Color',[1 0 1]);
    end
end


tarID  = unique(tarlist);
sampleM1ID = 87;
sampleNC = find( tarlist == tarID(sampleM1ID));
for i = 1 : length(sampleNC)
    nc = sampleNC(i);
        plot3([NNloc_VL(srclist(nc),2) NNloc_M1(tarlist(nc),2)],[ NNloc_VL(srclist(nc),3) NNloc_M1(tarlist(nc),3)],[ NNloc_VL(srclist(nc),4) NNloc_M1(tarlist(nc),4)],'Color','k');
end        
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
set(gca,'ZTickLabel','')

% %differentiate type
% for i = 1:N_EE
%         plot3([NNloc_VL(EE_NC(i,1),2) NNloc_VL(EE_NC(i,2),2)],[ NNloc_VL(EE_NC(i,1),3) NNloc_VL(EE_NC(i,2),3)],[ NNloc_VL(EE_NC(i,1),4) NNloc_VL(EE_NC(i,2),4)],'Color',[1 204/255 204/255]);
% end
% 
% for i = 1:N_EI
%         plot3([NNloc_VL(EI_NC(i,1),2) NNloc_VL(EI_NC(i,2),2)],[ NNloc_VL(EI_NC(i,1),3) NNloc_VL(EI_NC(i,2),3)],[ NNloc_VL(EI_NC(i,1),4) NNloc_VL(EI_NC(i,2),4)],'Color',[255/255 128/255 0/255]);
% end
% 
% for i = 1:N_IE
%         plot3([NNloc_VL(IE_NC(i,1),2) NNloc_VL(IE_NC(i,2),2)],[ NNloc_VL(IE_NC(i,1),3) NNloc_VL(IE_NC(i,2),3)],[ NNloc_VL(IE_NC(i,1),4) NNloc_VL(IE_NC(i,2),4)],'Color',[0 125/255 255/255]);
% end
% 
% for i = 1:N_II
%         plot3([NNloc_VL(II_NC(i,1),2) NNloc_VL(II_NC(i,2),2)],[ NNloc_VL(II_NC(i,1),3) NNloc_VL(II_NC(i,2),3)],[ NNloc_VL(II_NC(i,1),4) NNloc_VL(II_NC(i,2),4)],'Color',[102/255 178/255 255/255]);
% end

