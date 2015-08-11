% Plot Response Function
subPARAM = [20 40];
  rf = figure;  set(gcf, 'position',[     372         -94        1160         753])
    ef = figure;  set(gcf, 'position',[     372         -94        1160         753])  
    LEG = cell(length(mainPARAM)*length(subPARAM),1);
    nR = length(subPARAM);
    
    cnt = 0; 
for pp = 1 :  length(subPARAM)
  
    for mm = 1 : length(mainPARAM)        
%         p4_ii = pp; p2_ii = ii;
%         p3_ii = 1;
        eval(sprintf('p%d_ii = mm; p%d_ii = pp; p%d_ii = dd; ',mainID, subID,dummy ));
        CCI = length(mainPARAM)*pp +mm-1;
        figure(rf);
        errorbar(avgInputFR(:,p2_ii,p3_ii,p4_ii), avgFR(:,p2_ii,p3_ii,p4_ii,1),stdFR(:,p2_ii,p3_ii,p4_ii,1),[getDotStyle(mm) getLineStyle(pp) getColorCode(mm)], 'LineWidth',2); hold on;
%         errorbar(avgInputFR(:,p2_ii,p3_ii,p4_ii), avgFR(:,p2_ii,p3_ii,p4_ii,2),stdFR(:,p2_ii,p3_ii,p4_ii,2),[getDotStyle(mm) getLineStyle(2) getColorCode(pp)], 'LineWidth',2);
        figure(ef);
        errorbar(avgInputFR(:,p2_ii,p3_ii,p4_ii), avgEff(:,p2_ii,p3_ii,p4_ii,1),stdEff(:,p2_ii,p3_ii,p4_ii,1),[getDotStyle(mm) getLineStyle(pp) getColorCode(mm)], 'LineWidth',2); hold on;
%         errorbar(avgInputFR(:,p2_ii,p3_ii,p4_ii), avgEff(:,p2_ii,p3_ii,p4_ii,2),stdEff(:,p2_ii,p3_ii,p4_ii,2),[getDotStyle(mm) getLineStyle(2) getColorCode(pp)], 'LineWidth',2);
        
%         if(mainID == 4)
%                 
%         TestW = PARAM4(p4_ii);  IGmean = IGmean_LST(1);
%         Wscale = 0.001;
%         if(IGmean == 0)
%             Wspk = Wscale*TestW;
%         else
%             Wspk =  IGmean*-10*Wscale*TestW;
%         end
%         
%         LEG{2*mm-1} =  ['[WT]   W_s_p_k  = ' num2str(Wspk)];
%         LEG{2*mm} = ['[KO2]  W_s_p_k  = ' num2str(Wspk)];
%         else
cnt = cnt + 1;
             tmpTxt =  get_Parameters_titleText(PARAMETERS,[subID,mainID],[pp, mm]);
             LEG{cnt} =  ['[WT] ' tmpTxt];
%              cnt = cnt +1;
%              LEG{cnt} = ['[KO2] ' tmpTxt];
% %         end
    end   
end
 figure(rf);
    xlim([0 PARAM1(end)+5]); ylim([0 max(avgFR(:))+5])
    ylabel('Average Output fring rate (Hz)' ); xlabel('Input Firing rate(Hz)');
    title({'Response Function of each cell type', get_Parameters_titleText(PARAMETERS, [dummy],[dd]) })
    legend(LEG,'location','best');
    
    figure(ef);
    xlim([0 PARAM1(end)+5]); 
    ylabel('Efficacy (output FR /input FR)' ); xlabel('Input Firing rate(Hz)');
    title({'Response Function of each cell type', get_Parameters_titleText(PARAMETERS, [dummy],[ dd]) })
    legend(LEG,'location','best');
    
%     if(SAVE_FIG)
%         tmpTxt = get_Parameters_saveText(PARAMETERS,[subID, dummy],[pp, dd]);
%         fg = rf;  figName = ['CombineResFunc' tmpTxt];
%         saveas(fg,[dirLoc dirFig figName '.fig'],'fig'); saveas(fg, [dirLoc dirFig figName '.jpg'],'jpg');
%         fg = ef;  figName = ['CombineEfficacy' tmpTxt];
%         saveas(fg,[dirLoc dirFig figName '.fig'],'fig'); saveas(fg, [dirLoc dirFig figName '.jpg'],'jpg')
%     end
    