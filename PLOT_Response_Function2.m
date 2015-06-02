% Plot Response Function

for pp = 1 :  length(subPARAM)
    rf = figure;  set(gcf, 'position',[     372         -94        1160         753])
    ef = figure;  set(gcf, 'position',[     372         -94        1160         753])  
    LEG = cell(length(mainPARAM)*2,1);
    nR = length(subPARAM);
    for mm = 1 : length(mainPARAM)        
%         p4_ii = pp; p2_ii = ii;
%         p3_ii = 1;
        eval(sprintf('p%d_ii = mm; p%d_ii = pp; p%d_ii = dd; ',mainID, subID,dummy ));
        
        figure(rf);
        errorbar(PARAM1, avgFR(:,p2_ii,p3_ii,p4_ii,1),stdFR(:,p2_ii,p3_ii,p4_ii,1),[getDotStyle(mm) getLineStyle(1) getColorCode(mm)], 'LineWidth',2); hold on;
        errorbar(PARAM1, avgFR(:,p2_ii,p3_ii,p4_ii,2),stdFR(:,p2_ii,p3_ii,p4_ii,2),[getDotStyle(mm) getLineStyle(2) getColorCode(mm)], 'LineWidth',2);
        figure(ef);
        errorbar(PARAM1, avgEff(:,p2_ii,p3_ii,p4_ii,1),stdEff(:,p2_ii,p3_ii,p4_ii,1),[getDotStyle(mm) getLineStyle(1) getColorCode(mm)], 'LineWidth',2); hold on;
        errorbar(PARAM1, avgEff(:,p2_ii,p3_ii,p4_ii,2),stdEff(:,p2_ii,p3_ii,p4_ii,2),[getDotStyle(mm) getLineStyle(2) getColorCode(mm)], 'LineWidth',2);
        
        if(mainID == 4)
                
        TestW = PARAM4(p4_ii);  IGmean = IGmean_LST(1);
        Wscale = 0.001;
        if(IGmean == 0)
            Wspk = Wscale*TestW;
        else
            Wspk =  IGmean*-10*Wscale*TestW;
        end
        
        LEG{2*mm-1} =  ['[WT]   W_s_p_k  = ' num2str(Wspk)];
        LEG{2*mm} = ['[KO2]  W_s_p_k  = ' num2str(Wspk)];
        else
             tmpTxt =  get_Parameters_titleText(PARAMETERS,[subID, dummy],[pp, dd]);
             LEG{2*mm-1} =  ['[WT] ' tmpTxt];
             LEG{2*mm} = ['[KO2] ' tmpTxt];
        end
    end
    figure(rf);
    xlim([0 PARAM1(end)+10]); ylim([0 max(avgFR(:))+10])
    ylabel('Average Output fring rate (Hz)' ); xlabel('Input Firing rate(Hz)');
    title({'Response Function of each cell type', get_Parameters_titleText(PARAMETERS, [subID, dummy],[pp, dd]) })
    legend(LEG,'location','best');
    
    figure(ef);
    xlim([0 PARAM1(end)+5]); 
    ylabel('Efficacy (output FR /input FR)' ); xlabel('Input Firing rate(Hz)');
    title({'Response Function of each cell type', get_Parameters_titleText(PARAMETERS, [subID, dummy],[pp, dd]) })
    legend(LEG,'location','best');
    
    if(SAVE_FIG)
        tmpTxt = get_Parameters_saveText(PARAMETERS,[subID, dummy],[pp, dd]);
        fg = rf;  figName = ['ResFunc' tmpTxt];
        saveas(fg,[dirLoc dirFig figName '.fig'],'fig'); saveas(fg, [dirLoc dirFig figName '.jpg'],'jpg');
        fg = ef;  figName = ['Efficacy' tmpTxt];
        saveas(fg,[dirLoc dirFig figName '.fig'],'fig'); saveas(fg, [dirLoc dirFig figName '.jpg'],'jpg')
    end
end