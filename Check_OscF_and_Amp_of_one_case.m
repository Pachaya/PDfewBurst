%% Check Osc F and Amp of one case
% OCS F -->PARAM4, OSC AMP --> PARAM5
p1_ii = 1; p2_ii =1; p3_ii = 1; p6_ii = 1;
nR = length(PARAM4); nC = length(PARAM5);
fosc = figure;  set(fosc, 'position', [   147         799        2094         471]); set(gcf,'PaperPositionMode','auto');
cnt1 =0;
fg_fft = figure;  set(fg_fft,'position',[ 458         600        1421         568]); set(gcf,'PaperPositionMode','auto');
cntcnt = 0;

for p4_ii = 1 : length(PARAM4)  % OSC F
    for p5_ii = 1 : length(PARAM5) % OSC Amp
        
        tmpTxt = get_Parameters_titleText(PARAMETERS,[4, 5],[p4_ii,p5_ii]);
        %%% for VL
        tmpWTtrain_VL =  ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii,p6_ii}.VL.WT.All.spktrain;
        tmpKOtrain_VL =  ACT_Record{p1_ii,p2_ii,p3_ii,p4_ii,p5_ii,p6_ii}.VL.KO.All.spktrain;
        
        % Check Oscillation Frequency of VL
        
        T = TSTOP;
        tt = 0:1:T;
        
        % Spike Train
        [KOci, KOti] = find(tmpKOtrain_VL);
        [WTci, WTti] = find(tmpWTtrain_VL);
        
        % mean firing rate
        KOfmu = mean(tmpKOtrain_VL, 1)*1000;
        WTfmu = mean(tmpWTtrain_VL, 1)*1000;
        %         fosc = figure;  set(fosc, 'position', [496         763        1042         420]); set(gcf,'PaperPositionMode','auto');
        %         subplot(211);
        %         plot(WTfmu,'k');   hold on;  plot(KOfmu,'r'); legend('WT','KO'); title('Raw instantaneous average firing rate')
        
        gsig = 3;   % filter [ms]
        gfil = exp( -[-round(4*gsig):round(4*gsig)].^2 /2/gsig^2);
        gfil = gfil/sum(gfil(:));
        
        KOfmus = conv(KOfmu, gfil, 'same');
        WTfmus = conv(WTfmu, gfil, 'same');
        %         subplot(212); plot(WTfmus,'k');   hold on;  plot(KOfmus,'r'); legend('WT','KO'); title('Smooth instantaneous average firing rate');
        %         suptitle(tmpTxt)
        cnt1 = cnt1+1;
        figure(fosc); subplot(nR,nC,cnt1);
        plot(WTfmus,'k');   hold on;  plot(KOfmus,'r'); legend('WT','KO','location','best'); title(tmpTxt);  %title('Smooth instantaneous average firing rate');
        
        % FFT
        KOfft = abs(fftshift(fft(KOfmus)));
        WTfft = abs(fftshift(fft(WTfmus)));
        
        
        dHz = 1/T * 1000;
        Hz_ind = [-round(T/2):round(T/2)]*dHz;
        %         Hz_ind = Hz_ind(1:T);
        %         fg_fft = figure;  set(fg_fft,'position',[ 458         600        1421         568]); set(gcf,'PaperPositionMode','auto');
        %         subplot(121);
        figure(fg_fft);  cntcnt = cntcnt +1;
        subplot(nR,nC,cntcnt);
        plot(Hz_ind, WTfft,'k'); hold on;  plot(Hz_ind, KOfft,'r');  legend('WT','KO');
        xlim([1 100]); ylabel('FFT amplitude (a.u.)'); xlabel('Frequency (Hz)'); title(tmpTxt); %title('raw FFT')
        
        [Hind, KO_FFT, WT_FFT] = analyze_fft(Hz_ind, WTfft, KOfft);
        %         subplot(122);
        %         plot(Hind, WT_FFT.FF,'k') ;     hold on; plot(Hind, KO_FFT.FF,'r');  legend('WT','KO');
        %         xlim([1 100]); ylabel('FFT amplitude (a.u.)'); xlabel('Frequency (Hz)'); title('Smooth FFT');
        %         suptitle(tmpTxt);
        
    end
end
tmpTxt = get_Parameters_titleText(PARAMETERS, [1:3, 6], [p1_ii p2_ii p3_ii p6_ii ]);
figure(fosc); suptitle({['Smooth instantaneous average firing rate (Gaussian filter size = ' num2str(gsig) ')'], tmpTxt});
figure(fg_fft); suptitle({'Raw FFT', tmpTxt});
if (SAVE_FIG)
    tmpTxt = get_Parameters_saveText(PARAMETERS, [1:3 6], [p1_ii p2_ii p3_ii p6_ii]);
    ffig = [ dirLoc dirFig 'ChecksampleVLInst_Avg_FR_' tmpTxt ];
    saveas(  fosc, [ffig '.jpg'], 'jpg');        saveas(  fosc, [ffig '.fig'], 'fig');
    ffig = [ dirLoc dirFig 'ChecksampleVL_Osc_fft_' tmpTxt ];
    saveas(  fg_fft, [ffig '.jpg'], 'jpg');        saveas(  fg_fft, [ffig '.fig'], 'fig');
end
