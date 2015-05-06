function [ fghandle ] = plotFigure_paramMat( data1,data2, tt1,tt2, xl,yl, xt,yt, suptitleTxt,figLoc )
%Return figure handle of figure with two subplots of parameter matrix
%Input
%   data1,dataK2, xl for xlabel() ,yl for ylabel(), xt,yt

caxis2  = max([ data1(:); data2(:)] );
fghandle = figure; set(fghandle, 'position',figLoc);   set(gcf,'PaperPositionMode','auto')
subplot(121); imagesc(data1);
set(gca, 'YTick', 1:length(yt), 'XTick', 1:length(xt))
set(gca, 'YTickLabel',  yt, 'XTickLabel', xt)
axis xy;
xlabel(xl); ylabel(yl); colorbar();  caxis([0 caxis2]); colormap('gray')
title(tt1);



txtshift = 0.3;
for ii = 1 : length(yt)
    for jj = 1 : length(xt)       
        tmpTxt = sprintf('%2.4f', data1(ii,jj));
        if(data1(ii,jj) > caxis2/2)
        text(jj-txtshift,ii,tmpTxt,'Color','k','FontSize', 12)
        else 
            text(jj-txtshift,ii,tmpTxt,'Color','w','FontSize', 12)
        end
    end
end



subplot(122); imagesc(data2);
set(gca, 'YTick', 1:length(yt), 'XTick', 1:length(xt))
set(gca, 'YTickLabel', yt, 'XTickLabel', xt)
axis xy; xlabel(xl); ylabel(yl); colorbar(); caxis([0 caxis2]); colormap('gray')
title(tt2);

for ii = 1 : length(yt)
    for jj = 1 : length(xt)       
         tmpTxt = sprintf('%2.4f', data2(ii,jj));
        if(data2(ii,jj) > caxis2/2)
        text(jj-txtshift,ii,tmpTxt,'Color','k','FontSize', 12)
        else 
            text(jj-txtshift,ii,tmpTxt,'Color','w','FontSize', 12)
        end
    end
end

suptitle(suptitleTxt)



% contour 
Ncontour = 4; 

caxis2  = max([ data1(:); data2(:)] );
fghandle2 = figure; set(fghandle2, 'position',figLoc);   set(gcf,'PaperPositionMode','auto')
subplot(121); contourf(data1,Ncontour)
set(gca, 'YTick', 1:length(yt), 'XTick', 1:length(xt))
set(gca, 'YTickLabel',  yt, 'XTickLabel', xt)
axis xy;
xlabel(xl); ylabel(yl); colorbar();  caxis([0 caxis2]); colormap('gray')
title(tt1);
subplot(122); contourf(data2,Ncontour)
set(gca, 'YTick', 1:length(yt), 'XTick', 1:length(xt))
set(gca, 'YTickLabel', yt, 'XTickLabel', xt)
axis xy; xlabel(xl); ylabel(yl); colorbar(); caxis([0 caxis2]); colormap('gray')
title(tt2);
suptitle(suptitleTxt)

end

