function [ ]= plot_paramMat( data1, tt1,xl,yl, xt,yt)
%Return figure handle of figure with two subplots of parameter matrix
%Input
%   data1,dataK2, xl for xlabel() ,yl for ylabel(), xt,yt

caxis2  = max( data1(:) );

%fghandle = figure; set(fghandle, 'position',figLoc);   set(gcf,'PaperPositionMode','auto')
%subplot(121); 

% plot to current axis
imagesc(data1);
set(gca, 'YTick', 1:length(yt), 'XTick', 1:length(xt))
set(gca, 'YTickLabel',  yt, 'XTickLabel', xt);
axis xy;
set(gca,'fontsize',11);
xlabel(xl); ylabel(yl); colorbar();  caxis([0 caxis2]); colormap('gray')
title(tt1);

flag = 0;
if(flag)
txtshift = 0.3;
for ii = 1 : length(yt)
    for jj = 1 : length(xt)       
        tmpTxt = sprintf('%2.2f', data1(ii,jj));
        if(data1(ii,jj) > caxis2/2)
        text(jj-txtshift,ii,tmpTxt,'Color','k','FontSize', 12)
        else 
            text(jj-txtshift,ii,tmpTxt,'Color','w','FontSize', 12)
        end
    end
end
end

end

