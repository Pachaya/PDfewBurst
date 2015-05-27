function [] = plot_contourparamMat( data1, tt1,xl,yl, xt,yt, Ncontour)
%Return figure handle of figure with two subplots of parameter matrix
%Input
%   data1,data2, xl for xlabel() ,yl for ylabel(), xt for XTicks ,yt for YTicks

caxis2  = max( data1(:) );
% contour 
%Ncontour = 4; 
if ~exist('Ncontour','var') || isempty(Ncontour)
  Ncontour=4;
end
contourf(data1,Ncontour)
set(gca, 'YTick', 1:length(yt), 'XTick', 1:length(xt))
set(gca, 'YTickLabel',  yt, 'XTickLabel', xt);
axis xy;
set(gca,'fontsize',11);
xlabel(xl); ylabel(yl); colorbar();  caxis([0 caxis2]); colormap('gray')
title(tt1);
end

