function [] = plot_scatter( dataX_sumW, dataY_set, dataYLeg, tt1,xl,yl)
%UNTITLED6 �� �Լ��� ��� ���� ��ġ
%   �ڼ��� ���� ��ġ


    for pp = 1 : length(dataY_set)
        scatY = dataY_set{pp};
        scatter(dataX_sumW(:),scatY(:));     hold on;
    end
    xlabel(xl); ylabel(yl); title(tt1);
    if ~isempty(dataYLeg)
    legend(dataYLeg, 'Location','best');
    end

end

