function [] = plot_scatter_dotStyl( dataX_sumW, dataY_set, dataYLeg, tt1,xl,yl, dotStyl)
%UNTITLED6 이 함수의 요약 설명 위치
%   자세한 설명 위치


    for pp = 1 : length(dataY_set)
        scatY = dataY_set{pp};
        scatter(dataX_sumW(:),scatY(:),dotStyl);     hold on;
    end
    xlabel(xl); ylabel(yl); title(tt1);
    if ~isempty(dataYLeg)
    legend(dataYLeg, 'Location','best');
    end

end

