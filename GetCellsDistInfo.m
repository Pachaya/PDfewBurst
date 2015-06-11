function [varargout ] = GetCellsDistInfo( locTxt , PLOT)
%Return Matrix for position of the specific neural mosaics
%   Detailed explanation goes here



if(nargin == 1)
    PLOT = 1;
end

NNloc= importdata([locTxt '.txt']); 
ncells = NNloc(1); NNloc = NNloc(2:end);

Cellpos = reshape(NNloc,2,ncells)';
if(PLOT)
figure; scatter(Cellpos(:,1),Cellpos(:,2),20); axis square;
xlabel('um');  ylabel('um'); title('Distribution of cell position');
end

Distmat = squareform (pdist(Cellpos));
tmpnn = sort(Distmat,2);
nnDist = tmpnn(:,2);

if(PLOT)
figure; hist(nnDist,25); 
xlabel('Nearest Neighbor Distance (um)');
ylabel('# of Cells'); 
title({'Distribution of nearest neighbor distance', ['Mean = ' num2str(mean(nnDist)) ', std = ' num2str(std(nnDist))]});
end


switch nargout
    case 1 
        varargout = {Cellpos};
    case 2
        varargout = {Cellpos, Distmat};
    case 3
        varargout = {Cellpos, Distmat, nnDist};
    otherwise 
          varargout = {Cellpos, Distmat, nnDist};
end

end

