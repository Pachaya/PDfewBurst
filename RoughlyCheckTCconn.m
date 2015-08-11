
Fname = {'TC_Conn_INFO_gaussPgaussW_11-Jun-2015', 'TC_Conn_INFO_avgPuniformW_11-Jun-2015', 'TC_Conn_INFO_avgPnegexpW_11-Jun-2015'};

ii =1;
load(Fname{ii})
GaussTC = TC_Conn_INFO;
    
ii =2;
load(Fname{ii})
UnifTC = TC_Conn_INFO;

ii =3;
load(Fname{ii})
NegExpTC = TC_Conn_INFO;



rid = 5; wid = 6; 

figure; 

% for oscF = 1 : 2
%     for oscA = 1:3

        
sumWMat = zeros(5,6);
for rid = 1 : 5
for wid = 1: 6
tmpU = UnifTC{rid,wid,1,1,1};
nM = length(tmpU.TC_basedOnM1);
sumWperCell  = zeros(nM,1);
for ii  =1:nM
    sumWperCell(ii) = sum(tmpU.TC_basedOnM1{ii}.Weight);
end
% disp(mean(sumWperCell))
sumWMat(rid,wid) = mean(sumWperCell);
end
end

figure;
imagesc(sumWMat); axis xy; colormap('gray'); colorbar()

suptitle('Uniform')




sumWMat = zeros(5,6);
for rid = 1 : 5
for wid = 1: 6
tmpU = NegExpTC{rid,wid,1,1,1};
nM = length(tmpU.TC_basedOnM1);
sumWperCell  = zeros(nM,1);
for ii  =1:nM
    sumWperCell(ii) = sum(tmpU.TC_basedOnM1{ii}.Weight);
end
% disp(mean(sumWperCell))
sumWMat(rid,wid) = mean(sumWperCell);
end
end

figure; imagesc(sumWMat); axis xy; colormap('gray'); colorbar()
title('-Exp')



sumWMat = zeros(5,6);
for rid = 1 : 5
for wid = 1: 6
tmpU = GaussTC{rid,wid,1,1,1};
nM = length(tmpU.TC_basedOnM1);
sumWperCell  = zeros(nM,1);
for ii  =1:nM
    sumWperCell(ii) = sum(tmpU.TC_basedOnM1{ii}.Weight);
end
% disp(mean(sumWperCell))
sumWMat(rid,wid) = mean(sumWperCell);
end
end

figure; imagesc(sumWMat); axis xy; colormap('gray'); colorbar()
title('Gauss')


%% 
    
sumWMat = zeros(5,6);
for rid = 1 : 5
for wid = 1: 6
tmpU = UnifTC{rid,wid,1,1,1};
nM = length(tmpU.TC_basedOnM1);
sumWperCell  = zeros(nM,1);
for ii  =1:nM
    sumWperCell(ii) = length(tmpU.TC_basedOnM1{ii}.VL_ID);
end
% disp(mean(sumWperCell))
sumWMat(rid,wid) = mean(sumWperCell);
end
end

figure;
imagesc(sumWMat); axis xy; colormap('gray'); colorbar()

suptitle('Uniform')




sumWMat = zeros(5,6);
for rid = 1 : 5
for wid = 1: 6
tmpU = NegExpTC{rid,wid,1,1,1};
nM = length(tmpU.TC_basedOnM1);
sumWperCell  = zeros(nM,1);
for ii  =1:nM
    sumWperCell(ii) = length(tmpU.TC_basedOnM1{ii}.VL_ID);
end
% disp(mean(sumWperCell))
sumWMat(rid,wid) = mean(sumWperCell);
end
end

figure; imagesc(sumWMat); axis xy; colormap('gray'); colorbar()
title('-Exp')



sumWMat = zeros(5,6);
for rid = 1 : 5
for wid = 1: 6
tmpU = GaussTC{rid,wid,1,1,1};
nM = length(tmpU.TC_basedOnM1);
sumWperCell  = zeros(nM,1);
for ii  =1:nM
    sumWperCell(ii) = length(tmpU.TC_basedOnM1{ii}.VL_ID);
end
% disp(mean(sumWperCell))
sumWMat(rid,wid) = mean(sumWperCell);
end
end

figure; imagesc(sumWMat); axis xy; colormap('gray'); colorbar()
title('Gauss')
