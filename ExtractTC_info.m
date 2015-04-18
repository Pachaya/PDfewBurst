function [ varargout ] = ExtractTC_info( dirLoc, Name_postfix, display  )
%Get information of TC connection of the specified Name_postfix 
%The output cell ID based on NEURON's system -->  ID of the cell start from 0 so the ID should be add by one, the M1 cell start from 3450 and get only center part of M1
%Note on output arguments
%   #1      [TC_basedOnM1 ]                                           = ExtractTC_info( dirLoc, Name_postfix, display  );
%   #2      [TC_basedOnM1, TCconn_raw ]                               = ExtractTC_info( dirLoc, Name_postfix, display  );
%   #3      [TC_basedOnM1, TC_sumW, TC_maxW ]                         = ExtractTC_info( dirLoc, Name_postfix, display  );
%   #4      [TC_basedOnM1, TC_sumW, TC_maxW, TC_numVL ]               = ExtractTC_info( dirLoc, Name_postfix, display  );
%   #5      [TC_basedOnM1, TC_sumW, TC_maxW, TC_numVL, TCconn_raw ]   = ExtractTC_info( dirLoc, Name_postfix, display  );

tmpData = importdata([dirLoc 'TC_ConnectionWDParam_' Name_postfix '.txt'] );
N_TC = tmpData(1); tmpData = tmpData(2:end);  tmpData = reshape(tmpData,[5,N_TC])';
srcVLcell = tmpData(:,1);  tarM1Cell = tmpData(:,2); connType = tmpData(:,3); connWeight = tmpData(:,4); connDelay = tmpData(:,5); % ID of the cell start from 0 so the ID should be add by one, the M1 cell start from 3450 and get only center part of M1

% Find the set of target cell
[tarID_M1cell, IA, IC] = unique(tarM1Cell,'stable'); % tarID_M1cell = tarM1Cell(IA) and tarM1Cell = tarID_M1cell(IC)

% for each M1 cell that has got feedforward connection from VL
TC_basedOnM1 = cell(1,length(tarID_M1cell));
TC_sumW = zeros(1,length(tarID_M1cell));
TC_maxW = zeros(1,length(tarID_M1cell));
TC_numVL = zeros(1,length(tarID_M1cell));
for ii = 1: length(tarID_M1cell)
    
    tmpConnID = find(tarM1Cell == tarID_M1cell(ii));
    tmpCell = [];
    tmpCell.M1_ID = tarID_M1cell(ii);
    tmpCell.Conn_ID = tmpConnID;
    TC_numVL(ii) = length( tmpConnID);
    tmpCell.VL_ID = srcVLcell(tmpConnID);     tmpCell.Weight = connWeight(tmpConnID); tmpCell.Delay = connDelay(tmpConnID);
    TC_sumW(ii) =  sum(tmpCell.Weight); TC_maxW(ii) =  max(tmpCell.Weight);
    TC_basedOnM1{ii} = tmpCell;
end


TCconn_raw.N_TC = N_TC; % raw data
TCconn_raw.srcVLcell= srcVLcell;
TCconn_raw.tarM1Cell= tarM1Cell;
TCconn_raw.connType= connType;
TCconn_raw.connWeight= connWeight;
TCconn_raw.connDelay= connDelay;


switch nargout
    case 1
        varargout{1} = TC_basedOnM1;
    case 2
        varargout{1} = TC_basedOnM1;
        varargout{2} = TCconn_raw;
    case 3
        varargout{1} = TC_basedOnM1;
        varargout{2} =  TC_sumW;
        varargout{3} = TC_maxW;
    case 4
        varargout{1} = TC_basedOnM1;
        varargout{2} = TC_sumW;
        varargout{3} = TC_maxW;
        varargout{4} = TC_numVL ;
        case 5
        varargout{1} = TC_basedOnM1;
        varargout{2} = TC_sumW;
        varargout{3} = TC_maxW;
        varargout{4} = TC_numVL ;
        varargout{5} = TCconn_raw;
    otherwise
        disp('Unknown number of output arguments')
end

if( nargin == 3) && (display)
    
        disp(['Number of M1 cell = ' num2str(length(TC_basedOnM1))]);
        disp(['Average # of VL per M1 = ' num2str(mean(TC_numVL))])
        disp(['average weight summation of thalamocortical connection = ' num2str(mean(TC_sumW) )])
        disp(['average maximum value of weight = ' num2str(mean(TC_maxW) )])
    
end
    
    


end

