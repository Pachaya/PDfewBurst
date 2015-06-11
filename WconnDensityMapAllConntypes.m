% Density plot for each connection types



gauss = load('TC_Conn_INFO_gaussPgaussW_03-Jun-2015');
unif = load('TC_Conn_INFO_avgPuniformW_03-Jun-2015');
negExp = load('TC_Conn_INFO_avgPnegecpW_03-Jun-2015');
PARAMETERS = gauss.PARAMETERS;
AllConnData = cell(3,1); 
AllConnData.gauss = gauss.TC_Conn_INFO;
tmp = AllConnData.gauss{1,1,1,1}; 
length(tmp.TC_basedOnM1)