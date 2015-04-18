function [ saveTxt ] = get_Parameters_saveText( PARAMETERS, paramLST, ID_LST )
%Return string of title text that include all specify parameter's number
%     PARAMETERS{ii}.PARAM 
%     PARAMETERS{ii}.lblTxt
%     PARAMETERS{ii}.titleTxt
%     PARAMETERS{ii}.saveTxt = eval(sprintf('saveTxt%d',ii));
saveTxt = '';
for ii = 1 : length(paramLST)
    prmNO = paramLST(ii);
    saveTxt = [saveTxt '_' PARAMETERS{prmNO}.saveTxt num2str(PARAMETERS{paramLST(ii)}.PARAM(ID_LST(ii))) ];
end          

end

