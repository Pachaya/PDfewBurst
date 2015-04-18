function [ titleTxt ] = get_Parameters_titleText( PARAMETERS, paramLST, ID_LST )
%Return string of title text that include all specify parameter's number
%     PARAMETERS{ii}.PARAM 
%     PARAMETERS{ii}.lblTxt
%     PARAMETERS{ii}.titleTxt
%     PARAMETERS{ii}.saveTxt = eval(sprintf('saveTxt%d',ii));
titleTxt = [PARAMETERS{paramLST(1)}.titleTxt ' = ' num2str(PARAMETERS{paramLST(1)}.PARAM(ID_LST(1))) ];
for ii = 2 : length(paramLST)
    prmNO = paramLST(ii);
    titleTxt = [titleTxt ', ' PARAMETERS{prmNO}.titleTxt ' = ' num2str(PARAMETERS{paramLST(ii)}.PARAM(ID_LST(ii))) ];
end          

end

