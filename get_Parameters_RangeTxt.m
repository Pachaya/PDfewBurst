function [ frangeTxt ] = get_Parameters_RangeTxt( PARAMETERS, paramLST )
%Return string of title text that include all specify parameter's number
%     PARAMETERS{ii}.PARAM 
%     PARAMETERS{ii}.lblTxt
%     PARAMETERS{ii}.titleTxt
%     PARAMETERS{ii}.saveTxt = eval(sprintf('saveTxt%d',ii));
frangeTxt = '';
for ii = 1 : length(paramLST)
    prmNO = paramLST(ii);
    frangeTxt = [frangeTxt  '_' PARAMETERS{prmNO}.saveTxt '_' num2str(PARAMETERS{paramLST(ii)}.PARAM(1))  '_' num2str(PARAMETERS{paramLST(ii)}.PARAM(end))];
end          

end

