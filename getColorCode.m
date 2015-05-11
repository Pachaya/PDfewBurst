function [cTxt] = getColorCode(id)
Ccode = 'yrbkmgc';
N = length(Ccode);
cTxt = Ccode(mod(id,N)+1);

end

