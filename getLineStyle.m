function [ lTxt ] = getLineStyle(id)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

LineStyle = { '-.','-', '--', ':'};
N = length(LineStyle); 
lTxt = LineStyle{mod(id,N)+1};

end

