function [ Dstyl] = getDotStyle( id )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

DotStyleList = 'ho*.xsd^v><p+';
N= length(DotStyleList);
Dstyl = DotStyleList(mod(id,N)+1);

end

