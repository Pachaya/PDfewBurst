function [ val] = funcGauss_byDistance(maxval, dist, sig )
val = maxval*exp(- (dist.*dist)/(2.*sig.*sig));
end

