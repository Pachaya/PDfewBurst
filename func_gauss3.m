function yy=func_gauss3 (maxval, dist, sig)

yy = maxval.*exp(- (dist.*dist)/(2.*sig.*sig));
end