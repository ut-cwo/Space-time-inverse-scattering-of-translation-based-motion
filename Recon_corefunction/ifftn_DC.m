function [x] = ifftn_DC(X)
x=fftshift(ifftn(ifftshift(X)));