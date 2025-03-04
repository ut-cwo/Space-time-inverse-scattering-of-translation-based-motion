function [X] = fftn_DC(x)
X=fftshift(fftn(ifftshift(x)));