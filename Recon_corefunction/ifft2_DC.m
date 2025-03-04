function [x] = ifft2_DC(X)
x=fftshift(ifft2(ifftshift(X)));