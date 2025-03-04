function [X] = fft2_DC(x)
X=fftshift(fft2(ifftshift(x)));