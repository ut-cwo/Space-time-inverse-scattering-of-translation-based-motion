function outputImg = subpixelshift(inputImg, shift_X_in_Pixel, shift_Y_in_Pixel)
% Applying lateral shifting to the input image by muliplying the FT of the 
% input image with an additional phase term 

[mRow, nCol, nLayers]=size(inputImg);


inputImg = real(inputImg);
minval = min(inputImg(:));
inputImg = inputImg - minval;

[xGrid, yGrid] = meshgrid(1:mRow, 1:nCol);

xGrid              = gpuArray(xGrid);
yGrid              = gpuArray(yGrid);


AdditionalPhase = xGrid.*(shift_X_in_Pixel)/mRow + ...
    yGrid.*(shift_Y_in_Pixel)/nCol;

H = exp(-1j*2*pi*AdditionalPhase);

% Perform FFT on all layers at once
inputImgFT = fft2(inputImg);
inputImgFT = fftshift(fftshift(inputImgFT, 1), 2); % Applying fftshift separately to each dimension
outputImgFT = inputImgFT .* H;

% Perform inverse FFT and return the absolute value
outputImgFT = ifftshift(ifftshift(outputImgFT, 1), 2); % Applying ifftshift separately to each dimension
outputImg = abs(ifft2(outputImgFT))+minval;

end


    
