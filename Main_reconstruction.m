%% Main reconstruction code

% Reference for MSBP:
% S. Chowdhury, M. Chen, R. Eckert, D. Ren, F. Wu, N. Repina, and L. 
% Waller, "High-resolution 3D refractive index microscopy of multiple-
% scattering samples from intensity images," Optica 6, 1211-1219 (2019) 
%
% Acknowledgements:
% This MSBP demonstration utilizes regularization using the 3D total 
% variation (TV) proximal operator. We gratefully acknowledge the work done 
% by the UNLocBoX team, which has provided an extensive open-source 
% compilation of several optimizers utilizing proximity operators. From
% the UNLocBoX toolbox, we use the 'prox_tv1d' and 'prox_tv3d' function for 
% 3D TV prox regularization. Full documentation of UNLocBoX can be found at:
%
% UNLocBoX toolbox:     https://epfl-lts2.github.io/unlocbox-html/
% download this toolbox and put into the same path as this .M file 
%
% Reference for 'prox_tv3d':  
% A. Beck and M. Teboulle. Fast gradient-based algorithms for constrained 
% total variation image denoising and deblurring problems. Image 
% Processing, IEEE Transactions on, 18(11):2419--2434, 2009

%%

clear; close all;
clc;

%% adding relevant paths

folderPath1 = 'Recon_corefunction';
folderPath2 = 'Recon_corefunction\unlocbox';

addpath(folderPath1);
addpath(folderPath2);
init_unlocbox();

%% Downloading scattering data and imaging parameters from .MAT file

load('Phantom_data.mat');       % Amplitude: Normalized scattering field amplitude data 

%% Setting parameters relevant to physical object volume (some are drawn directly from data .MAT file)

save_recon  = true;                 % True: save results; FALSE: do not save results
use_gpu     = true;                 % TRUE: use GPU device; FALSE: use computer CPU
use_field   = false;                % TRUE: uses field-component for reconstruction; FALSE: uses amplitude-only component for reconstruction

ps          = ps;                   % pixel size (x,y,z) in object space (micron)
lambda      = lambda;               % central wavelength (micron)
NA          = NA;                   % numerical aperture of imaging and detection lens
n_m         = n_m;                  % refractive index of immersion media
n_imm       = n_imm;                % refractive index of media (for reconstruction purposes)
pdar        = 0;                    % padding size to avoid edge artifacts
z_plane     = ps*0;                 % center plane of reconstruction volume, where 0 um is object volume center

%% Setting up the field-of-views (FOV) to divide up complete field-of-view

FOV_size   = 550;

x_start = 80;
y_start = 100;

rows = x_start:x_start+FOV_size-1;
cols = y_start:y_start+FOV_size-1;

figure(1)
imagesc(Amplitude(:,:,1)); axis equal; colormap gray; axis tight; clim([0 2]);
hold on;
rectangle('Position',[cols(1),rows(1),FOV_size,FOV_size], 'LineWidth',3, 'EdgeColor','r')
hold off;

amp_acqs    = Amplitude(rows,cols,:);

%% Setting spatial and frequency axes and propagation kernels

N           = size(amp_acqs,1)+2*pdar;  % lateral pixel dimension of a padded patch within the object
x           = ps*(-N/2:N/2-1);          % 1D padded axis in x
[xx,yy]     = meshgrid(x,x);            % 2D padded grid in x/y

dfx         = 1/(N*ps);                 % Fourier spacing of padded axis
fx          = dfx*(-N/2:N/2-1);         % 1D padded axis in fx
[fxx,fyy]   = meshgrid(fx,fx);          % 2D padded grid in fx/fy

fx          = ifftshift(fx);            % FFT shifting Fourier axes
fxx         = ifftshift(fxx);           % FFT shifting Fourier axes
fyy         = ifftshift(fyy);           % FFT shifting Fourier axes

% setting propagation kernels and pupil support
prop_phs            = 1i*2*pi*sqrt((n_imm/lambda)^2-(fxx.^2+fyy.^2));
NA_crop             = (fxx.^2 + fyy.^2 > (NA/lambda)^2);

% converting into GPU arrays if user targets gpu-enabling
if use_gpu
    xx              = gpuArray(xx);
    yy              = gpuArray(yy);
    fyy             = gpuArray(fyy);
    fyy             = gpuArray(fyy);
    prop_phs        = gpuArray(prop_phs);
    amp_acqs        = gpuArray(amp_acqs);
end

%% Downloading illumination k-vectors from .MAT file and accounting for system scan-angle orientation

fx_in           = -fx_illum;
fy_in           = -fy_illum;

%% Translation trajectory initialization

% Choose translation trajectory initialization method (ImageCentroid OR ImageRegistration)
% ImageCentroid: Initialize translation trajectory using centroid of image
% ImageRegistration: Initialize translation trajectory using standard image registration
% StageFeedback: Initialize translation trajectory using motorized stage feedback

Translation_initial = 'ImageCentroid';  

switch Translation_initial

    case 'ImageCentroid'

        amp_acqs_backsub  = gather(abs(amp_acqs-1));

        xposition_recon = zeros(1,size(amp_acqs,3));
        yposition_recon = zeros(1,size(amp_acqs,3));

        for i = 1:size(amp_acqs,3)

            Binarized_img = imbinarize(amp_acqs_backsub(:,:,i));

            totalMass = sum(Binarized_img(:));
            [rows, cols] = size(Binarized_img);
            [X, Y] = meshgrid(1:cols, 1:rows);

            % Calculate the centroid of image
            xposition_recon(i) = sum(sum(X.* Binarized_img))/totalMass;
            yposition_recon(i) = sum(sum(Y.* Binarized_img))/totalMass;

        end

        xposition_recon = xposition_recon-xposition_recon(1);
        yposition_recon = yposition_recon-yposition_recon(1);

    case 'ImageRegistration'

        usfac = 10;

        x_shift = zeros(1,size(amp_acqs,3)-1);
        y_shift = zeros(1,size(amp_acqs,3)-1);

        for i = 1:size(amp_acqs,3)-1

            image1 = amp_acqs(:,:,i);
            image2 = amp_acqs(:,:,i+1);

            [output, Greg] = dftregistration(fft2(image1),fft2(image2),usfac);
            x_shift(i) = -output(4);
            y_shift(i) = -output(3);
        end

        xposition_recon = zeros(1,size(amp_acqs,3));
        yposition_recon = zeros(1,size(amp_acqs,3));

        for i = 1:length(x_shift)
            if i == 1
                xposition_recon(i+1) = x_shift(i);
                yposition_recon(i+1) = y_shift(i);
            else
                xposition_recon(i+1) = xposition_recon(i)+x_shift(i);
                yposition_recon(i+1) = yposition_recon(i)+y_shift(i);
            end
        end

    case 'StageFeedback'

        xposition_recon = xposition_stage;
        yposition_recon = yposition_stage;
end

% save the initial guess for comparison with the reconstruction results
xposition_guess_initial = xposition_recon;
yposition_guess_initial = yposition_recon;

figure, plot(xposition_guess_initial, yposition_guess_initial); axis equal; axis tight;

%% initializing forward model measurements and initial guess of reconstructed object

O               = 120;                  % axial dimension size of reconstruction space
psz             = 0.4;                  % pixel size (z) in reconstructed object space(micron)

reconObj        = single(0*randn([N,... % initialization of guess of reconstructed object (deltaRI, not RI), to be updated iteratively
    N,...
    O,]));

if use_gpu
    reconObj     = gpuArray(reconObj);
end

%% optimization params for iterative reconstruction

maxiter         = 200;                  % number of iterations to run optimization protocol for

step_size       = 1e-4;                 % step size for gradient-based optimization protocol
Step_size_move  = 2e2;                  % setp size for translation trajectory optimization protocol

plot_range      = [-0.02,0.04];         % contrast to be used to show the reconstruction at each iteration
cost            = zeros(maxiter,1);     % cost function to evaluate convergence

reconObj_prox   = reconObj;             % used for Nesterov acceleration protocol for faster convergence
t_k             = 1;                    % parameter used for Nesterov acceleration

regParam_spatial= 1.5e-3;               % regularization parameter for 3D proxTV across spatial domain
regParam_time   = 1e-2;                 % regularization parameter for 1D proxTV across time 

%% initializing Figure windows to observe iterative process

close all;

% triframe cross-sectional views of the reconstructed object, as it undergoes
% iterative updates

figure('Name','Reconstruction result');
figNum = 1;
MSBP_progview(real(reconObj),figNum,plot_range,cost, 0)

pause(0.01);

%% Running iterative optimization of object volume. 
% Variable 'reconObj' is the final 3D refractive-index reconstruction and xposition_recon and yposition_recons are final translation trajectory

gpu_fail_num = 0;

tic;

seq = 1:1:length(fx_in);

for iter = 1:maxiter

    pause(0.01);


    for illum_angle = 1:length(fx_in)

        % sample translation
        reconObj = subpixelshift(real(reconObj), xposition_recon(seq(illum_angle)), yposition_recon(seq(illum_angle)));

        % compute estimated exit field on the camera plane
        [efield,efield_vol,U_in2]     = MultiSlice_Forward(reconObj, psz, xx, yy, dfx, prop_phs, NA_crop, lambda, fx_in(seq(illum_angle)), fy_in(seq(illum_angle)), z_plane, pdar, use_gpu);

        % compute gradient (and update refractive index at each layer and translation trajectory)
        [reconObj,funcVal,xposition_update,yposition_update]      = BPM_update_move(reconObj, psz, efield, efield_vol, amp_acqs(:,:,seq(illum_angle)), prop_phs, NA_crop, lambda, z_plane, step_size, pdar, use_field, U_in2,fxx,fyy,xposition_recon(seq(illum_angle)),yposition_recon(seq(illum_angle)),seq(illum_angle),Step_size_move);

        % move sample back to center
        reconObj = subpixelshift(real(reconObj), -xposition_recon(seq(illum_angle)), -yposition_recon(seq(illum_angle)));

        % replace translation trajectory value
        xposition_recon(seq(illum_angle)) = xposition_update;
        yposition_recon(seq(illum_angle)) = yposition_update;


        % compute accumulated error for current iteration
        cost(iter)         = cost(iter) + gather(funcVal);
        fprintf('illum_angle: %1.0d  iteration: %1.0d\n',illum_angle,iter)
    end

    % apply non-negativity constraint
    Positive_mask = reconObj<0;
    reconObj(Positive_mask) = 0;

    % Prox operator is a memory-intensiver operator. If GPU crashes due to
    % memory requirements, use CPU instead. It will be slower but the
    % program won't crash.

    % Regularization across time domain
    xposition_guess_TV = prox_tv1d(xposition_recon', regParam_time);
    yposition_guess_TV = prox_tv1d(yposition_recon', regParam_time);

    xposition_recon = xposition_guess_TV';
    yposition_recon = yposition_guess_TV';

    % Regularization across spatial domain
    try
        reconObj_prox1 = prox_tv3d(real(reconObj), regParam_spatial);
    catch
        disp('running regularizer on CPU because GPU ran out of memory for this memory-intensive procedure');
        reconObj_prox1  = prox_tv3d(gather(real(reconObj)), regParam_spatial);
        reconObj_prox1  = gpuArray(reconObj_prox1);
        gpu_fail_num    = gpu_fail_num+1;   % counter to keep track of how many times GPU failed to regularize
    end

    if iter>1
        if cost(end) > cost(end-1)
            t_k   = 1;
            reconObj = reconObj_prox;
            continue;
        end
    end

    % Nesterov's update
    t_k1       = 0.5 * (1 + sqrt(1 + 4 * t_k^2));
    beta       = (t_k - 1)/t_k1;
    reconObj   = reconObj_prox1 + beta*(reconObj_prox1 - reconObj_prox);
    t_k        = t_k1;
    reconObj_prox = reconObj_prox1;
    fprintf('iteration: %d, error: %5.5e, elapsed time: %5.2f seconds\n',iter, cost(iter));

    MSBP_progview(real(reconObj), figNum, plot_range, cost, iter)
    pause(0.01);

end
toc;
close_unlocbox();


%% In case you want to save reconstruted data and relevant parameters
if save_recon
    save('Recon_results.mat',...
        'reconObj','xposition_guess_initial','yposition_guess_initial','xposition_recon','yposition_recon');
    disp('done saving file');
end