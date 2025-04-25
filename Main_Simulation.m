%% Main simulation reconstruction code

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
%% Downloading phantom from .TIF file

n_media     = 1.5;         % refractive index of media
n_max       = 1.6;         % refractive index of max density feature
pdar        = 0;           % padding size before running forward model (to avoid edge artifacts)
objRaw      = makePhantom('phantom.tif',n_media,n_max);     % making phantom
objRaw      = padarray(objRaw,[50,50],n_media,"both");


%% Setting parameters relevant to physical object volume
save_recon  = true;
use_gpu     = true;                 % TRUE: use GPU device; FALSE: use computer CPU
use_field   = false;                % TRUE: uses field-component for reconstruction; FALSE: uses amplitude-only component for reconstruction

ps          = 0.1;                  % pixel size (x,y,z) in object space (micron)
lambda      = 0.635;                % central wavelength (micron)
NA          = 1.42;                 % numerical aperture of imaging and detection lens
n_imm       = n_media;              % refractive index of immersion media
n_m         = n_media;              % refractive index of media (for reconstruction purposes)
z_plane     = ps*0;                 % center plane of reconstruction volume, where 0 um is object volume center


obj         = single(objRaw-n_m);   % subtracting out global RI difference due to media
%% Setting spatial and frequency axes and propagation kernels

N           = size(obj,1)+2*pdar;   % lateral pixel dimension of padded object
x           = ps*[-N/2:N/2-1];      % 1D padded axis in x
[xx,yy]     = meshgrid(x,x);        % 2D padded grid in x/y

dfx         = 1/(N*ps);             % Fourier spacing of padded axis
fx          = dfx*[-N/2:N/2-1];     % 1D padded axis in fx
[fxx,fyy]   = meshgrid(fx,fx);      % 2D padded grid in fx/fy

fx          = ifftshift(fx);        % FFT shifting Fourier axes
fxx         = ifftshift(fxx);       % FFT shifting Fourier axes
fyy         = ifftshift(fyy);       % FFT shifting Fourier axes

% setting propagation kernels and pupil support
prop_phs            = 1i*2*pi*sqrt((n_imm/lambda)^2-(fxx.^2+fyy.^2));
NA_crop             = (fxx.^2 + fyy.^2 > (NA/lambda)^2);
prop_crop           = (fxx.^2 + fyy.^2 > (n_imm/lambda)^2);

% converting into GPU arrays if user targets gpu-enabling
if use_gpu
    obj             = gpuArray(obj);
    xx              = gpuArray(xx);
    yy              = gpuArray(yy);
    fxx             = gpuArray(fxx);
    fyy             = gpuArray(fyy);
    prop_phs        = gpuArray(prop_phs);
end

%% Setting parameters relevant to illumination angles

N_k             = 100;          % number of illumination-angle acquisitions
revs            = 6;            % number of revolutions the spiral takes
outerCirc       = true;         % boolean on whether to include outer circle points
N_o             = 30;           % if outerCirc == true, number of points in outer circle

% freq initialization of illumination angles
[kx_in,ky_in]   = generateSpiralPath(N_k, revs, outerCirc, N_o);    % generating normalized spiral coordinates

kx_in           = NA/lambda*kx_in;                                  % scaling fx by NA/lambda to span pupil function
ky_in           = NA/lambda*ky_in;                                  % scaling fy by NA/lambda to span pupil function

fx_in = ky_in;
fy_in = kx_in;

plot(fx_in,fy_in,'o:'); axis equal; axis tight;
title('illumination angle trajectory');
%% Running forward model on object phantom to simulate instrument measurements

obj_pad         = padarray(obj,[pdar,pdar,0],0);            % padding 3D phantom array
illumAngles     = 1:length(fx_in);
efield_acqs     = zeros([size(obj,1),size(obj,2),length(illumAngles)]);

if use_gpu
    efield_acqs    = gpuArray(efield_acqs);
end



[xposition,yposition]   = generateSpiralPath(100, 5, false, 0);    % generating normalized spiral coordinates

xposition = round(xposition*15,1);
yposition = round(yposition*15,1);

figure, plot(xposition,-yposition,'o:'); axis equal; axis tight;


tic;
for idx = illumAngles
    [efield,~]              = MultiSlice_Forward_acqdat_move(obj_pad, ps, xx, yy, dfx, ...
        prop_phs, NA_crop, lambda, ...
        fx_in(idx), fy_in(idx), z_plane, ...
        pdar,use_gpu,xposition(idx),yposition(idx));     % Multi-slice forward model


    efield_acqs(:,:,idx)    = efield;

    disp(['simulate data: ', num2str(idx)]);
end

disp('Simulated acquisitions via MSBP forward model is complete, and ready for reconstruction');
toc;


amp_acqs            = single(efield_acqs);          % incorporating both amp and phase

%% Translation trajectory initialization

% Choose translation trajectory initialization method (ImageCentroid OR ImageRegistration)
% ImageCentroid: Initialize translation trajectory using centroid of image
% ImageRegistration: Initialize translation trajectory using standard image registration

Translation_initial = 'ImageRegistration';  

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

end

% save the initial guess for comparison with the reconstruction results
xposition_guess_initial = xposition_recon;
yposition_guess_initial = yposition_recon;

figure, plot(xposition_guess_initial, yposition_guess_initial); axis equal; axis tight;
%% initializing forward model measurements and initial guess of reconstructed object

O               = 100;                   % axial dimension size of reconstruction space
psz             = ps;                  % pixel size (z) in reconstructed object space(micron)
% lateral pixel size is assumed to be same as variable 'ps'

reconObj        = single(0*randn([N,... % initialization of guess of reconstructed object (deltaRI, not RI), to be updated iteratively
    N,...
    O,]));

if use_gpu
    reconObj     = gpuArray(reconObj);
end
%% optimization params for iterative reconstruction

maxiter         = 200;                  % number of iterations to run optimization protocol for
step_size       = 0.5e-3;                % step size for gradient-based optimization protocol
Step_size_move  = 1e2;                  % setp size for translation trajectory optimization protocol

plot_range      = [-0.02,0.1];          % contrast to be used to show the reconstruction at each iteration
cost            = zeros(maxiter,1);     % cost function to evaluate convergence


reconObj_prox   = reconObj;             % used for Nesterov acceleration protocol for faster convergence
t_k             = 1;                    % parameter used for Nesterov acceleration
regParam        = 0.1e-3; %1e-3;                 % regularization parameter for 3D proxTV
%% initializing Figure windows to observe iterative process
close all;

% triframe cross-sectional views of the true phantom (to be used as a visual
% benchmark to evaluate convergence accuracy)
figure('Name','True Phantom (padded) RI difference');
MSBP_progview(real(obj_pad),1,plot_range);

% triframe cross-sectional views of the reconstructed object, as it goes
% iterative updates
figure('Name','Reconstruction result RI difference');
MSBP_progview(real(reconObj),2,plot_range,cost, 0)

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


    % Regularization across spatial domain
    try
        reconObj_prox1 = prox_tv3d(real(reconObj), regParam);
    catch
        disp('running regularizer on CPU because GPU ran out of memory for this memory-intensive procedure');
        reconObj_prox1  = prox_tv3d(gather(real(reconObj)), regParam);
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


%%
if save_recon
    save('recon.mat',...
        'reconObj','cost','xposition', 'yposition','xposition_recon','yposition_recon');
    disp('done saving file');
end
