close all
clear;clc

% !!! Download 'bruker_data.mat' from the BassData folder on Box so you
% don't have to run RussRecon

% This script loads in k space data and displays a slice of the image data.
% It then uses a sampling mask to generate undersampled k space data and
% displays three different types of reconstructions: zero-filled recon
% (i.e. directly reconstructing the data), SENSE reconstruction
% (autocalibrated), and ESPIRiT reconstruction (using ESPIRiT calibration)

% BART must be installed for bart commands to work
% %! means these settings must be changed for each user depending on where
% your code/BART path is
% %* means you can toggle these options

main_dir = '/Users/janetchen/Documents/Bass Connections'; %!
bartpath = '/Applications/bart/matlab'; %!*
setenv('TOOLBOX_PATH','/Applications/bart') %!

% The following packages are needed in your main directory:
paths = {main_dir;[main_dir '/STI_Suite_v2/STI_Suite_v2.1/Source Code_v2.1'];...
    [main_dir '/RussRecon'];[main_dir '/ImageProcessing'];...
    [main_dir '/NIfTI_20140122'];[main_dir '/MultipleEchoRecon'];...
    [main_dir '/XCalc']};

% If you need to extract the image data, you need to have the .work folder
% in the main directory. If you're loading it in (and bruker_img.MAT should
% be in the addecode main directory), you don't need this folder
workpath = [main_dir '/B04027.work'];

use_espirit_data = false; %* try out the example image or use our image?
% If you want to use the example image, the BART folder at
% https://github.com/mikgroup/espirit-matlab-examples must be downloaded
% (or just the 'full.cfl' from the 'data' folder). Put the correct path
% name here (don't include the '.cfl' at the end):
fullpath = 'BART/espirit-matlab-examples-master/data/full';

load_data = true; %* 'true' if you have the k space data saved and don't
% want to run get_RussRecon_img again
save_data = false; %* 'true' if you want to save the data
data_path = 'bruker_data.mat'; % image path if it's already been saved

bart_mask = true; %* if generating mask, use variable density Poisson mask
% from BART code (true), or code from BJ Anderson (false)? Right now, the
% undersampling factor can be better controlled by the Anderson code
load_mask = false; %* load sampling mask in? Currently, only BJ mask is saved
num_iter = 1; % Number of sampling patterns/reconstructions to generate
% (to get mean RLNE, PSNR values)

addpath(bartpath)
addpath(main_dir)
if ~use_espirit_data && ~load_data % Need to be in workpath for RussRecon
    cd(workpath)
    % Add paths so that RussRecon can be used
    for ii = 1:length(paths)
        addpath(paths{ii,1})
    end
end

if use_espirit_data
    % Read in k space data from BART's example (fully sampled)
    ksp_data = readcfl(sprintf('%s/%s',main_dir,fullpath));
else
    if load_data
        % Load in previously-extracted k space data
        load(sprintf('%s/%s',main_dir,data_path))
    else
        % Get k space data from .work directory
        ksp_data=get_RussRecon_img('bruker','center','save');
        % If data should be saved, save it
        if save_data
            save(sprintf('%s/%s',main_dir,data_path),'bruker_data')
        end
    end
end

% If you want to write the data to a cfl file for BART
%cflpath = sprintf('%s/img',paths{2,1});
%writecfl(cflpath,img)

% The MRI-specific BART commands assume that the first three dimensions
% represent space/spatial k-space, the next two dimensions represent coils
% and ESPIRiT maps, and the 10th dimension (zero-indexing) represents
% temporal phases.

num_coils = size(ksp_data,4);
if use_espirit_data
    % img is 1x230x180x8
    % Swap x, y, z dimensions so 3rd dimension is 1
    ksp_data_echo1 = permute(ksp_data(1,:,:,1:num_coils),[2,3,1,4]); 
    % Take inverse Fourier transform
    coilimg = bart('fft -i 6', ksp_data);
    % Take the root sum of squares to produce one image from multiple coils
    rss = bart('rss 4',squeeze(coilimg));
    slice = 1; % Original X dimension is 1 in length, so use that 'slice'
    sz_x = size(ksp_data,2); sz_y = size(ksp_data,3);
else
    % Z-axis slice to look at
    slice = 52;
    % 8 echoes, 2 coils
    % Original k space data dimensions: x, y, z, # echoes, # coils
    % ksp_data matrix is 192x192x90x8x2
    sz_x = size(ksp_data,1); sz_y = size(ksp_data,2);
    
    % To match espirit 'full' img dimensions (XxYxcoils, 230x180x8), slice
    % in z-dimension IN IMAGE domain
    
    % Inverse Fourier transform the k space data to get the image data
    % 'squeeze' just removes singleton dimensions (i.e. dimensions of
    % length 1)
    temp = bart('fft -i 7',squeeze(ksp_data(:,:,:,1,:)));
    % Another method to get the ifft: ifftnc(<matrix>);
    
    % Take the z-axis slice
    coilimg = temp(:,:,slice,:);
    % Fourier transform back to k space
    % ksp_data_echo1 here is 192x192x1x2 (XxYxZxcoils). This will be used
    % later to get the undersampled data
    ksp_data_echo1 = bart('fft 7',coilimg);
    
    % Take the root sum of squares to produce one image from multiple coils
    rss = bart('rss 4', squeeze(coilimg));
end

%% Generate sampling pattern and undersample K-space

% Coil # for which to plot k space data and undersampling
coil = 1;

if bart_mask
    % bart poisson -Y $dim -Z $dim -y $yaccel -z $zaccel -C $caldim -v -e mask
    % -Y: dim 1. -Z: dim 2. -v: variable density (which the workshop uses)
    % -e: elliptical scanning
    % Default values: Y: 192, Z: 192, y: 1.5, z: 1.5, C: 32
    % Change C to 20?
    % ~25% of elements in und2x2 in espirit.m were nonzero
    % !!! Although this should reflect 9216 points, only 8005 points are
    % generated
    num_points = sz_x*sz_y*0.25;
    
    % Sampling pattern. Binary mask (ones and zeros), where ones are the
    % only points that we sample
    % Original output of BART command is 3 dimensional, with first
    % dimension of length 1, so squeeze to remove this dimension
    sampling_pattern = squeeze(bart(sprintf('poisson -Y %d -Z %d -y 1 -z 1 -C 60 -R %d -v -e',sz_x,sz_y,num_points)));% -v -e');
else
    % Mask from BJ Anderson's code
    if load_mask
        load('bj_sampling_pattern.mat')
    else
        % Inputs
        acceleration = 4; %* Undersample by a factor of <acceleration>
        pa = 2.3;
        pb = 5.6;
        
        % Sampling pattern. Logical binary mask (ones and zeros), where ones
        % indicate the only points that we sample
        sampling_pattern = sampling_mask(acceleration,sz_x,sz_y,pa,pb);
    end
end

% Show how much of the data is kept
fprintf(sprintf('\tMask non-zero percentage: %.2f%%',length(find(sampling_pattern ~= 0))/numel(sampling_pattern)*100))
fprintf('\n')

% Zero out unsampled areas through element-by-element multiplication with
% sampling pattern
% 'us' means 'undersampled'
us_ksp_data_echo1 = bart('fmac',squeeze(ksp_data_echo1),sampling_pattern); % img_echo1.*sampling_pattern;

if use_espirit_data
    n = num_coils;
else
    n = 2; % 2 x n espirit maps, we have 2 coils
end

%% Zero-filled reconstruction versus ESPIRiT

% Zero-filled reconstruction (i.e. direct rss from undersampled data)
% Do IFFT to get image data
us_coilimg = bart('fft -i 7',us_ksp_data_echo1);
% Combine coil data using root sum of squares
us_rss = bart('rss 4', us_coilimg);

% Add singleton dimension so first 3 dimensions reflect x, y, z
us_ksp_data_echo1_slice = permute(us_ksp_data_echo1,[1,2,4,3]);

% SENSE reconstruction
% Maps generated using autocalibration
sensemaps = bart('caldir 20', us_ksp_data_echo1_slice);
% 'PICS' - parallel-imaging compressed-sensing reconstruction
% Use sliced, undersampled k space data and generated maps to reconstruct
% image
sensereco = bart('pics -l1 -r0.01', us_ksp_data_echo1_slice, sensemaps);

% ESPIRiT reconstruction

% Maps generated using ESPIRiT calibration
% % Dimensions from espirit.m:
% % Input to this was 1x230x180x8 (breaks if input is 3D). Expects 3 k space
% % dimensions, 1 coil dimension
% % calib dims: 1x230x180x8x2 (two sets of maps)
% % emaps dims: 1x230x180x1x2
% % [calib emaps] = bart('ecalib -r 20', <image slice>);
% % r: cal size. k: kernel size. m: # maps.
espiritmaps = bart('ecalib -r 60 -k 5 -m 2', us_ksp_data_echo1_slice); %bart('ecalib -S', us_img_slice_dim4);
% size of input 1 = size of input 2
% l1-wavelet, l2 regularization
% r: regularization parameter
espiritreco = bart('pics -l1 -r0.01', us_ksp_data_echo1_slice, espiritmaps);
% espiritreco_rss dims in espirit.m were 1x320x252
espiritreco_rss = bart('rss 16', espiritreco);

%% Figures

fig = figure;
ax = gca;

imshow(abs(rss),[])
% For some reason, when you've used imshow and you don't put in the current
% axis, the title appears on the previous axis
title(ax,'Original image','FontSize',15)
set(fig,'Position',[50 600 300 250])

% Plot the original k space, sampling mask, and undersampled k space for
% previously sliced data
% If data were generated multiple times, the k space/calib maps from the
% last iteration will be plotted
fig2 = figure;
s1 = subplot(1,3,1);
imagesc(abs(squeeze(ksp_data_echo1(:,:,1,coil))))
title(sprintf('Original K-space for coil %d',coil),'FontSize',15)

s2 = subplot(1,3,2);
imagesc(sampling_pattern);
title('Sampling pattern','FontSize',15)


s3 = subplot(1,3,3);
imagesc(abs(us_ksp_data_echo1(:,:,1,coil)))
title(sprintf('Undersampled K-space for coil %d',coil),'FontSize',15)

s = suptitle(sprintf('Slice %d of z-axis',slice));
set(s,'FontSize',16,'FontWeight','bold')
set(fig2,'Position',[300 600 800 300])

% View ESPIRiT maps
fig3 = figure;
ax = gca;
imshow3(abs(squeeze(espiritmaps)),[],[2,n])
title(ax,'ESPIRiT maps','FontSize',15)
set(fig3,'Position',[50 100 400 400])

% Comparisons of different reconstruction methods
% Plot images from last iteration
% Zero-filled reconstruction
fig4 = figure;
subplot(2,2,1)
ax = gca;
imshow(us_rss,[])
title(ax,'Zero-filled reconstruction','FontSize',15);

subplot(2,2,2)
ax = gca;
espiritreco_map1 = abs(espiritreco(:,:,:,1,1));
imshow(espiritreco_map1, [])
title(ax,'ESPIRiT recon (map 1)','FontSize',15)

% This image should be the most accurate/closest to image from fully
% sampled k space
subplot(2,2,3)
ax = gca;
imshow(espiritreco_rss,[]);
title(ax,'ESPIRiT rss','FontSize',15)

subplot(2,2,4)
ax = gca;
imshow(abs(sensereco),[])
title(ax,'SENSE reconstruction','FontSize',15)
set(fig4,'Position',[450 100 400 400])

% Assess quality of reconstruction
% https://onlinelibrary.wiley.com/doi/full/10.1002/ima.22260
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5240324/
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4917755/
% Relative L2-norm error (RLNE), high frequency error norm (HFEN) (?),
% peak signal-to-noise ratio (PSNR), structural similarity (SSIM)
% For these metrics, normalize data to range from 0-255

% RLNE: ||x_hat-x||/||x||.
% Lower RLNE: reconstruction closer to fully sampled image

rss_s = permute(rss,[3 1 2]); % Add singleton dimension to match other images

% Fully sampled image
x = abs(squeeze(gs_normalize(rss_s,255)));
% Zero-filled reconstruction
x_hat_zero = abs(gs_normalize(us_rss,255));
% ESPIRiT reconstruction (2 combined maps)
x_hat_espirit = abs(gs_normalize(espiritreco_rss,255));
% SENSE reconstruction
x_hat_sense = abs(gs_normalize(sensereco,255));

rlne_zero = get_rlne(x,x_hat_zero);
rlne_espirit = get_rlne(x,x_hat_espirit);
rlne_sense = get_rlne(x,x_hat_sense);

% PSNR: 20*log10(255/sqrt(MSE)).
% Higher PSNR: reconstructed pixel value more consistent to original image
psnr_zero = get_psnr(x,x_hat_zero,255);
psnr_espirit = get_psnr(x,x_hat_espirit,255);
psnr_sense = get_psnr(x,x_hat_sense,255);

% Modified RLNE/PSNR metrics

% % Adjusting for difference in mean intensity (luminance):
% 
% % Method (1): normalize by mean intensity
% zero_diff = abs(mean(x_hat_zero(:))-mean(x(:)));
% espirit_diff = abs(mean(x_hat_espirit(:))-mean(x(:)));
% sense_diff = abs(mean(x_hat_sense(:))-mean(x(:)));
% 
% norm_rlne_zero = rlne_zero/zero_diff;
% norm_rlne_espirit = rlne_espirit/espirit_diff;
% norm_rlne_sense = rlne_sense/sense_diff;
% 
% norm_psnr_zero = psnr_zero*zero_diff;
% norm_psnr_espirit = psnr_espirit*espirit_diff;
% norm_psnr_sense = psnr_sense*sense_diff;
% 
% % Method (2): Make all images zero mean, unit variance
% % Fully sampled image
% x_zmuv = abs(zmuv(squeeze(rss_s)));
% % Zero-filled reconstruction
% x_hat_zero_zmuv = abs(zmuv(us_rss));
% % ESPIRiT reconstruction (2 combined maps)
% x_hat_espirit_zmuv = abs(zmuv(espiritreco_rss));
% % SENSE reconstruction
% x_hat_sense_zmuv = abs(zmuv(sensereco));
% 
% rlne_zero_zmuv = get_rlne(x_zmuv,x_hat_zero_zmuv);
% rlne_espirit_zmuv = get_rlne(x_zmuv,x_hat_espirit_zmuv);
% rlne_sense_zmuv = get_rlne(x_zmuv,x_hat_sense_zmuv);
% 
% psnr_zero_zmuv = get_psnr(x_zmuv,x_hat_zero_zmuv,1);
% psnr_espirit_zmuv = get_psnr(x_zmuv,x_hat_espirit_zmuv,1);
% psnr_sense_zmuv = get_psnr(x_zmuv,x_hat_sense_zmuv,1);

fprintf(sprintf('\nMean intensity diff:\n\tZero-filled recon: %.4f\n\tESPIRiT recon: %.4f\n\tSENSE recon: %.4f\n',...
    zero_diff,espirit_diff,sense_diff));

fprintf(sprintf('\nRLNE:\n\tZero-filled recon: %.4f\n\tESPIRiT recon: %.4f\n\tSENSE recon: %.4f\n',...
    rlne_zero,rlne_espirit,rlne_sense));
% fprintf(sprintf('Normalized RLNE:\n\tZero-filled recon: %.4f\n\tESPIRiT recon: %.4f\n\tSENSE recon: %.4f\n',...
%     norm_rlne_zero,norm_rlne_espirit,norm_rlne_sense));
% fprintf(sprintf('RLNE with ZMUV data:\n\tZero-filled recon: %.4f\n\tESPIRiT recon: %.4f\n\tSENSE recon: %.4f\n',...
%     rlne_zero_zmuv,rlne_espirit_zmuv,rlne_sense_zmuv));

fprintf(sprintf('\nPSNR:\n\tZero-filled recon: %.4f\n\tESPIRiT recon: %.4f\n\tSENSE recon: %.4f\n',...
    psnr_zero,psnr_espirit,psnr_sense));
% fprintf(sprintf('Normalized PSNR:\n\tZero-filled recon: %.4f\n\tESPIRiT recon: %.4f\n\tSENSE recon: %.4f\n',...
%     norm_psnr_zero,norm_psnr_espirit,norm_psnr_sense));
% fprintf(sprintf('PSNR with ZMUV data:\n\tZero-filled recon: %.4f\n\tESPIRiT recon: %.4f\n\tSENSE recon: %.4f\n',...
%     psnr_zero_zmuv,psnr_espirit_zmuv,psnr_sense_zmuv));


% Structural similarity (SSIM): equation 33 in
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4917755/
% SSIM: better image structures preserved

% Also assess using difference maps (?)
figure;ax = gca;imshow([abs(x-x_hat_zero),...
    abs(x-x_hat_espirit),...
    abs(x-x_hat_sense)],[]);
% For imshow, have to manually change current axes, otherwise the titles
% plot on the prev axes, and suptitle doesn't accept axes as an input
axes(ax);ax.FontSize = 13;
title(sprintf('Zero-filled recon\t\t\t\t\t\tESPIRiT recon (2 maps)\t\t\t\t\t\tSENSE recon'));
t = suptitle('Difference maps');set(t,'FontSize',16,'FontWeight','bold')

function y = gs_normalize(x,scale)
% Args:
% % x: grayscale images (dimensions X x Y)
% Returns:
% % y: normalized grayscale images
% Normalize value of x to 0-255 range.
    y = (x-min(x(:)))*scale/(max(x(:))-min(x(:)));
end

function rlne = get_rlne(x,x_hat)
% Args:
% % x: grayscale image constructed from fully sampled k space
% % x_hat: grayscale imageS constructed from undersampled k space
% % x_hat dimensions are X x Y
% Returns:
% % rlne: relative L2-norm error

rlne = norm(x_hat-x)./norm(x);

end

function psnr = get_psnr(x,x_hat,scale)
% Args:
% % x: grayscale image constructed from fully sampled k space
% % x_hat: grayscale imageS constructed from undersampled k space
% % x_hat dimensions are X x Y
% Returns:
% % psnr: peak signal-to-noise ratio

psnr = 20*log10(scale/sqrt(immse(x,x_hat)));

end

function normalized_x = zmuv(x)
% Args:
% % x: grayscale images (dimensions X x Y)
% Returns:
% % normalized_x: x normalized so each image has zero mean, unit variance

normalized_x = reshape(zscore(x(:)),size(x,1),size(x,2));

end