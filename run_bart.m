lose all;clear;clc;

% If you don't want to run RussRecon (and install all the folders), set
% load_data to false, download 'bruker_data.mat' and comment out all
% 'addpath' lines except addpath(bartpath)

% This script loads in k space data and displays a slice of the image data.
% It then uses a sampling mask to generate undersampled k space data and
% displays three different types of reconstructions: zero-filled recon
% (i.e. directly reconstructing the data), SENSE reconstruction
% (autocalibrated), and ESPIRiT reconstruction (using ESPIRiT calibration)

% BART must be installed for bart commands to work
% %! means these settings must be changed for each user depending on where
% your code/BART path is
% %* means you can toggle these options

main_dir = '/Users/alex/janet/Toolboxes'; %!

%bartpath = '/Users/alex/janet/Toolboxes/bart'; %!
bartpath = '/Users/alex/alex_code/bart-0.5.00'; % Compiled path for kefalonia


%bartmatlabpath = '/Users/alex/janet/Toolboxes/bart/matlab'; %!*
bartmatlabpath= [bartpath '/matlab'];

[addecode_dir,~,~] = fileparts(matlab.desktop.editor.getActiveFilename)
setenv('TOOLBOX_PATH',bartpath)

% The following packages are needed in your main directory:
paths = {main_dir;[main_dir '/STI_Suite_v2/STI_Suite_v2.1/Source Code_v2.1'];...
    [main_dir '/RussRecon'];[main_dir '/ImageProcessing'];...
    [main_dir '/NIfTI_20140122'];[main_dir '/MultipleEchoRecon'];...
    [main_dir '/XCalc'];addecode_dir};


 recon_along_x = 1; % 14 February 2020--our protocols will recon along x instead of z.

% If you need to extract the image data, you need to have the .work folder
% in the main directory. If you're loading it in (and bruker_img.MAT should
% be in the addecode main directory), you don't need this folder
workpath = [main_dir '/B04027.work'];

use_espirit_data = false; %* try out the example image or use our image?
% If you want to use the example image, the BART folder at
% https://github.com/mikgroup/espirit-matlab-examples must be downloaded
% (or just the 'full.cfl' from the 'data' folder). Put the correct path
% name here (don't include the '.cfl' at the end):
fullpath = 'BART/espirit-matlab-examples-master/data/full'; %!

load_data = false; %* 'true' if you have the k space data saved and don't
% want to run get_RussRecon_img again
save_data = false; %* 'true' if you want to save the data
data_path = 'bruker_data.mat'; % image path if it's already been saved

bart_mask = true; %* if generating mask, use variable density Poisson mask
% from BART code (true), or code from BJ Anderson (false)? Right now, the
% undersampling factor can be better controlled by the Anderson code
load_mask = false; %* load sampling mask in? Currently, only BJ mask is saved
num_iter = 1; % Number of sampling patterns/reconstructions to generate
% (to get mean RLNE, PSNR values)

acceleration = 4; %* Undersample by a factor of <acceleration>

% !!! ^ currently implementing
show_recon_steps = true; %* Show sampling pattern, calibration maps, etc.?
% Width and height of images in plot comparing reconstruction quality
wh = [0.2582,0.4128];
gen_recon_plots = true; %* Generate plot comparing reconstruction quality?
%* Save plot comparing reconstruction quality? Only saves plots if they were actually generated.
save_recon_plots = false;
gen_diff_maps = true; %* Generate difference maps

% Directory in which to save image slices
plot_dirpath = '/Users/janetchen/Documents/Bass Connections/Reconstructed images';
%* Save 3D images (fully-sampled and undersampled-ESPIRiT recon images) to NIfTI file
save_nifti = false;

addpath(bartmatlabpath)
addpath(main_dir)
%if ~use_espirit_data && ~load_data % Need to be in workpath for RussRecon
    %cd(workpath)
    % Add paths so that RussRecon can be used
    for ii = 1:length(paths)
        addpath(paths{ii,1})
    end
%end

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

if use_espirit_data
    num_coils = size(ksp_data,4);
    sz_x = size(ksp_data,2); sz_y = size(ksp_data,3);
else
    num_echoes = size(ksp_data,4); % Will need this if reconning all echoes.
    num_coils = size(ksp_data,5);
    % Z-axis slice to look at
    % 8 echoes, 2 coils
    % Original k space data dimensions: x, y, z, # echoes, # coils
    % ksp_data matrix is 192x192x90x8x2
    
    if recon_along_x
        sz_x = size(ksp_data,2); sz_y = size(ksp_data,3);
    else
        sz_x = size(ksp_data,1); sz_y = size(ksp_data,2);
    end
end

%% Generate sampling pattern and undersample K-space

% Coil # for which to plot k space data and undersampling
coil = 1;

if bart_mask
    % bart poisson -Y $dim -Z $dim -y $yaccel -z $zaccel -C $caldim -v -e mask
    % -Y: dim 1. -Z: dim 2. -v: variable density (which the workshop uses)
    % -e: elliptical scanning
    % Default values: Y: 192, Z: 192, y: 1.5, z: 1.5, C: 32
    % ~25% of elements in und2x2 in espirit.m were nonzero
    % !!! Although this should reflect 9216 points, only 8005 points are
    % generated
    num_points = sz_x*sz_y/acceleration;
    
    % Sampling pattern. Binary mask (ones and zeros), where ones are the
    % only points that we sample
    % Original output of BART command is 3 dimensional, with first
    % dimension of length 1, so squeeze to remove this dimension
    % C manually tuned to 60 for best reconstruction results
    sampling_pattern = squeeze(bart(sprintf('poisson -Y %d -Z %d -y 1 -z 1 -C 60 -R %d -v -e',sz_x,sz_y,num_points)));
else
    % Mask from BJ Anderson's code
    if load_mask
        load('bj_sampling_pattern.mat')
    else
        % Inputs
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

ksp_dims = size(ksp_data);
% Collect original and reconstructed images in 3D arrays
fs_3d_img = zeros(ksp_dims(1:3));
espirit_3d_img = zeros(ksp_dims(1:3));

%% Create simulated data from arbitrary images.
% Most of this work belongs outside of the slice loop.
% We may have multi-echo data, but for now let's just work on echo 1.
% We want to make echo_num a variable so we can generalize later.
echo_num=1;

if ~use_espirit_data
    % To match espirit 'full' img dimensions (XxYxcoils, 230x180x8), slice
    % in z-dimension IN IMAGE domain
    % Inverse Fourier transform the k space data to get the image data
    % 'squeeze' just removes singleton dimensions (i.e. dimensions of
    % length 1)
    
    ksp_echo_n_data=squeeze(ksp_data(:,:,:,echo_num,:));
    
    % Only going to transform readout direction
    if recon_along_x
        readout_dim=1;
    else
        readout_dim=3;
    end
    
    %ksp_echo_n_ifft = bart('fft -i 7',squeeze(ksp_data(:,:,:,1,:)));
    
	ksp_echo_n_ifft = fftshift(ifft(fftshift(ksp_echo_n_data,readout_dim),[],readout_dim),readout_dim); 
end

%% Iterate through z-slices

% To run through all slices, set below to 1:size(ksp_data,3)
slices_to_generate = 52; %* Reconstruct images for these slices
%* Z-axis slice to keep for reconstruction quality comparison (say you're
% generating images for all slices but want to compare reconstruction
% quality for one slice). Shows difference maps, quality metrics
slices_to_compare = 52;

%% Iterate through x-slices
% 14 February 2020: We are resetting the variables directly above, as we
% want to iterate through the x-slices.

% To run through all slices, set below to 1:size(ksp_data,1) <--note
% difference compared to code above (dim 1 vs dim 3)
slices_to_generate = 52; %* Reconstruct images for these slices
%* X-axis slice to keep for reconstruction quality comparison (say you're
% generating images for all slices but want to compare reconstruction
% quality for one slice). Shows difference maps, quality metrics
slices_to_compare = 52;


% !!! Slices 68-72 have issues
for slice = slices_to_generate
    if use_espirit_data
        % img is 1x230x180x8
        % Swap x, y, z dimensions so 3rd dimension is 1
        ksp_data_echo_n_z1 = permute(ksp_data(1,:,:,1:num_coils),[2,3,1,4]);
        % Take inverse Fourier transform
        coilimg = bart('fft -i 6', ksp_data);
        % Take the root sum of squares to produce one image from multiple coils
        fs_finalimg = bart('rss 4',squeeze(coilimg));
        slice = 1; % Original X dimension is 1 in length, so use that 'slice'
    else
  
        % Take the z-axis slice
        % coilimg = ksp_echo_n_ifft(:,:,slice,:);
 
        if (recon_along_x)
            ksp_data_echo_n_z1 = squeeze(ksp_echo_n_ifft(slice,:,:,:));
        else
            ksp_data_echo_n_z1 = squeeze(ksp_echo_n_ifft(:,:,slice,:));
        end
        
        coilimg = bart('fft -i 7',ksp_data_echo_n_z1);

        % Take the root sum of squares to produce one image from multiple coils
        fs_finalimg = bart('rss 4', squeeze(coilimg));
    end
    
    % Zero out unsampled areas through element-by-element multiplication with
    % sampling pattern
    % 'us' means 'undersampled'
    us_ksp_echo_n_z1 = bart('fmac',squeeze(ksp_data_echo_n_z1),sampling_pattern);
    % Equivalently: squeeze(ksp_data_echo_n_z1).*sampling_pattern;
    
    %% Zero-filled reconstruction versus ESPIRiT
    
    % Zero-filled reconstruction (i.e. direct rss from undersampled data)
    % Do IFFT to get image data
    us_coilimg = bart('fft -i 7',us_ksp_echo_n_z1);
    % Combine coil data using root sum of squares. Reconstruction from under-
    % sampled data without use of calibration/PICS is known as zero-filled
    % econstruction
    zerofilled_finalimg = bart('rss 4', us_coilimg);
    
    % Add singleton dimension so first 3 dimensions reflect x, y, z
    us_ksp_data_echo_n_slice = permute(us_ksp_echo_n_z1,[1,2,4,3]);
    
    % SENSE reconstruction (direct calibration from k-space center)
    % Maps generated using autocalibration
    % !!!Increasing caldir size seems to improve recon? Attempted up to 60
    sense_maps = bart('caldir 20', us_ksp_data_echo_n_slice);
    % 'PICS' - parallel-imaging compressed-sensing reconstruction
    % Use sliced, undersampled k space data and generated maps to reconstruct
    % image
    % r: regularization parameter, set to 0.001 in Uecker 2014 paper:
    % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4142121/
    sense_finalimg = bart('pics -l1 -r0.01', us_ksp_data_echo_n_slice, sense_maps);
    
    % ESPIRiT reconstruction
    
    % Maps generated using ESPIRiT calibration
    % % Dimensions from espirit.m:
    % % Input to this was 1x230x180x8 (breaks if input is 3D). Expects 3 k space
    % % dimensions, 1 coil dimension
    % % calib dims: 1x230x180x8x2 (two sets of maps)
    % % emaps dims: 1x230x180x1x2
    % % r: cal size. k: kernel size. m: # maps.
    % Don't increase kernel size past 5 - details start to disappear
    espirit_maps = bart('ecalib -r 20 -k 5 -m 2', us_ksp_data_echo_n_slice); %bart('ecalib -S', us_img_slice_dim4);
    % size of input 1 = size of input 2
    % l1-wavelet, l2 regularization
    % r: regularization parameter (0.001 in Uecker 2014 paper:
    % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4142121/)
    % 'pics -l1 -r0.01
    espirit_coilimg = squeeze(bart('pics -l1 -r0.01', us_ksp_data_echo_n_slice, espirit_maps));
    espirit_finalimg = bart('rss 4', espirit_coilimg); % 16 if 'squeeze' not used above
    
    %% Figures
    
    if show_recon_steps
        fig = figure;
        ax = gca;
        
        imshow(abs(fs_finalimg),[])
        % For some reason, when you've used imshow and you don't put in the current
        % axis, the title appears on the previous axis
        title(ax,'(a) Original image','FontSize',15)
        set(fig,'Position',[50 600 300 250])
        
        % Plot the original k space, sampling mask, and undersampled k space for
        % previously sliced data
        % If data were generated multiple times, the k space/calib maps from the
        % last iteration will be plotted
        fig2 = figure;
        s1 = subplot(1,3,1);

        imagesc(abs(squeeze(ksp_data_echo_n_z1(:,:,coil))))

        title(sprintf('(a) Original K-space (coil %d)',coil),'FontSize',15)
        
        s2 = subplot(1,3,2);
        imagesc(sampling_pattern);
        title('(b) Sampling pattern','FontSize',15)
        
        
        s3 = subplot(1,3,3);

        title(sprintf('(c) Undersampled K-space (coil %d)',coil),'FontSize',15)
        imagesc(squeeze(abs(us_ksp_echo_n_z1(:,:,1,coil))))
        if (recon_along_x) 
            s = suptitle(sprintf('Slice %d of x-axis',slice));  
        else  
            s = suptitle(sprintf('Slice %d of z-axis',slice));
        end
        
        
        set(s,'FontSize',16,'FontWeight','bold')
        set(fig2,'Position',[300 600 800 300])
        
        % View ESPIRiT maps
        fig3 = figure;
        ax = gca;
        
        % imshow3 seems to be breaking right now, even though it is in the
        % MATLAB path...commenting out for now.
        % imshow3(abs(squeeze(espirit_maps)),[],[2,num_coils])
        title(ax,'ESPIRiT maps','FontSize',15)
        set(fig3,'Position',[50 100 400 400])
    end
    
    if gen_recon_plots
        % Comparisons of different reconstruction methods
        % Plot images from last iteration
        % Zero-filled reconstruction
        fig4 = figure;
        subplot(2,3,[1 4])
        fs_ax = gca;
        imagesc(fs_finalimg);colormap('gray')
        title('(a) Original image','FontSize',15);
        
        subplot(2,3,2)
        ax = gca;
        imshow(zerofilled_finalimg,[])
        title(ax,'(b) Zero-filled recon','FontSize',15);
        posn(1,:) = get(ax,'Position');
        
        subplot(2,3,3)
        ax(2) = gca;
        espirit_map1_finalimg = abs(espirit_coilimg(:,:,1));
        imshow(espirit_map1_finalimg, [])
        title(ax(2),'(c) ESPIRiT recon (map 1)','FontSize',15)
        posn(2,:) = get(ax(2),'Position');
        
        % This image should be the most accurate/closest to image from fully
        % sampled k space
        subplot(2,3,5)
        ax(3) = gca;
        imshow(espirit_finalimg,[]);
        title(ax(3),'(d) ESPIRiT rss','FontSize',15)
        posn(3,:) = get(ax(3),'Position');
        
        subplot(2,3,6)
        ax(4) = gca;
        imshow(abs(sense_finalimg),[])
        title(ax(4),'(e) SENSE recon','FontSize',15)
        set(fig4,'Position',[400 100 500 300])
        posn(4,:) = get(ax(4),'Position');
        
        % Reposition image axes
        posn([1,3],1) = 0.34;posn([2,4],1) = 0.68;
        posn([1,2],2) = posn([1,2],2)-0.07;posn([3,4],2) = posn([3,4],2)-0.1;
        posn(:,[3,4]) = repmat(wh,4,1);
        fs_posn = get(fs_ax,'Position'); fs_posn=[0.05 0.33 wh];
        set(fs_ax,'Position',fs_posn)
        for ii = 1:size(posn,1)
            set(ax(ii),'Position',posn(ii,:))
        end
        
        if save_recon_plots
            %set(gcf,'Color','none')
            saveas(fig4,sprintf('%s/Slice %d.png',plot_dirpath,slice))
        end
    end
    
   
    if (recon_along_x)
        % Add current x-axis slice to 3D arrays
        fs_3d_img(slice,:,:) = fs_finalimg;
        espirit_3d_img(slice,:,:) = espirit_finalimg;         
    else
        % Add current z-axis slice to 3D arrays
        fs_3d_img(:,:,slice) = fs_finalimg;
        espirit_3d_img(:,:,slice) = espirit_finalimg;
    end
    
    
    
    %close all;%close(fig4)
    
    if gen_diff_maps && ismember(slice,slices_to_compare)
        
        fs_rss = permute(fs_finalimg,[3 1 2]); % Add singleton dimension to match other images
        
        % Fully sampled image
        x = abs(squeeze(gs_normalize(fs_rss,255)));
        % Zero-filled reconstruction
        x_hat_zero = abs(gs_normalize(zerofilled_finalimg,255));
        % ESPIRiT reconstruction (2 combined maps)
        x_hat_espirit = abs(gs_normalize(espirit_finalimg,255));
        % SENSE reconstruction
        x_hat_sense = abs(gs_normalize(sense_finalimg,255));
        % Also assess using difference maps (?)
        figure;ax = gca;imshow([abs(x-x_hat_zero),...
            abs(x-x_hat_espirit),...
            abs(x-x_hat_sense)],[]);
        % For imshow, have to manually change current axes, otherwise the titles
        % plot on the prev axes, and suptitle doesn't accept axes as an input
        axes(ax);ax.FontSize = 13;
        title(sprintf('(a) Zero-filled recon\t\t\t\t\t(b) ESPIRiT recon (2 maps)\t\t\t\t\t\t(c) SENSE recon'));
        t = suptitle('Difference maps');set(t,'FontSize',16,'FontWeight','bold')
    end
    fprintf(sprintf('Reconstructed slice %d\n',slice))
end

if save_nifti
    fs_nii = make_nii(fs_3d_img,[],ksp_dims(1:3),64);
    save_nii(fs_nii,sprintf('%s/fs_img.nii',main_dir));
    
    espirit_nii = make_nii(espirit_3d_img, [], ksp_dims(1:3), 64);
    save_nii(espirit_nii,sprintf('%s/espirit_recon.nii',main_dir));
end

% Assess quality of reconstruction for a slice
% https://onlinelibrary.wiley.com/doi/full/10.1002/ima.22260
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5240324/
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4917755/
% Relative L2-norm error (RLNE), high frequency error norm (HFEN) (?),
% peak signal-to-noise ratio (PSNR), structural similarity (SSIM)
% For these metrics, normalize data to range from 0-255

fs_rss = permute(fs_finalimg,[3 1 2]); % Add singleton dimension to match other images

% !!! Update code to calculate these even if this statement isn't
% entered, because RLNE/PSNR is calculated outside of this for loop
% Fully sampled image
x = abs(squeeze(gs_normalize(fs_rss,255)));
% Zero-filled reconstruction
x_hat_zero = abs(gs_normalize(zerofilled_finalimg,255));
% ESPIRiT reconstruction (2 combined maps)
x_hat_espirit = abs(gs_normalize(espirit_finalimg,255));
% SENSE reconstruction
x_hat_sense = abs(gs_normalize(sense_finalimg,255));

% RLNE: ||x_hat-x||/||x||.
% Lower RLNE: reconstruction closer to fully sampled image

rlne_zero = get_rlne(x,x_hat_zero);
rlne_espirit = get_rlne(x,x_hat_espirit);
rlne_sense = get_rlne(x,x_hat_sense);

% PSNR: 20*log10(255/sqrt(MSE)).
% Higher PSNR: reconstructed pixel value more consistent to original image
psnr_zero = get_psnr(x,x_hat_zero,255);
psnr_espirit = get_psnr(x,x_hat_espirit,255);
psnr_sense = get_psnr(x,x_hat_sense,255);

% % Modified RLNE/PSNR metrics
%
% % % Adjusting for difference in mean intensity (luminance):
% %
% % % Method (1): normalize by mean intensity
% % zero_diff = abs(mean(x_hat_zero(:))-mean(x(:)));
% % espirit_diff = abs(mean(x_hat_espirit(:))-mean(x(:)));
% % sense_diff = abs(mean(x_hat_sense(:))-mean(x(:)));
% %
% % norm_rlne_zero = rlne_zero/zero_diff;
% % norm_rlne_espirit = rlne_espirit/espirit_diff;
% % norm_rlne_sense = rlne_sense/sense_diff;
% %
% % norm_psnr_zero = psnr_zero*zero_diff;
% % norm_psnr_espirit = psnr_espirit*espirit_diff;
% % norm_psnr_sense = psnr_sense*sense_diff;
% %
% % % Method (2): Make all images zero mean, unit variance
% % % Fully sampled image
% % x_zmuv = abs(zmuv(squeeze(rss_s)));
% % % Zero-filled reconstruction
% % x_hat_zero_zmuv = abs(zmuv(us_rss));
% % % ESPIRiT reconstruction (2 combined maps)
% % x_hat_espirit_zmuv = abs(zmuv(espiritreco_rss));
% % % SENSE reconstruction
% % x_hat_sense_zmuv = abs(zmuv(sensereco));
% %
% % rlne_zero_zmuv = get_rlne(x_zmuv,x_hat_zero_zmuv);
% % rlne_espirit_zmuv = get_rlne(x_zmuv,x_hat_espirit_zmuv);
% % rlne_sense_zmuv = get_rlne(x_zmuv,x_hat_sense_zmuv);
% %
% % psnr_zero_zmuv = get_psnr(x_zmuv,x_hat_zero_zmuv,1);
% % psnr_espirit_zmuv = get_psnr(x_zmuv,x_hat_espirit_zmuv,1);
% % psnr_sense_zmuv = get_psnr(x_zmuv,x_hat_sense_zmuv,1);

% fprintf(sprintf('\nMean intensity diff:\n\tZero-filled recon: %.4f\n\tESPIRiT recon: %.4f\n\tSENSE recon: %.4f\n',...
%     zero_diff,espirit_diff,sense_diff));
%
fprintf(sprintf('\nRLNE:\n\tZero-filled recon: %.4f\n\tESPIRiT recon: %.4f\n\tSENSE recon: %.4f\n',...
    rlne_zero,rlne_espirit,rlne_sense));
% % fprintf(sprintf('Normalized RLNE:\n\tZero-filled recon: %.4f\n\tESPIRiT recon: %.4f\n\tSENSE recon: %.4f\n',...
% %     norm_rlne_zero,norm_rlne_espirit,norm_rlne_sense));
% % fprintf(sprintf('RLNE with ZMUV data:\n\tZero-filled recon: %.4f\n\tESPIRiT recon: %.4f\n\tSENSE recon: %.4f\n',...
% %     rlne_zero_zmuv,rlne_espirit_zmuv,rlne_sense_zmuv));
%
% fprintf(sprintf('\nPSNR:\n\tZero-filled recon: %.4f\n\tESPIRiT recon: %.4f\n\tSENSE recon: %.4f\n',...
%     psnr_zero,psnr_espirit,psnr_sense));
% % fprintf(sprintf('Normalized PSNR:\n\tZero-filled recon: %.4f\n\tESPIRiT recon: %.4f\n\tSENSE recon: %.4f\n',...
% %     norm_psnr_zero,norm_psnr_espirit,norm_psnr_sense));
% % fprintf(sprintf('PSNR with ZMUV data:\n\tZero-filled recon: %.4f\n\tESPIRiT recon: %.4f\n\tSENSE recon: %.4f\n',...
% %     psnr_zero_zmuv,psnr_espirit_zmuv,psnr_sense_zmuv));


% Structural similarity (SSIM): equation 33 in
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4917755/
% SSIM: better image structures preserved

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