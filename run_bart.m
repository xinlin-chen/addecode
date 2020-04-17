close all;clear;clc;
% Xinlin Chen
% version dated 2020/04/17

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
% your code/data/toolbox/BART path is, or where you want to save figures
% %* means you can toggle these options

main_dir = '/Users/janetchen/Documents/Bass Connections';
data_dir = '/Users/janetchen/Documents/Bass Connections';
bartpath = '/Applications/bart';

%% Set, add directories

% main_dir = '/Users/alex/janet/Toolboxes'; %!
% data_dir = '/Users/alex/janet/Data'; %!
% bartpath = '/Users/alex/alex_code/bart-0.5.00'; %! Compiled path for kefalonia

bartmatlabpath= [bartpath '/matlab'];

[addecode_dir,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
setenv('TOOLBOX_PATH',bartpath)

% The following packages are needed in your main directory:
paths = {main_dir;[main_dir '/STI_Suite_v2/STI_Suite_v2.1/Source Code_v2.1'];...
    [main_dir '/RussRecon'];[main_dir '/ImageProcessing'];...
    [main_dir '/NIfTI_20140122'];[main_dir '/MultipleEchoRecon'];...
    [main_dir '/XCalc'];addecode_dir};

recon_along_x = true; %!! 14 February 2020--our protocols will recon along x instead of z.

% If you need to extract the image data, you need to have the .work folder
% in the data directory. If you're loading it in (and bruker_img.MAT should
% be in the addecode main directory), you don't need this folder
%workpath = [data_dir '/B04027.work'];
data_strs = {'/B04027.work','/127','/115'}; %! Ensure you have correct image paths
data_ind = 1; % Pick the image
data_str = data_strs{data_ind};
workpath = [data_dir data_str];
%! Re-center k-space. Set to false for B04027, which is already centered
if strcmp(data_str,'/B04027.work')
    do_fft_shift = false; % This FFT shift re-centers k-space
else
    do_fft_shift = true;
end
y_shift = 10;

addpath(bartmatlabpath)
addpath(main_dir)
% Add paths so that RussRecon can be used
for ii = 1:length(paths)
    addpath(paths{ii,1})
end

%% General settings
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

% Whether to generate quality metrics (QM: RLNE/PSNR/SSIM) for each slice
gen_quality_metrics = true;

% Whether to plot QM (gen_quality_metrics should be set to true if this is)
plot_quality_metrics = false;
% When saving QM figure generated with BJ's sampling pattern, note this in
% the filename
if ~bart_mask && gen_quality_metrics
    mask = 'bj_';
else
    mask = '';
end

%% Plot settings

ksp_bwh = [0.08 0.29 0.8]; %* Bottom position and size of plots of K-space
%* Width and height of images in plot comparing reconstruction quality
recon_wh = [0.325,0.35];
poster_wh = [0.9,0.32]; %0.28];
cols = 'brgmc';
show_recon_steps = false; %* Show sampling pattern, calibration maps, etc.?
plot_recon_figs = false; %* Generate plot comparing reconstruction quality?
%* Save plot comparing reconstruction quality? Only saves plots if they were actually generated.
save_recon_plots = false;
plot_poster_figs = false; % Generate plots for Neuro poster session?
plot_diff_maps = false; %* Generate difference maps

% Directory in which to save image slices
recon_img_dirpath = '/Users/janetchen/Documents/Bass Connections/Reconstructed images'; %!
%* Save 3D images (fully-sampled and undersampled-ESPIRiT recon images) to NIfTI file
save_nifti = false;

% Directory in which to save other figures
fig_dirpath = '/Users/janetchen/Dropbox/Bass Connections/addecode/Figures'; %!

%% Load or read in data

if use_espirit_data
    % Read in k space data from BART's example (fully sampled)
    ksp_data = readcfl(sprintf('%s/%s',main_dir,fullpath));
else
    if load_data
        % Load in previously-extracted k space data
        load(sprintf('%s/%s',main_dir,data_path))
    else
        cd(workpath);
        % Get k space data from .work directory
        ksp_data=get_RussRecon_img('bruker','center','save');
        % kspace may be centered in along some axes but not others.
        
        if do_fft_shift % Don't recenter if using B04027.work
            % Recenter RARE k-space in z-direction
            ksp_data=fftshift(ksp_data, 3);
            % MANUALLY recenter in y-direction
            ksp_data = circshift(ksp_data,y_shift,2);
        end
        
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

ksp_dims = size(ksp_data);

%% Reconstruction settings

% To run through all slices, set below to 1:size(ksp_data,<dimension of interest>)
if recon_along_x
    slices_to_generate = 1:ksp_dims(1); % 57 %* Reconstruct images for these slices
    xlbl = 'X-axis slice';
else
    slices_to_generate = 1:ksp_dims(3);
    xlbl = 'Z-axis slice';
end
% !!! Slices 68-72 have issues
% slices_to_generate = 57;
%* Z-axis slice to keep for reconstruction quality comparison (say you're
% generating images for all slices but want to compare reconstruction
% quality for one slice). Shows difference maps, quality metrics
slices_to_compare = slices_to_generate;

%% Undersampling/recon settings

% Undersampling parameters
accel = 4; %* Undersample by a factor of <acceleration>. NOTE: this
% is often not the 'true' acceleration for some reason, so look at
% 'actual_accel' too (calculated directly from generated sampling pattern)
% For images 127 and B04027, setting accel to 4 with a calibration region
% size (cal_reg) of 40 results in a true acceleration close to 8 (~7.94 and
% 7.81, respectively)
cal_reg = 40; % Size of calibration region. 40. Previously 60

% Reconstruction parameters
reg_stepsize = 1e-2; %* Size of regularization parameter. 0.01
reg_method = 'l1'; %* Regularization method ('l1' or 'l2'). l1
% caldir.c, ecalib.c
recon_cal_size = cal_reg; %* Upper limit of calibration region size
if strcmp(data_str,'/B04027.work')
    espirit_num_maps = 2; %* Number of sensitivity maps to calculate.
    espirit_kernel_size = 4; %* Kernel size
else
    espirit_num_maps = 4;
    espirit_kernel_size = 5; %* Kernel size
end
% ^ Cannot be more than the number of coils used to acquire the image
gen_sense_recon = false; % Generate SENSE reconstruction (as well as ESPIRiT)?

%% Grid search
% Choose variable to iterate over
iter_vars = {'cal_reg','accel','reg_stepsize',...
    'espirit_kernel_size','espirit_num_maps'};
% Set iteration type from options in cell above
iter_type = 0; % Set to 0 to not iterate

% !!! Note: 'accel' is the input to the function generating the
% sampling pattern. For BART, the 'true' acceleration will depend on the
% size of the calibration region ('cal_reg') and the value of 'accel'
% e.g. for image 127, when cal_reg = 50 and acceleration = 1, the 'true'
% acceleration is 4.09
iter_var_titles = {'Calibration Region Size','Acceleration','Step-size',...
    'Kernel Size','Number of maps'};
switch iter_type
    case 0 % Don't iterate through anything
        iter_values = NaN;
    case 1 % Calibration region size
        iter_values = 20:10:50;
        if strcmp(data_str,'/127')
            % Image 1 (127):
            % For true acceleration of ~5.1:
            iter_values2 = [1,1.2,1.66,2.8];
            % For true acceleration of ~6.48:
            % iter_values2 = [1.58,1.9,2.67,7];
            % For true acceleration of ~7.95 (only for cal reg 20-40)
            %iter_values2 = [2.03,2.51,4];
            % Image 2 (B04027)
        else
            iter_values2 = ones(1,length(iter_values));
        end
    case 2 % Acceleration
        % These are nominal acceleration values
        if strcmp(data_str,'/127')
            iter_values = [1,1.8,2.4,3,3.8,4.4]; % for CRS 40
        else % Naive acceleration values
            iter_values = 1:5;
        end
    case 3 % Regularization step-size
        % % fewer stepsizes for L2 regularization
        % iter_values = [5e-3,1e-2,5e-2,1e-1];
        % % higher stepsizes
        %iter_values = [2.5e-2,5e-2,7.5e-2,1e-1]; %stepsize (slicing in x) -
        % % lower stepsizes
        iter_values = [5e-3,7.5e-3,1e-2,2.5e-2,5e-2]; %stepsize (slicing in x)
        % iter_values = [1e-4,1e-3,5e-3,1e-2,5e-2,1e-1]; %stepsize (slicing in z)
    case 4 % Kernel size
        iter_values = 2:5; % k size
    case 5 % Number of ESPiRIT maps to calculate
        % Max # maps is # coils (2 for image B04027, 4 for images 127, 115)
        iter_values = 1:espirit_num_maps;
    
end

% If acceleration or calibration region is being changed, the sampling
% pattern changes and therefore zero-filled reconstruction will also be
% affected (meaning more than one zero-filled plot is needed per quality
% metric for visualization purposes)
if iter_type > 0 && (strcmp(iter_vars{iter_type},'accel') || ...
        strcmp(iter_vars{iter_type},'cal_reg'))
    regenerate_mask = true;
else
    % If something like regularization step-size is being changed,
    % zero-filled reconstruction is not affected, and therefore you only
    % need to compare one ZF recon with all the ESPIRiT/SENSE recons
    regenerate_mask = false;
end

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

%% Iterate through desired recon settings (only allows you to search
% one parameter at a time
for iter = 1:length(iter_values)
    if iter_type > 0
        % Step through values for iterating variable
        eval(sprintf('%s = %d;',iter_vars{iter_type},iter_values(iter)))
        if strcmp(iter_vars{iter_type},'cal_reg')
            accel = iter_values2(iter);
        end
    end
    %% Generate sampling pattern and undersample K-space
    if iter == 1 || regenerate_mask % Re-generate if acceleration or calibration region is changed
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
            num_points = sz_x*sz_y/accel;
            
            % Sampling pattern. Binary mask (ones and zeros), where ones are the
            % only points that we sample
            % Original output of BART command is 3 dimensional, with first
            % dimension of length 1, so squeeze to remove this dimension
            % C manually tuned to 60 for best reconstruction results
            % Method 1: acceleration set by y, z
            sampling_pattern = squeeze(bart(sprintf('poisson -Y %d -Z %d -y %d -z %d -C %d -v -e',sz_x,sz_y,accel,accel,cal_reg)));
            % Method 2: acceleration set by R
%             sampling_pattern = squeeze(bart(sprintf('poisson -Y %d -Z %d -y 1 -z 1 -C %d -R %d -v -e',sz_x,sz_y,cal_reg,num_points)));
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
                sampling_pattern = sampling_mask(accel,sz_x,sz_y,pa,pb);
            end
        end
        
        % Show how much of the data is kept. Even though acceleration can
        % be set, it doesn't necessarily reflect the true acceleration
        % Also, the size of the calibration region affects the true
        % ('effective') acceleration
        actual_accel(iter) = 1/(length(find(sampling_pattern ~= 0))/numel(sampling_pattern));
        fprintf(sprintf('\tMask non-zero percentage: %.2f%%',100/actual_accel(iter)))
        fprintf('\n')
    end
    
    % !!! Now, peg ESPIRiT calibration region size to sampling region size
    % if it changes
    recon_cal_size = cal_reg;
    
    % Collect original and reconstructed images in 3D arrays
    fs_3d_img = zeros(ksp_dims(1:3)); % Fully sampled
    zf_3d_img = zeros(ksp_dims(1:3)); % Undersampled, zero-filled
    espirit_3d_img = zeros(ksp_dims(1:3)); % Recon using espirit
    if gen_sense_recon
        sense_3d_img = zeros(ksp_dims(1:3));
    end
    
    time_start(iter) = tic;
    %% Reconstruct every slice
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
            % Take the x/z-axis slice
            % coilimg = ksp_echo_n_ifft(:,:,slice,:);
            
            if recon_along_x
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
        zf_finalimg = bart('rss 4', us_coilimg); % Zero-filled final img
        
        % Add singleton dimension so first 3 dimensions reflect x, y, z
        if recon_along_x
            us_ksp_data_echo_n_slice = permute(us_ksp_echo_n_z1,[4,1,2,3]);
        else
            us_ksp_data_echo_n_slice = permute(us_ksp_echo_n_z1,[1,2,4,3]);
        end
        
        if gen_sense_recon || plot_recon_figs
            % SENSE reconstruction (direct calibration from k-space center)
            % Maps generated using autocalibration
            % !!!Increasing caldir size seems to improve recon? Attempted up to 60
            sense_maps = bart(sprintf('caldir %d',recon_cal_size), us_ksp_data_echo_n_slice);
            % 'PICS' - parallel-imaging compressed-sensing reconstruction
            % Use sliced, undersampled k space data and generated maps to reconstruct
            % image
            % r: regularization parameter, set to 0.001 in Uecker 2014 paper:
            % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4142121/
            sense_finalimg = squeeze(bart(sprintf('pics -%s -r%d',reg_method,reg_stepsize), us_ksp_data_echo_n_slice, sense_maps));
        end
        
        % ESPIRiT reconstruction
        
        % Maps generated using ESPIRiT calibration
        % % Dimensions from espirit.m:
        % % Input to this was 1x230x180x8 (breaks if input is 3D). Expects 3 k space
        % % dimensions, 1 coil dimension
        % % calib dims: 1x230x180x8x2 (two sets of maps)
        % % emaps dims: 1x230x180x1x2
        % % r: cal size. k: kernel size. m: # maps.
        % Don't increase kernel size past 5 - details start to disappear
        espirit_maps = bart(sprintf('ecalib -r %d -k %d -m %d',recon_cal_size,...
            espirit_kernel_size,espirit_num_maps),us_ksp_data_echo_n_slice); %bart('ecalib -S', us_img_slice_dim4);
        % size of input 1 = size of input 2
        % l1-wavelet, l2 regularization
        % r: regularization parameter (0.001 in Uecker 2014 paper:
        % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4142121/)
        % 'pics -l1 -r0.01
        espirit_coilimg = squeeze(bart(sprintf('pics -%s -r%d',reg_method,...
            reg_stepsize), us_ksp_data_echo_n_slice, espirit_maps));
        espirit_finalimg = bart('rss 4', espirit_coilimg); % 16 if 'squeeze' not used above
        
        if recon_along_x && plot_recon_figs
            espirit_map1_finalimg = abs(espirit_coilimg(:,:,1));
            % Transpose for images
            fs_finalimg = fs_finalimg';
            zf_finalimg = zf_finalimg';
            espirit_map1_finalimg = espirit_map1_finalimg';
            espirit_finalimg = espirit_finalimg';
            sense_finalimg = sense_finalimg';
        end
        
        %% Figures
        if ismember(slice,slices_to_compare)
            if show_recon_steps
                fig = figure;
                ax = gca;
                
                imshow(abs(fs_finalimg'),[])
                % For some reason, when you've used imshow and you don't put in the current
                % axis, the title appears on the previous axis
                title(ax,'(a) Original image','FontSize',16)
                set(fig,'Position',[50 600 300 250])
                %axis image ij
                % Plot the original k space, sampling mask, and undersampled k space for
                % previously sliced data
                % If data were generated multiple times, the k space/calib maps from the
                % last iteration will be plotted
                fig2 = figure;
                s1 = subplot(1,3,1);
                
                %imagesc(abs(squeeze(ksp_data_echo_n_z1(:,:,coil))))
                imagesc(abs(squeeze(ksp_data_echo_n_z1(:,:,coil)))')
                
                title(sprintf('(a) Original K-space (coil %d)',coil),'FontSize',16)
                set(s1,'Position',[0.03 ksp_bwh])
                
                s2 = subplot(1,3,2);
                imagesc(sampling_pattern');
                title('(b) Sampling pattern','FontSize',16)
                set(s2,'Position',[0.36 ksp_bwh])
                
                s3 = subplot(1,3,3);
                
                imagesc(squeeze(abs(us_ksp_echo_n_z1(:,:,1,coil)))')
                title(sprintf('(c) Undersampled K-space (coil %d)',coil),'FontSize',16)
                set(s3,'Position',[0.7 ksp_bwh])
                
                if recon_along_x
                    s = suptitle(sprintf('Slice %d of x-axis',slice));
                else
                    s = suptitle(sprintf('Slice %d of z-axis',slice));
                end
                
                set(s,'FontSize',16,'FontWeight','bold')
                set(fig2,'Position',[300 600 800 270])
                
                % View ESPIRiT maps
                fig3 = figure;
                ax = gca;
                
                imshow3(abs(squeeze(espirit_maps)),[],[espirit_num_maps,num_coils])
                title(ax,'ESPIRiT maps','FontSize',15)
                set(fig3,'Position',[50 100 400 130*espirit_num_maps])
                set(ax,'Position',[0.05 0.04 0.9 0.88])
            end
            
            if plot_recon_figs
                posn = zeros(4,4);
                % Comparisons of different reconstruction methods
                % Plot images from last iteration
                % Zero-filled reconstruction
                fig4 = figure;
                subplot(2,3,[1 4])
                fs_ax = gca;
                imshow(fs_finalimg,[]);colormap('gray')
                title('(a) Original image','FontSize',16);
                
                subplot(2,3,2)
                ax = gca;
                imshow(zf_finalimg,[])
                title(ax,'(b) Zero-filled recon','FontSize',16);
                posn(1,:) = get(ax,'Position');
                
                subplot(2,3,3)
                ax(2) = gca;
                imshow(espirit_map1_finalimg', [])
                title(ax(2),'(c) ESPIRiT recon (map 1)','FontSize',16)
                posn(2,:) = get(ax(2),'Position');
                
                % This image should be the most accurate/closest to image from fully
                % sampled k space
                subplot(2,3,5)
                ax(3) = gca;
                imshow(espirit_finalimg,[]);
                title(ax(3),'(d) ESPIRiT recon (all maps)','FontSize',16)
                posn(3,:) = get(ax(3),'Position');
                
                if gen_sense_recon
                    
                    subplot(2,3,6)
                    ax(4) = gca;
                    imshow(abs(sense_finalimg),[])
                    title(ax(4),'(e) SENSE recon','FontSize',16)
                    set(fig4,'Position',[400 100 1000 625])
                    posn(4,:) = get(ax(4),'Position');
                end
                
                % Reposition image axes
                % Shift to same column position
                posn([1,3],1) = 0.335;posn([2,4],1) = 0.67;
                % Shift row position
                posn([1,2],2) = posn([1,2],2)-0.04;posn([3,4],2) = posn([3,4],2)-0.04;
                posn(:,[3,4]) = repmat(recon_wh,4,1);
                fs_posn = get(fs_ax,'Position'); fs_posn=[0 0.31 recon_wh];
                set(fs_ax,'Position',fs_posn)
                for ii = 1:length(ax)
                    set(ax(ii),'Position',posn(ii,:))
                end
                
                if save_recon_plots
                    %set(gcf,'Color','none')
                    saveas(fig4,sprintf('%s/Slice %d.png',recon_img_dirpath,slice))
                end
            end
            
            if plot_poster_figs
                posn = zeros(3,4);
                % Comparisons of different reconstruction methods
                % Top: image from fully sampled k space
                % Middle: image from undersampled k space, zero-filled
                % reconstruction
                % Bottom: image from undersampled k space, ESPIRiT
                % reconstruction
                % Plot images from last iteration
                % Zero-filled reconstruction
                poster_fig = figure;
                subplot(3,1,1)
                ax = gca;
                imshow(fs_finalimg',[]);colormap('gray')
                %             ylabel('(a) Original image','FontSize',22,'FontWeight',...
                %                 'bold','Position',[-5 96 0])
                posn(1,:) = get(ax(1),'Position');
                
                subplot(3,1,2)
                ax(2) = gca;
                imshow(zf_finalimg',[])
                %             ylabel({'(b) Zero-filled','reconstruction'},'FontSize',...
                %                 22,'FontWeight','bold','Position',[-5 96 0])
                posn(2,:) = get(ax(2),'Position');
                
                % This image should be the most accurate/closest to image from fully
                % sampled k space
                subplot(3,1,3)
                ax(3) = gca;
                imshow(espirit_finalimg',[]);
                %             ylabel('(c) CS reconstruction','FontSize',22,...
                %                 'FontWeight','bold','Position',[-5 96 0])
                posn(3,:) = get(ax(3),'Position');
                
                % Shift to same column position
                posn(:,1) = 0.05;
                % Shift to same row position
                posn(1,2) = 0.67; posn(2,2) = 0.34; posn(3,2) = 0.01;
                posn(:,[3,4]) = repmat(poster_wh,3,1);
                
                if recon_along_x
                    set(poster_fig,'Position',[400 100 500 900])
                else
                    set(poster_fig,'Position',[400 100 400 900])
                end
                for ii = 1:size(posn,1)
                    set(ax(ii),'Position',posn(ii,:))
                end
            end
        end
        
        if recon_along_x
            if plot_recon_figs % Transpose back
                espirit_map1_finalimg = abs(espirit_coilimg(:,:,1));
                fs_finalimg = fs_finalimg';
                zf_finalimg = zf_finalimg';
                espirit_map1_finalimg = espirit_map1_finalimg';
                espirit_finalimg = espirit_finalimg';
                sense_finalimg = sense_finalimg';
            end
            % Add current x-axis slice to 3D arrays
            fs_3d_img(slice,:,:) = fs_finalimg;
            zf_3d_img(slice,:,:) = zf_finalimg;
            espirit_3d_img(slice,:,:) = espirit_finalimg;
            if gen_sense_recon
                sense_3d_img(slice,:,:) = sense_finalimg;
            end
        else
            % Add current z-axis slice to 3D arrays
            fs_3d_img(:,:,slice) = fs_finalimg;
            zf_3d_img(:,:,slice) = zf_finalimg;
            espirit_3d_img(:,:,slice) = espirit_finalimg;
            if gen_sense_recon
                sense_3d_img(:,:,slice) = sense_finalimg;
            end
        end
        
        
        if plot_diff_maps && ismember(slice,slices_to_compare)
            
            fs_rss = permute(fs_finalimg,[3 1 2]); % Add singleton dimension to match other images
            
            % Fully sampled image
            x = abs(squeeze(gs_normalize(fs_rss,255)));
            % Zero-filled reconstruction
            x_hat_zero = abs(gs_normalize(zf_finalimg,255));
            % ESPIRiT reconstruction (all combined maps)
            x_hat_espirit = abs(gs_normalize(espirit_finalimg,255));
            if gen_sense_recon
                % SENSE reconstruction
                x_hat_sense = abs(gs_normalize(sense_finalimg,255));
                all_images = [abs(x-x_hat_zero)',...
                abs(x-x_hat_espirit)',...
                abs(x-x_hat_sense)'];
                title_str = sprintf('(a) Zero-filled recon\t\t\t\t\t(b) ESPIRiT recon (%d maps)\t\t\t\t\t\t(c) SENSE recon',...
                    espirit_num_maps);
            else
                all_images = [abs(x-x_hat_zero)',...
                abs(x-x_hat_espirit)'];
                title_str = sprintf('(a) Zero-filled recon\t\t\t\t\t(b) ESPIRiT recon (%d maps)',...
                    espirit_num_maps);
            end
            % Also assess using difference maps
            fig5 = figure;ax = gca;imshow(all_images,[]);
            % For imshow, have to manually change current axes, otherwise the titles
            % plot on the prev axes, and suptitle doesn't accept axes as an input
            axes(ax);ax.FontSize = 13;
            title(title_str);
            t = suptitle('Difference maps');set(t,'FontSize',16,'FontWeight','bold')
            set(ax,'Position',[0.01 0.1 0.98 0.6])
            set(fig5,'Position',[150 400 500 150])
            
            if plot_poster_figs
                clear posn ax
                % 2-panel plot:
                % Top: image from undersampled, zero-filled reconstruction
                % Bottom: image from undersampled ESPIRiT reconstruction
                % Plot images from last iteration
                % Zero-filled reconstruction
                poster_fig2 = figure;
                subplot(2,1,1)
                ax = gca;
                if recon_along_x
                    imshow(abs(x-x_hat_zero)',[]);
                else
                    imshow(abs(x-x_hat_zero),[]);
                end
                colormap('gray')
                %             ylabel({'(a) Zero-filled','reconstruction'},'FontSize',22,...
                %                 'FontWeight','bold','Position',[-5 96 0])
                posn(1,:) = get(ax(1),'Position');
                
                subplot(2,1,2)
                ax(2) = gca;
                if recon_along_x
                    imshow(abs(x-x_hat_espirit)',[])
                else
                    imshow(abs(x-x_hat_espirit),[])
                end
                %             ylabel('(b) CS reconstruction','FontSize',22,...
                %                 'FontWeight','bold','Position',[-5 96 0])
                posn(2,:) = get(ax(2),'Position');
                
                % Shift to same column position
                posn(:,1) = 0.08;
                % Shift row position
                posn(2,2) = 0.03;
                if recon_along_x
                    set(poster_fig2,'Position',[100 100 500 500])
                    posn(:,[3,4]) = repmat([0.9,0.45],2,1);
                    posn(1,2) = 0.51;
                else
                    set(poster_fig2,'Position',[100 100 350 500])
                    posn(:,[3,4]) = repmat([0.9,0.39],2,1);
                    posn(1,2) = 0.61;
                end
                
                for ii = 1:size(posn,1)
                    set(ax(ii),'Position',posn(ii,:))
                end
            end
        end
        fprintf(sprintf('Reconstructed slice %d\n',slice))
    end
    time_stop(iter) = toc(time_start(iter));
    
    if save_nifti
        fs_nii = make_nii(fs_3d_img,[],ksp_dims(1:3),64);
        save_nii(fs_nii,sprintf('%s/fs_img.nii',main_dir));
        
        espirit_nii = make_nii(espirit_3d_img, [], ksp_dims(1:3), 64);
        save_nii(espirit_nii,sprintf('%s/espirit_recon.nii',main_dir));
    end
    
    if gen_quality_metrics
        for slice = slices_to_compare
            % Assess quality of reconstruction for a slice
            % https://onlinelibrary.wiley.com/doi/full/10.1002/ima.22260
            % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5240324/
            % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4917755/
            % Relative L2-norm error (RLNE), high frequency error norm (HFEN) (?),
            % peak signal-to-noise ratio (PSNR), structural similarity (SSIM)
            % For these metrics, normalize data to range from 0-255
            
            % fs_rss = permute(fs_img(slice,:,:),[3 1 2]); % Add singleton dimension to match other images
            
            if recon_along_x
                fs_rss = squeeze(fs_3d_img(slice,:,:)); % Add singleton dimension to match other images
                zf_img = squeeze(zf_3d_img(slice,:,:));
                espirit_img = squeeze(espirit_3d_img(slice,:,:));
            else
                fs_rss = squeeze(fs_3d_img(:,:,slice)); % Add singleton dimension to match other images
                zf_img = squeeze(zf_3d_img(:,:,slice));
                espirit_img = squeeze(espirit_3d_img(:,:,slice));
            end
            
            % !!! Update code to calculate these even if this statement isn't
            % entered, because RLNE/PSNR is calculated outside of this for loop
            % Fully sampled image
            x = abs(squeeze(gs_normalize(fs_rss,255)));
            % Zero-filled reconstruction
            x_hat_zero = abs(gs_normalize(zf_img,255));
            % ESPIRiT reconstruction (2 combined maps)
            x_hat_espirit = abs(gs_normalize(espirit_img,255));
            
            % RLNE: ||x_hat-x||/||x||.
            % Lower RLNE: reconstruction closer to fully sampled image
            rlne_espirit(slice,iter) = get_rlne(x,x_hat_espirit);
            
            % PSNR: 20*log10(255/sqrt(MSE)).
            % Higher PSNR: reconstructed pixel value more consistent to original image
            psnr_espirit(slice,iter) = get_psnr(x,x_hat_espirit,255);
            
            % SSIM from MATLAB function
            ssim_espirit(slice,iter) = ssim(x_hat_espirit,x);
            
            if iter == 1 || regenerate_mask
                rlne_zero(slice,iter) = get_rlne(x,x_hat_zero);
                psnr_zero(slice,iter) = get_psnr(x,x_hat_zero,255);
                ssim_zero(slice,iter) = ssim(x_hat_zero,x);
            end
            
            % SENSE recon optional
            if gen_sense_recon
                if recon_along_x
                    sense_img = squeeze(sense_3d_img(slice,:,:));
                else
                    sense_img = squeeze(sense_3d_img(:,:,slice));
                end
                % SENSE reconstruction
                x_hat_sense = abs(gs_normalize(sense_img,255));
                rlne_sense(slice,iter) = get_rlne(x,x_hat_sense);
                psnr_sense(slice,iter) = get_psnr(x,x_hat_sense,255);
                ssim_sense(slice,iter) = ssim(x_hat_sense,x);
            end
        end
    end
end

% Plot reconstruction quality performance
if plot_quality_metrics
    qm_fig = figure;
    
    % Save true acceleration rate to filename if not shown on plot
    if ~strcmp(iter_vars{iter_type},'accel')
        accel = actual_accel;
    else
        % Replace 'nominal' acceleration with 'true' acceleration for x
        % axis labels
        iter_values = actual_accel;
    end
    
    if ~regenerate_mask
        num_plots = 3;
        num_rows = 1;
        zf_cols = 'k';
        set(qm_fig,'Position',[50 300 950 300])
        legend_str = {'Zero-filled'};
        if strcmp(iter_vars{iter_type},'reg_stepsize')
            legend_str(end+1:end+length(iter_values)) = arrayfun(@(x) ...
                sprintf('ESPIRiT, %1.0e',x),iter_values,'UniformOutput',false);
        else
            legend_str(end+1:end+length(iter_values)) = arrayfun(@(x) ...
                sprintf('ESPIRiT, %d',x),iter_values,'UniformOutput',false);
        end
    else
        num_plots = 3*length(iter_values)*regenerate_mask;
        num_rows = ceil(num_plots/3);
        zf_cols = repmat('k',1,size(rlne_zero,2));
        cols = repmat('b',1,size(rlne_espirit,2));
        set(qm_fig,'Position',[50 300 950 min([280*num_rows,700])])
        if strcmp(iter_vars{iter_type},'accel')
            % True acceleration values are generally not integers, so
            % control output format
            legend_str = arrayfun(@(x) ...
                sprintf('ZF, %.2f',x),iter_values,'UniformOutput',false);
            legend_str(end+1:end+length(iter_values)) = arrayfun(@(x) ...
                sprintf('ESPIRiT, %.2f',x),iter_values,'UniformOutput',false);
        else
            legend_str = arrayfun(@(x) ...
                sprintf('ZF, %d',x),iter_values,'UniformOutput',false);
            legend_str(end+1:end+length(iter_values)) = arrayfun(@(x) ...
                sprintf('ESPIRiT, %d',x),iter_values,'UniformOutput',false);
        end
    end
    
    if gen_sense_recon
        legend_str2 = {'Zero-filled','ESPIRiT','SENSE'};
    else
        legend_str2 = {'Zero-filled','ESPIRiT'};
    end
    
    for ii = 1:num_plots
        s(ii) = subplot(num_rows,3,ii);hold(s(ii),'on')
    end
    
    iter_var_str = strrep(iter_vars{iter_type},'_',' ');
    % Plot individually so that line colors can be custom set
    % Plot each type of reconstruction together (for grouping in legend)
    for ii = 1:size(rlne_espirit,2)
        plot_ind = 1+3*regenerate_mask*(ii-1);
        col = cols(mod(ii-1,length(cols))+1);
        if floor((ii-1)/length(cols))>0
            line_style = ':';
        else
            line_style = '--';
        end
        if ii == 1 || regenerate_mask
            zf_col = zf_cols(mod(ii-1,length(zf_cols))+1);
            plot(s(plot_ind),rlne_zero(:,ii),zf_col,'LineWidth',1.5)
            plot(s(plot_ind+1),psnr_zero(:,ii),zf_col,'LineWidth',1.5)
            plot(s(plot_ind+2),ssim_zero(:,ii),zf_col,'LineWidth',1.5)
        end
        
        plot(s(plot_ind),rlne_espirit(:,ii),[line_style,col],'LineWidth',1.5)
        plot(s(plot_ind+1),psnr_espirit(:,ii),[line_style,col],'LineWidth',1.5)
        plot(s(plot_ind+2),ssim_espirit(:,ii),[line_style,col],'LineWidth',1.5)
        if strcmp(iter_vars{iter_type},'accel')
            ylabel(s(plot_ind),{sprintf('%s: %.2f',iter_var_str,iter_values(ii)),...
                '\leftarrow RLNE'})
        elseif strcmp(iter_vars{iter_type},'cal_reg')
            ylabel(s(plot_ind),{sprintf('%s: %d',iter_var_str,iter_values(ii)),...
                '\leftarrow RLNE'})
        else
            ylabel(s(plot_ind),'\leftarrow RLNE')
        end
        ylabel(s(plot_ind+1),'PSNR \rightarrow')
        ylabel(s(plot_ind+2),'SSIM \rightarrow')
        
%         ylabel(s(plot_ind),'\leftarrow Relative L2-Norm Error')
%         ylabel(s(plot_ind+1),'Peak Signal to Noise Ratio \rightarrow')
%         ylabel(s(plot_ind+2),'Structural Similarity Index \rightarrow')
    end
    
    for ii = 1:num_plots
        set(s(ii),'FontSize',13)
        if ii>num_plots-3
            xlabel(s(ii),xlbl)
        end
    end
    
    % Add values of all other variables to the filename (for notekeeping)
    non_iter = setdiff(1:length(iter_vars),iter_type);
    recon_str = '';
    
    for ii = non_iter
        eval(['curr_var = ',iter_vars{ii},';'])
        recon_ind = [1,(regexp(iter_vars{ii},'_')+1)];
        if strcmp(iter_vars{ii},'reg_stepsize')
            recon_str = [recon_str,'_',iter_vars{ii}(recon_ind),sprintf('%1.0g',curr_var)];
        elseif strcmp(iter_vars{ii},'accel')
            % If calibration region is changed, don't add acceleration to
            % filename (because that will also change)
            recon_str = [recon_str,'_',iter_vars{ii}(recon_ind),sprintf('%.2f',mean(curr_var))];
        else
            recon_str = [recon_str,'_',iter_vars{ii}(recon_ind),sprintf('%d',curr_var)];
        end
    end
    
    if gen_sense_recon
        if regenerate_mask % Line style for plot
            sense_ls = '--';
            sense_cols = repmat('r',1,size(rlne_espirit,2));
        else
            sense_ls = ':';
            sense_cols = cols;
        end
        for ii = 1:size(rlne_sense,2)
            % Plot in a new row if zero-filled data is plotted across
            % multiple rows
            plot_ind = 1+3*regenerate_mask*(ii-1);
            plot(s(plot_ind),rlne_sense(:,ii),[sense_ls,sense_cols(ii)],'LineWidth',1.5)
            plot(s(plot_ind+1),psnr_sense(:,ii),[sense_ls,sense_cols(ii)],'LineWidth',1.5)
            plot(s(plot_ind+2),ssim_sense(:,ii),[sense_ls,sense_cols(ii)],'LineWidth',1.5)
        end
        if strcmp(iter_vars{iter_type},'reg_stepsize')
            legend_str(end+1:end+length(iter_values)) = arrayfun(@(x) ...
                sprintf('SENSE, %1.0e',x),iter_values,'UniformOutput',false);
        elseif strcmp(iter_vars{iter_type},'accel')
            legend_str(end+1:end+length(iter_values)) = arrayfun(@(x) ...
                sprintf('SENSE, %.2f',x),iter_values,'UniformOutput',false);
        else
            legend_str(end+1:end+length(iter_values)) = arrayfun(@(x) ...
                sprintf('SENSE, %d',x),iter_values,'UniformOutput',false);
        end
        [~,min_rlne_ind] = min(mean([rlne_zero,rlne_espirit,rlne_sense]));
        title(s(1),['RLNE (best: ',legend_str{min_rlne_ind},')'])
        [~,max_psnr_ind] = max(mean([psnr_zero,psnr_espirit,psnr_sense]));
        title(s(2),['PSNR (best: ',legend_str{max_psnr_ind},')'])
        [~,max_ssim_ind] = max(mean([ssim_zero,ssim_espirit,ssim_sense]));
        title(s(3),['SSIM (best: ',legend_str{max_ssim_ind},')'])
        sense_flag = 'sense_';
    else
        [~,min_rlne_ind] = min(mean([rlne_zero,rlne_espirit]));
        title(s(1),['RLNE (best: ',legend_str{min_rlne_ind},')'])
        [~,max_psnr_ind] = max(mean([psnr_zero,psnr_espirit]));
        title(s(2),['PSNR (best: ',legend_str{max_psnr_ind},')'])
        [~,max_ssim_ind] = max(mean([ssim_zero,ssim_espirit]));
        title(s(3),['SSIM (best: ',legend_str{max_ssim_ind},')'])
        sense_flag = '';
    end
    
    if regenerate_mask
        l = legend(s(3),legend_str2);
        l.Location = 'northeast';
        l.Position(1) = 0.9;
    else
        l = legend(legend_str);
        if ~gen_sense_recon
            l.Position([1,2]) = [0.885,0.44];
        end
        for ii = 1:length(s)
            set(s(ii),'Position',[0.045+(ii-1)*0.295,0.13,0.24,0.7])
        end
        title(l,iter_var_str)
    end
    l.ItemTokenSize = [25,15];
    s_ttl = suptitle(['Reconstruction performance over ',lower(iter_var_titles{iter_type})]);
    set(s_ttl,'FontSize',19,'FontWeight','bold')%,'Position',[mean(s2_axl(1:2)),1.13*s2_axl(4)-0.13*s2_axl(3)])
    
    % Show effective acceleration when calibration region is changed in the
    % form of a caption in the first plot of each row
    % Set caption after suptitle; otherwise, they go off the edge of the
    % plot
    if strcmp(iter_vars{iter_type},'cal_reg')
        plot_ind = 1;
        for ii = 1:length(iter_values)
            % Add caption to each row
            axl = axis(s(plot_ind));
            text(s(plot_ind),0.05*axl(2),0.05*axl(3)+0.95*axl(4),...
                sprintf('accel: %.2f',actual_accel(ii)),'FontSize',12)
            plot_ind = plot_ind + 3;
        end
    end
    
    saveas(qm_fig,sprintf('%s/recon_qm_%s%s%s%s.fig',fig_dirpath,mask,sense_flag,xlbl(1),recon_str))
    saveas(qm_fig,sprintf('%s/recon_qm_%s%s%s%s.png',fig_dirpath,mask,sense_flag,xlbl(1),recon_str))
    
    % Generate plot of MEAN QM for each value over which you iterate
    if regenerate_mask
        num_repeats = 1;
    else
        num_repeats = length(iter_values);
    end
    
    qm_mean_fig = figure; set(qm_mean_fig,'Position',[50 300 700 300])
    s_ttl2 = suptitle(['Mean reconstruction performance over ',lower(iter_var_titles{iter_type})]);
    set(s_ttl2,'FontSize',19,'FontWeight','bold')
    s2(1) = subplot(1,3,1); plot(iter_values,repmat(mean(rlne_zero),1,...
        num_repeats),'k-x');hold on;
    plot(iter_values,mean(rlne_espirit),'-x');title('Relative L2-Norm Error')
    s2(2) = subplot(1,3,2); plot(iter_values,repmat(mean(psnr_zero),1,...
        num_repeats),'k-x');hold on; plot(iter_values,mean(psnr_espirit),'-x')
    title('Peak Signal to Noise Ratio')
    s2(3) = subplot(1,3,3); plot(iter_values,repmat(mean(ssim_zero),1,...
        num_repeats),'k-x');hold on; plot(iter_values,mean(ssim_espirit),'-x')
    title('Structural Similarity Index')
    if gen_sense_recon
        plot(s2(1),iter_values,mean(rlne_sense),'-x')
        plot(s2(2),iter_values,mean(psnr_sense),'-x')
        plot(s2(3),iter_values,mean(ssim_sense),'-x')
    end
    s2_axl = axis(s2(2));
    
    ylabel(s2(1),'\leftarrow RLNE')
    ylabel(s2(2),'PSNR \rightarrow')
    ylabel(s2(3),'SSIM \rightarrow')
    for ii = 1:3
        xlabel(s2(ii),iter_var_titles{iter_type});set(s2(ii),'FontSize',13)
        xlim(s2(ii),[min(iter_values) max(iter_values)])
        set(s2(ii),'Position',[0.07+(ii-1)*0.33,0.13,0.25,0.7])
    end
    legend(legend_str2)
    saveas(qm_mean_fig,sprintf('%s/mean_recon_qm_%s%s%s%s.fig',fig_dirpath,mask,sense_flag,xlbl(1),recon_str))
    saveas(qm_mean_fig,sprintf('%s/mean_recon_qm_%s%s%s%s.png',fig_dirpath,mask,sense_flag,xlbl(1),recon_str))
end

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

% fprintf(sprintf('\nRLNE:\n\tZero-filled recon: %.4f\n\tESPIRiT recon: %.4f\n\tSENSE recon: %.4f\n',...
%     rlne_zero,rlne_espirit,rlne_sense));
% fprintf(sprintf('\nSSIM:\n\tZero-filled recon: %.4f\n\tESPIRiT recon: %.4f\n\tSENSE recon: %.4f\n',...
%     ssim_zero,ssim_espirit,ssim_sense));
% fprintf(sprintf('\nPSNR:\n\tZero-filled recon: %.4f\n\tESPIRiT recon: %.4f\n\tSENSE recon: %.4f\n',...
%     psnr_zero,psnr_espirit,psnr_sense));

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