function img = RussReconADMask_ME(varargin)
%%function raw = Pseries(varargin)
% Before calling the function, put the pfiles (single- or multi-echo) in a
% separate folder. Call the function in the created folder. 
%
%   
% varargin is a list of options in 'single quotes' separated by commas
%   (no input)  - output is complex, raw, frequency space data
%   'raw'       - output is complex data in k-space
%   'echofill'  - replaces faulty data in last echo when phase1*phase2 is
%                 not divisible by baseline_spacing (usually 2048)
%   'shift'     - linearly corrects even echo shift to match odd echoes
%   [-1 -1 -1]  - flips readout, phase1, and phase2 array orientations
%                 (good for getting data into WHS if you scanned it weird)
%   'mask'      - creates a binary, skull-stripped mask
%   'BIAC#'     - saves data in format for susceptibility computation by
%                 the BIAC cluster
%                   
%
% Example: imagedata = Pseries4('img',[-1,1,-1]);
% imagedata will be an image space volume of complex data with the readout
% and phase2 directions reversed.
%
% img = Pseries4('img','shift','mask');
%
%
% Based on the function m-file, MapFreq.m, created by Chunlei Liu
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ORIG Input:
% function [Freq] = MapFreq()
% Input:     The function will read in all pfiles from the current
%            directory.
% Output:
% Freq  -    Frequency maps in Hz, stored in binary float format 
%            nrows x ncols x nslices with the file name of echo#.freq
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Chunlei Liu, Duke University, 06/2009
%
% slg - 100715 V2 removed Dim as param MapFreq, use
%       dimensions from pfile header (as stored by civm 3D cartesian scans)
%       Changed raw data reader to use header dimensions and bl skip.
%       V3 added fft2c and fftnc and inverses, provided by Chunlei.
%       Wrapped GE header reader to notice missing Pfile.
% slg - 100716 V4 Now able to discern raw format and read short and int/EDR 
%       Pfiles.
% rmd - 111025 V5 Now mimics the effect of radish reconstruction to pull
%       data from multiple, multi-echo pfiles and sort it into echoes
%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% User-defined parameters
SHORT_SIZE = 2; % in bytes
rawswitch = false;
cryomask = false;
echoflag = false;
shift9flag = false;
maskinflag = false;
fermi_flag = false;
% if matlabpool('size') > 0
%     pswitch = true;
% end
for narg = 1:nargin
    if ischar(varargin{narg})
        if strcmp(varargin{narg},'theworks')
            imgswitch = true;
            shiftswitch = true;
            phaseswitch = true;
            maskswitch = true;
        elseif strcmp(varargin{narg},'bruker')
            bruker_flag = true;
        elseif strcmp(varargin{narg},'cryomask')
            cryomask_flag = true;
        elseif strcmp(varargin{narg},'multimask')
            multimask_flag = true;
            maskswitch = true;
        elseif strcmp(varargin{narg},'Ashift9')
            shift9flag = true;
        elseif strcmp(varargin{narg},'agilent')
            agilent_postrecon = true;
        elseif strcmp(varargin{narg},'nomask')
            imgswitch = true;
            shiftswitch = true;
            phaseswitch = true;
        elseif strcmp(varargin{narg},'raw')
            rawswitch = true;
        elseif strcmp(varargin{narg},'X')
            Xswitch = true;
        elseif strcmp(varargin{narg},'post')
            postprocessingonly = true;
        elseif strcmp(varargin{narg},'zcorrect')
            zcorrect_flag = true;
        elseif strcmp(varargin{narg},'fermi')
            fermi_flag = true;
        elseif strcmp(varargin{narg},'parallel')
            pswitch = true;
            if matlabpool('size') == 0
                matlabpool open
            end
        elseif strcmp(varargin{narg},'rp')
            rp_flag = true;
            rp_file = dir('*.rp');
        elseif strcmp(varargin{narg},'T2')
            T2switch = true;
        elseif strcmp(varargin{narg},'shift')
            shiftswitch = true;
        elseif strcmp(varargin{narg},'phaseshift')
            phaseswitch = true;
        elseif strcmp(varargin{narg},'mask')
            maskswitch = true;
        elseif strcmp(varargin{narg},'center')
            centerswitch = true;
            imgswitch = true;
        elseif strcmp(varargin{narg},'int')
            pfile_type = 'int';
        elseif strcmp(varargin{narg},'save')
            matsave = true;
        elseif strcmp(varargin{narg},'short')
            pfile_type = 'short';    
        elseif strcmp(varargin{narg},'convert')
            SHORT_SIZE = 2; % in bytes
            convertflag = true;
        elseif strcmp(varargin{narg}(1:4),'mask')
            maskinflag = true;
            mask = open_nii(varargin{narg});
        elseif length(varargin{narg}) < 5
            scan_name = varargin{narg};
        elseif strcmp(varargin{narg}(1:4),'echo')
            echoflag = true;
            iecho = str2num(varargin{narg}(5:end));
        elseif strcmp(varargin{narg}(1:4),'BIAC')
            BIACswitch = true;
            maskswitch = true;
%             iBIAC = str2num(varargin{narg}(5:end));
            runno = str2num(varargin{narg}(5:end));
        elseif ~strcmp(varargin{narg}(1:4),'BIAC')
            scan_name = varargin{narg};
        end
    else
        if length(varargin{narg}) == 3    
            if length(varargin{narg}(abs(varargin{narg}) > 1)) > 1
                permuteswitch = true;
                permute_vector = varargin{narg};
                save permute_vector permute_vector
                permute_vector = [permute_vector 4];
            else
                orientswitch = true;
                orientation = varargin{narg};
            end
        end
    end   
end

%% Recon summary



%% Convert GE p-file series by removing baselines and creating echo files

if ~exist('scan_name','var')
    scan_name = '';
end

if exist('convertflag','var')
    Pconvert(scan_name);
    scan_name = [scan_name 'k'];
%     pfile_type = 'short';
elseif exist('rp_flag','var')
    [hdr,header_bytes] = WrapGEheader(rp_file(1).name);
    bp_half_element = hdr.rdb.point_size;
%     baseline_spacing = hdr.rdb.nframes;
    SHORT_SIZE = 2;  
%     flag='l';  % default byte swap
    if bp_half_element == SHORT_SIZE
        pfile_type = 'short';
    else
        if bp_half_element ~= 4
            error('\nPfile rdb.point_size unknown.\n');
        end
        pfile_type = 'int';
    end
    % Save scan paramters
    dims = [hdr.rdb.da_xres, hdr.rdb.user7, ...
        hdr.rdb.user8, hdr.rdb.nechoes]; save dims dims
    fov = [hdr.rdb.user16, hdr.rdb.user17, hdr.rdb.user18]; save fov fov
    vox = fov./dims(1:3); save vox vox   % voxel size, mm
    if  ((hdr.rdb.user33 == 0) && dims(4) > 1)
        hdr.rdb.user33 = waitinput('What is the echo spacing (in ms)? ',10);
        if isnan(hdr.rdb.user33)
            hdr.rdb.user33 = 2.896;
        end
    end
    echo_spacing = hdr.rdb.user33/1e3; % Only valid in GE multi-gradient echo
    TE = (0:dims(4)-1)*echo_spacing+hdr.rdb.te/1e6; save TE TE
    save hdr hdr header_bytes

    raw=zeros(dims(1),dims(2),dims(3));
    back = 0; time = 0;
    fid = fopen(rp_file(1).name,'r','l');
    [~]=fread(fid,[1 header_bytes/SHORT_SIZE],'short');  % skip header
    for z=1:dims(3)
        [back,time] = progress(z,dims(3),'Reading slice',back,time);
            for y=1:dims(2)    
                dump=fread(fid,[2 dims(1)], pfile_type);
                raw(:,y,z)=squeeze(dump(2,1:end)+1i*dump(1,:));
            end
    end
    fclose(fid);   
    
end


%% Read in raw data, (Agilent, Bruker or GE)

if exist('bruker_flag','var')   %%% Bruker Style
    
    B0 = 7;
    
    if ~exist('postprocessingonly','var')
        % [raw,hdr] = BrukerRecon(pwd);
       [raw,hdr] = BrukerReconRawMGEMS;
    else
        hdr.acqp = readBrukerHeader('acqp');
        hdr.method = readBrukerHeader('method');
%         raw = readBrukerFID('',hdr.method);
        hdr.subject = readBrukerHeader('../subject');
        load img2channels
        raw = permute(img,[1 2 3 5 4]); clear img;
     
    end
    
    if strcmp(hdr.method.PVM_SPackArrSliceOrient,'coronal')
        B0dir = [1 0 0]; 
        BIACpermute = [2 3 1];
    elseif strcmp(hdr.method.PVM_SPackArrSliceOrient,'axial')
        B0dir = [0 0 1];
        BIACpermute = [1 2 3];
    elseif strcmp(hdr.method.PVM_SPackArrSliceOrient,'sagittal')
        B0dir = [1 0 0];
        BIACpermute = [3 2 1];
    end
    scan_orientation = hdr.method.PVM_SPackArrSliceOrient;
    save B0dir B0dir BIACpermute scan_orientation
    
    
     
    FOV = hdr.method.PVM_Fov;
    raw = single(raw);

%     vox = hdr.method.PVM_SpatResol; 
     
    %%% Pad the array here if desired, e.g.: raw = padarray(raw,[0 0 88 0],'both');
%     dims = [hdr.method.PVM_Matrix hdr.method.PVM_NEchoImages];
    dims = size(raw); 
    dims = dims(1:4);
    vox = FOV./dims(1:3);
    if length(vox) < 3
        if strcmp(hdr.method.Method,'MDEFT');
            FOV(3) = hdr.method.PVM_SliceThick*hdr.method.PVM_SPackArrNSlices;
            TE = hdr.method.PVM_EchoTime1/1000; 
        else
            FOV(3) = hdr.method.PVM_SliceThick;
            TE = hdr.method.EffectiveTE/1000;
        end
        vox(3) = hdr.method.PVM_SliceThick;
    else
        TE = hdr.method.EffectiveTE/1000;
    end
    save vox vox
    % Corrects Bruker's shift in phase encode 2 due to bad pulses.
    if dims(4) > 1 
        if exist('zcorrect_flag','var')
            raw = shiftMEbruker(raw); 
        end
    end  
    
    % first echo, first coil
    % !!! set breakpoint here
    img = raw;
end
