function [res]=civm_to_nii(img_dir, file_name, file_ext, out_nii, dimx, dimy, dimz, format, vsizex, vsizey,vsizez, flipx, flipz, startzcrop, endzcrop, roll_string)
% function [res]=civm_to_nii(img_dir, file_name, file_ext, out_nii,
% dimx, dimy, dimz, format, vsizex, vsizey,vsizez, 
% flipx, flipz, startzcrop, endzcrop, roll_string)
% All arguments are required, however when used to reprocess nifti
% files, 
% the following fields will be auto-populated from the input nifti.
%   format
%   vsizex
%   vsizey
%   vsizez
%   dimx
%   dimy
%   dimz
%example call
% [res]=civm_to_nii('/Volumes/atlas2/15.abb.07/B04114/', 'B04114', 'raw', 'tmp.nii', 200,180,100, 512, .1, .1,.1,0,0,0);
%[res]=civm_to_nii('/Volumes/atlas2/15.abb.07/B04114/', 'B04114', 'raw', '/Users/omega/Natalie/recons/B04114.nii', 200,180,100, 512, .1, .1,.1,0,0,0);

addpath('/Users/omega/Documents/MATLAB/NIfTI_20140122');
if nargin<11
    error('not enought arguments, need at least 11 arguments');
end
if ~exist('roll_string','var')
  roll_string='0:0';
end
fprintf('using niitools at %s\n',which('save_nii'));
%niitools='/Volumes/pipe_home/matlab_library/NIFTI_20090909';
%if exist(niitools,'dir')
%    addpath('/Volumes/pipe_home/matlab_library/NIFTI_20090909')
%end

tic
% sample call for full res set
% civm_to_nii_feb ( '/cretespace/N33054Labels-inputs/N33055', 'N33055t9imx.0', 'raw', '/cretespace/N33054Labels-work/T2star-N33055.nii', 1024, 512, 512, 4, 1, 1 )
% sample call for 43 res dayta set
% civm_to_nii_may ( '/cretespace/N33820Labels-inputs/N33820', 'N33820t9imx.0', 'raw', '/cretespace/N33820Labels-work/T2star-N33820_flipz.nii', 512, 256, 256, 4, 2, 0,1 )


% from civmtonii_atlas2
%{
091118 slg fix
LOG: matlab mfile contains: civmtonii_atlas2 ( '/cretespace/N32817Labels-inputs/N32817', 'N32817t9imx.0', 'raw', '/cretespace/N32817Labels-work/T1.nii', 1024, 512, 512, 4, 1 )
LOG: Executing nifti conversion: /usr/bin/matlab < /cretespace/N32817Labels-work/T1civmtonii_atlas2 > /tmp/matlab_pipe_stuff
??? Error using ==> reshape
To RESHAPE the number of elements must not change.

Error in ==> civmtonii_atlas2 at 66
    vol(:, :, i)=reshape(im, dimx, dimy);
%}

%function to save civm image files  to nii
%sample call:
%[res]=civm_to_nii_feb('/Volumes/atlas1/07.soderling.01/N31300',
%'N31300t9imx.0', 'raw', 'tmp.nii', 512, 256, 256, 512, 2, 1);

%options for format parameter
%		2 - uint8,  4 - int16,  8 - int32,  16 - float32,
%		32 - complex64,  64 - float64,  128 - RGB24,
%		256 - int8,  512 - uint16,  768 - uint32, 
%		1792 - complex128
 %change voxel spacing for different resolution images...
 %voxel size for T2W is half the voxel size in T1, T2s therefore vsize should be set to 2

 
 %example call
 %[res]=civm_to_nii_feb('07.soderling.01', 'N31300', 'raw', 'tmp.nii', 512,256, 256, 512, 2, 1);
% [res]=civm_to_nii('15.abb.07', 'B04114', 'raw', 'tmp.nii', 200,180,100, 512, .1, .1,.1);

%% load
if strcmp(file_ext,'rawl')
% CIVM data that is little endian is notated .rawl (i.e. our ct images) vs .raw
% this script assumes all other input formats are big endian
        error('nifti creator only handles ieee-be data right now, not .rawl images'); 
elseif strcmp(file_ext,'raw')
    %disp(format)
    switch format
        case 2
            datatype='uint8';
        case 4
            datatype='int16';
        case 8
            datatype='int32';
        case 16
            datatype='single';
        case 256
            datatype='int8';
        case 512
            datatype='uint16';
        case 768
            datatype='uint32';
        otherwise
            res=-1;
            error('nifti: you have to specify data types as: 2 - uint8,  4 - int16,  8 - int32,  16 - float32');
    end
    
    % --- read raw files and make nifti volume
    
    vol=zeros(dimx, dimy, dimz);
    
    npix = dimx * dimy;
    if exist('read3d','file')==2
        pf=strsplit(file_name,'.');
        pf=pf{1};
        vol=read3d_debugging(dimx,dimy,dimz,[img_dir '/' pf ],datatype,'b');
        %bt7imx
    else
        for i=1:dimz
            switch 1
                case {i <= 9}
                    fname=strcat(img_dir, '/',  file_name, 'bt7imx.', '000', num2str(i), '.', file_ext);
                case {i > 9 & i <= 99 }
                    fname=strcat(img_dir, '/',  file_name,'bt7imx.', '00' , num2str(i), '.', file_ext);
                case {i > 99 & i <= 999}
                    fname=strcat(img_dir, '/',  file_name,'bt7imx.', '0'  , num2str(i), '.', file_ext);
                otherwise
                    'you are out of my file numbering range dear Alice'
            end
            %disp(fname);
            fid = fopen(fname, 'r', 'ieee-be');
            
            im=fread(fid, npix, datatype);
            [x,y] = size(im);
            if (x ~= npix)
                dbug (1,1, 'slice %d, got %d pix, need %d pix', i, x, npix);
                'WARNING input image slice read in nifti converter did not give enough data!'
                %make_matlab_error_file(code, 'error in to civm_to_nii: %s read %d pix, need %d pix', fname,  x, npix);
                return;
            end
            
            im=reshape(im, dimx, dimy);
            vol(:,:,i)=im;
            fclose(fid);
        end
    end
elseif strcmp(file_ext,'nii')
   %%% must load niis
%   error('NOt loading niis yet');
   fname=strcat(img_dir, '/',  file_name, '.', file_ext);
   nii=load_untouch_nii(fname);
   pixsize=nii.hdr.dime.pixdim(2:4);
   format=nii.hdr.dime.datatype;
   vsizex=pixsize(1);
   vsizey=pixsize(2);
   vsizez=pixsize(3);
   vol=nii.img;
   dimx=size(vol,1);
   dimy=size(vol,2);
   dimz=size(vol,3);
   % because start and end crop can be used to reverse the z direction this cant be a simple check.
   if  max([startzcrop endzcrop]) - min([startzcrop endzcrop])==0 
     % was a check for 1, but that was wrong, perhaps there was some reason for this.
        clear starzcrop endzcrop;
   end
else
    error('Unrecognized file_ext %s',file_ext);
end
%% roll
rolls=strsplit(roll_string,':');
vol=circshift(vol,[str2double(rolls{1}) str2double(rolls{2})]);

%% flip
if flipx==1
    disp('rotating to "flip" x');
    vol=imrotate(vol,180);
else
end
if flipz==1
    disp('rotating to "flip" z');
    %vol=flipdim(vol,1);
    vol=flipdim(vol,3);
    vol=flipdim(vol,2);
    
end
%% crop
if exist('startzcrop','var') && exist('endzcrop','var')
    if and((startzcrop>endzcrop) , startzcrop ~= dimz)
        inc=-1;
        %ind=[startzcrop:sz,startzcrop-1:inc:endzcrop,0:endzcrop-1]
        ind=[startzcrop:dimz,1:endzcrop]
        vol=vol(:,:,ind);
    elseif and((startzcrop>endzcrop) , startzcrop == dimz)
        inc=-1;
        vol=vol(:,:,startzcrop:inc:endzcrop);
    else
        inc=1;
        vol=vol(:,:,startzcrop:inc:endzcrop);
    end
    %vol=vol(:,:,startzcrop:inc:endzcrop);
    
end
%% save
% it appears  nii format indicates data endianness by the endianness of the nifti header (?!),
% so no header field to set to specify endianness
if format==0 
  error('out data type set incorrect');
end

nii = make_nii(vol, [vsizex vsizey vsizez], [round(dimx/2) round(dimy/2) round(dimz/2)], format);

% nii = make_nii(vol, [vsize vsize vsize], [0 0 0], format);
try
    save_nii(nii, out_nii);
    disp('saved nii');
catch me
    %# report the problematic image, and the reason for failure
    error('Nii creation failure % : %s\n',out_nii,me.message)
end
res=1;
fprintf('Total civm_to_nii time:%f\n', toc);    
