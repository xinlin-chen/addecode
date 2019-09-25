function im=read3d(dim1,dim2,dim3,filename,datatype,endian)
% Syntax:  im = read3d(dim1,dim2,dim3, filename,'datatype','endian');
% This function reads from disk the specified 2/3D magnitude image,
%
% === Inputs ===
% dim1
% dim2
% dim3
% filename: a the full filepath or a fullpath with a prefix
% datatype: number precision, 'float', 'uint16','real*4','ubit16', 'ubit8', etc.
% see doc fread for full list of precisions allowed.
% endian = l|ieee-le|b|ieee-be
%
% === Outputs ===
% im: a 3-dimensional matrix of size dim1 x dim2 x dim3 in the
% workspace ready to be used
%
% === Examples ===
% single slice data
% read3d(880,440,440,'/S50914usimx','uint16','b');
% volume data
% read3d(880,440,440,'/large3dvolume.raw','uint16','b');

if(nargin < 6) % assume big endian(mac default) if not specified
    endian = 'b';
    postfix = 'raw';
end

if (strcmp(endian,'l')||strcmp(endian,'ieee-le')) 
    postfix='rawl';
elseif (strcmp(endian,'b')||strcmp(endian,'ieee-be'))
    postfix='raw';
end
if strcmp(datatype,'float')||strcmp(datatype,'single') %||strcmp(datatype,'double')
    postfix='fp32';
    datatype='single';
end
supported_types={'single' 'int16' 'uint16'};

if ~regexp(datatype,strjoin(supported_types,'|'))
    warning('you chose an unexpected/unsupported type will try to just continue')
end

% interpret datatype
% not sure why ubit was used instead of uint, 16 and 8 are standard uint sizes.
% even more confusing, why take the first char of data type and then
% recreate data type from that?, lets just assume that bogus and comment
% these lines
% if (datatype(1)=='f') datatype='float'; end
% if (datatype(1)=='1') datatype='ubit16'; end
% if (datatype(1)=='8') datatype='ubit8'; end


%checkfor bad filename
fd=fopen(filename,'r',endian);% is this volume at a time
civmslicepattern='%s.%04i.%s';
ijslicepattern='%s%04i.%s';
%

% test for start at 0
% test for civmpattern
% test for ijpattern
% test for start at 1
% test for civmpattern
% test for ijpattern

slices=false;
volume=false;
attemptedfiles={};
if (fd<0)
    z=0;% slices starting at 0
    slicename=sprintf(civmslicepattern,filename,z,postfix)
    if(exist(slicename,'file') && ~slices)
        slicestart=0;
        slices=true;
        volume=false;
        pattern=civmslicepattern;
    end
    attemptedfiles{end+1}=slicename;
    slicename=sprintf(ijslicepattern,filename,z,postfix)
    if( exist(slicename,'file') && ~slices)
        slicestart=0;
        slices=true;
        volume=false;
        pattern=ijslicepattern;
    end
    attemptedfiles{end+1}=slicename;
    z=1;% slices starting at 1
    slicename=sprintf(civmslicepattern,filename,z,postfix)
    if(exist(slicename,'file') && ~slices)
        slicestart=1
        slices=true;
        volume=false;
        pattern=civmslicepattern;
    end
    attemptedfiles{end+1}=slicename;
    slicename=sprintf(ijslicepattern,filename,z,postfix);
    if(exist(slicename,'file') && ~slices)
        slicestart=1;
        slices=true;
        volume=false;
        pattern=ijslicepattern;
    end
    attemptedfiles{end+1}=slicename;
    if(~slices)
        string=sprintf('ERROR: cant find file %s or slice(s)\n',filename);
        for i=1:length(attemptedfiles)
            string=[string sprintf('\t%s\n',attemptedfiles{i}) ];
        end
        error('%s',string);
    end
else
    slices=false;
    volume=true;
    %fclose(fd);
end


im=zeros(dim1,dim2,dim3,datatype); % not sure classing this as uint16 will ultimately be useful given plan to write/re-write array
if(volume)
    %fd=fopen(filename,'r',endian);
    [im, countread]=fread(fd,dim1*dim2*dim3,datatype,0,endian);
    if(countread~=dim1*dim2*dim3)
        error('Dimension mismatch, make sure you are supplying proper dimensions to load.');
    end
    
    im=reshape(im, [dim1 dim2 dim3]); %since we did zeroes do we even
    %need reashape?
    fclose(fd);
end
if(slices) % this looks like a job for parfor! when tested in windows matlabpool failed.
%     if(matlabpool('size')==0)
%         matlabpool('open');
%     end
    for z=1:dim3 % why isnt this a straight fread, reshape operation? was
        % this code adapted from slice reading code.... Maybe we
        % should add that option.
        % filename would have to be a pattern read here...
        if(slicestart==0)
            slicename=sprintf(pattern,filename,z-1,postfix);
        else
            slicename=sprintf(pattern,filename,z,postfix);
        end
        %slicename;
        fd=fopen(slicename,'r',endian);
        im(:,:,z)=fread(fd,[dim1 dim2],datatype,0,endian);
        fclose(fd);
    end
%     matlabpool('close');
end
%im=im0;