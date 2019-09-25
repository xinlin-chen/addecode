antspath='/Applications/ants_201410/';
RARE = '/Users/omega/APOE_Study/RARE_APOE/B05220.nii';
GRE = '/Users/omega/APOE_Study/Raw_Data_Folders/B05221.work/B05221_FE.nii.gz';
mask = '/Users/omega/APOE_Study/RARE_Masks/B05220_mask.nii';
mymask='/Users/omega/APOE_Study/RARE_Masks/B05220_mask2.nii';
outputTransformPrefix= '/Users/omega/APOE_Study/RARE_Masks/B05221transf'
myAffine='/Users/omega/APOE_Study/RARE_Masks/B05221transf0GenericAffine.mat'
outputWarpedImage='/Users/omega/APOE_Study/RARE_Masks/B05220Warped2RARE.nii.gz'
outputInverseWarpedImage='/Users/omega/APOE_Study/RARE_Masks/B05221bInvWarped.nii.gz'

outputm2f='/Users/omega/APOE_Study/RARE_Masks/B05221RARE2GRE.nii.gz'

GREmasked='/Users/omega/APOE_Study/RARE_Masks/B05221GREmasked.nii.gz'
RAREmasked='/Users/omega/APOE_Study/RARE_Masks/B05220RAREmasked.nii.gz'
GREmask='/Users/omega/APOE_Study/RARE_Masks/B05221GREmask.nii.gz';

mymask='/Users/omega/APOE_Study/RARE_Masks/B05220_mask2.nii';
RAREmask=mymask;

output_dir= '/Users/omega/APOE_Study/RARE_Masks/Corrected/';

input_image=RAREmasked; 
runnos={'B05220';'B05221'};
radius=1;

for i=2:2
    runno=char(runnos{i})
    if i==1 
        input_image=RAREmasked;
        in_mask=mymask; %
    elseif i==2
        input_image=GREmasked
        in_mask=GREmask; % 
    end
% mask_dir=???;
% runno=???;
% input_image =???;



dilated_mask=[output_dir runno '_dilated_mask_r' radius '.nii.gz'];
bfc_image = [ output_dir runno '_bfc.nii.gz'];
rescaled_bfc_image = [ output_dir runno '_bfc_rs.nii.gz'];
bfc_field = [ output_dir runno '_bfc_field.nii.gz'];
bfc_masked_image = [ output_dir runno '_BFmasked.nii.gz'];

convergence = '[2000x2000,1e-6]';
spline = '[2x2x2,3]';
shrink = '4';

%antscmd_p5 = [antspath 'CopyImageHeaderInformation ' mgfile ' ' in_mask ' ' in_mask ' 1 1 1'];
antscmd=[antspath 'ImageMath 3 ' dilated_mask  ' MD ' in_mask ' ' radius];
[res,msg]=system(antscmd,'-echo')
antscmd2=[antspath 'N4BiasFieldCorrection 3' ' -v 1 -i ' input_image  ' -x ' dilated_mask  ' -s ' shrink ' -c ' convergence  ' -b ' spline ' -o [' bfc_image ','  bfc_field ']'];
[res,msg]=system(antscmd2,'-echo')
rescale_cmd = [antspath 'ImageMath 3 ' rescaled_bfc_image ' RescaleImage ' bfc_image ' 0 65536'];
[res,msg]=system(rescale_cmd,'-echo')
convert_cmd = [antspath 'ConvertImagePixelType ' rescaled_bfc_image ' ' rescaled_bfc_image ' 3'];
[res,msg]=system(convert_cmd,'-echo');
cmd= [antspath 'ImageMath 3 ' bfc_masked_image  ' m ' in_mask ' ' rescaled_bfc_image]; % Run this once you have a satisfactory mask
[res,msg]=system(cmd,'-echo');
end


%second step registration for better masking

outputTransformPrefix2= '/Users/omega/APOE_Study/RARE_Masks/Corrected/B05221transf2'
myAffine2='/Users/omega/APOE_Study/RARE_Masks/Corrected/B05221transf20GenericAffine.mat'
outputWarpedImage2='/Users/omega/APOE_Study/RARE_Masks/Corrected/B05220Warped2RARE2.nii.gz'
outputInverseWarpedImage2='/Users/omega/APOE_Study/RARE_Masks/Corrected/B05221bInvWarped2.nii.gz'
GREmasked=[output_dir char(runnos{2}) '_BFmasked.nii.gz'] 
RAREmasked=[output_dir char(runnos{1}) '_BFmasked.nii.gz'] 


GREmask2=[output_dir char(runnos{2}) '_dilated_mask_r.nii.gz']
GREmasked2=[output_dir char(runnos{2}) '_BFmasked2_undil.nii.gz'] ; 

RAREmasked2=[output_dir char(runnos{1}) '_BFmasked2.nii.gz'] ; 
outputm2f2=[output_dir char(runnos{1}) '_2GRE2.nii.gz'] ; 

cmd12=['/Applications/ants_201410/antsRegistration -d 3  -o [ ' outputTransformPrefix2 ',' outputm2f2 ']' ' -m MI[' GREmasked ',' RAREmasked ',1,32,Regular,1] -t Affine[0.1] -c [ 100x100x20,1e-6,20] ' ...
    ' -s 4x4x4vox -f 8x4x1 -r ['  GREmasked ',' RAREmasked ',1 ] -u -v'  ];

cmd22=['/Applications/ants_201410/' 'antsApplyTransforms -d 3 -e 0 -i ' mymask  ' -r ' GREmasked  ' -o ' GREmask2 ' -n NearestNeighbor -t ' myAffine2 ' --float'];
%cmd2b=['/Applications/ants_201410/' 'antsApplyTransforms -d 3 -e 0 -i ' RAREmasked  ' -r ' GREmasked  ' -o ' outputm2f2 ' -n BSpline -t ' myAffine2 ' --float'];


cmd32=['/Applications/ants_201410/ImageMath 3 ' GREmasked2 ' m ' GREmask ' ' GRE ]

[res,msg]=system(cmd12,'-echo')
[res,msg]=system(cmd22,'-echo')
[res,msg]=system(cmd32,'-echo')



%% Use system(cmd), etc  to actually run all the above commands.