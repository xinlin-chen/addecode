mGRE_refmask = '/Users/omega/APOE_Study/Raw_Data_Folders/B05216.work/B05216_mask.nii.gz' ;
mGRE_outmask = '/Users/omega/APOE_Study/Raw_Data_Folders/B05221.work/B05221_mask_test.nii.gz' ;
mGRE_out = '/Users/omega/APOE_Study/Raw_Data_Folders/B05221.work/B05221_FE.nii.gz'; 
mGRE_ref = '/Users/omega/APOE_Study/Raw_Data_Folders/B05216.work/B05216_FE.nii.gz';
outputTransformPrefix= '/Users/omega/APOE_Study/Raw_Data_Folders/B05221.work/B05221transf';
myAffine='/Users/omega/APOE_Study/Raw_Data_Folders/B05221.work/B05221transf0GenericAffine.mat';
GREmasked='/Users/omega/APOE_Study/Raw_Data_Folders/B05221.work/B05221GREmasked.nii.gz';
outputm2f='/Users/omega/APOE_Study/Raw_Data_Folders/B05221.work/B05221GREmask2NewMGRE.nii.gz';
%%
ref_nii = load_nii('Users/omega/APOE_Study/Raw_Data_Folders/B05216.work/B05216mGRE.nii.gz');
myimg2=abs(ref_nii.img(:,:,:,2));
myGRE = make_nii(myimg2,[0.1,0.1,0.1],[0, 0 ,0],16);
save_nii(myGRE, 'B05216_FE.nii.gz');
%%
cmd1=['/Applications/ants_201410/antsRegistration -d 3  -o [ ' outputTransformPrefix ',' outputm2f ']' ' -m MI[' mGRE_out ',' mGRE_ref ',1,32,Regular,1] -t Affine[0.1] -c [ 100x100x20,1e-6,20] ' ...
    ' -s 4x4x4vox -f 8x4x1 -r ['  mGRE_out ',' mGRE_ref ',1 ] -u -v'  ];
cmd2=['/Applications/ants_201410/' 'antsApplyTransforms -d 3 -e 0 -i ' mGRE_refmask  ' -r ' mGRE_out  ' -o ' mGRE_outmask ' -n NearestNeighbor -t ' myAffine ' --float'];
cmd3=['/Applications/ants_201410/ImageMath 3 ' GREmasked ' m ' mGRE_outmask ' ' mGRE_out];

[res1,msg1]=system(cmd1,'-echo');
[res2,msg2]=system(cmd2,'-echo');
[res3,msg3]=system(cmd3,'-echo');
%%
output_dir= '/Users/omega/APOE_Study/Raw_Data_Folders/B05221.work/';

in_mask=mGRE_outmask; 
input_image = GREmasked;
runnos={'B05221'};
runno=char(runnos{1});
radius=num2str(2);

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
