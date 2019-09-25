%Making Nifti's for RARE Recons
addpath('/Users/omega/APOE_Study/')
postfix = {'bt7'};
runnums = zeros(1,44);
runnums(1,1) = 5265;
for i = 2:44
    runnums(1,i) = runnums(1,i-1) + 5;
end
    
for i = 1:length(runnums)
    file_name = ['B0' num2str(runnums(i))];
    file_name_images = [file_name 'images'];
    img_dir = [ '/Users/omega/APOE_Study/' file_name '/' file_name_images ];
    savednifti=['/Users/omega/APOE_Study/RARE_APOE/' file_name '.nii'];
   [res]=civm_to_nii(img_dir, file_name, 'raw' , savednifti, 200,180,100, 512, .1, .1, .1, 0, 0, 0);
end
%% SHIFTING OF RARE
addpath('/Users/omega/Documents/MATLAB/NIfTI_20140122') 
rare =load_nii('/Users/omega/APOE_Study/RARE_APOE/B05465.nii');
rareimg = rare.img;
shift_rareimg = circshift(rareimg,80,3);
new_rareimg = make_nii(shift_rareimg);
save_nii(new_rareimg,'/Users/omega/APOE_Study/RARE_APOE/B05465_shift.nii');

%% Skull Stripping
nums = 5480;
img_dir = [ '/Users/omega/APOE_Study/RARE_APOE/'  ['B0' num2str(nums)] '.nii' ];
skullsave = ['/Users/omega/APOE_Study/skull_stripped_RARE/'  ['B0' num2str(nums)] 'strip_mask.nii'];
strip_mask_alex_mm(img_dir, 1, 9500, skullsave, 2 , 5,0 );