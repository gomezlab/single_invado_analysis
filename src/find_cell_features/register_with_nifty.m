function register_with_nifty(exp_dir,varargin)
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

addpath(genpath('matlab_scripts'));
filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base_dir = fullfile(exp_dir,'individual_pictures');

image_dirs = dir(base_dir);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collect all the images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%registration will go in a step-wise fashion, register image 2 to image 1,
%then image 3 to 2, then image 4 to 3, ...
for i_num = 1:(size(image_dirs,1)-1)
    % for i_num = 1:
    file1 = fullfile(base_dir,image_dirs(i_num).name,filenames.puncta);
    file2 = fullfile(base_dir,image_dirs(i_num+1).name,filenames.puncta);
    
    img1 = int16(imread(file1));
    img2 = int16(imread(file2));
    
    img1_nii = make_nii(img1);
    img2_nii = make_nii(img2);
    
    file_nii_temp_1 = fullfile(base_dir,image_dirs(i_num).name,[filenames.puncta,'.nii']);
    file_nii_temp_2 = fullfile(base_dir,image_dirs(i_num+1).name,[filenames.puncta,'.nii']);
    
    save_nii(img1_nii,file_nii_temp_1);
    save_nii(img2_nii,file_nii_temp_2);
    
    reg_result_file = fullfile(base_dir,image_dirs(i_num+1).name,'reg_result.nii');
    affine_file = fullfile(base_dir,image_dirs(i_num+1).name,'reg_values.txt');
    
    system(['reg_aladin -source ', file_nii_temp_2, ' -target ',file_nii_temp_1, ...
        ' -result ',reg_result_file, ' -aff ',affine_file,' -noRot -rigOnly >/dev/null 2>/dev/null']);
    reg_result = load_nii(reg_result_file);
    
    imwrite(uint16(reg_result.img),file2);
    
    delete(file_nii_temp_1,file_nii_temp_2,reg_result_file);

    gel_image = imread(fullfile(base_dir,image_dirs(i_num+1).name,filenames.gel));
    nifti_transform = dlmread(affine_file,' ');
    transform = maketform('affine',[1,0;0,1;-1*nifti_transform(2,4),-1*nifti_transform(1,4)]);
    
    gel_image = imtransform(gel_image,transform, 'XData',[1 size(gel_image,2)], 'YData', [1 size(gel_image,1)]);
    imwrite(gel_image,fullfile(base_dir,image_dirs(i_num+1).name,filenames.gel));
    
    if (mod(i_num,round(size(image_dirs,1)/10)) == 0)
        fprintf('Done with %d/%d\n',i_num,size(image_dirs,1));
    end
end

%The output images are the same size as the input images, with the regions
%outside the image after registration filled in with zeros or near zeros.
%We want to determine where these filled in regions are and remove them
%from all the images to avoid edge effects. First we build a mask that only
%contains the image regions.
image_mask = [];
for i_num = 2:size(image_dirs,1)
    FA_file = fullfile(base_dir,image_dirs(i_num).name,filenames.puncta);
    
    FA_img = imread(FA_file);
    FA_img_binary = FA_img >= 5;
    if (size(image_mask,1) == 0)
        image_mask = FA_img_binary;
    else
        image_mask = FA_img_binary & image_mask;
    end
end
lims = [find(sum(image_mask,2),1,'first')+1,find(sum(image_mask,2),1,'last')-1;
    find(sum(image_mask),1,'first')+1,find(sum(image_mask),1,'last')-1];

%Now we apply the mask to the images in the set
for i_num = 1:size(image_dirs,1)
    FA_file = fullfile(base_dir,image_dirs(i_num).name,filenames.puncta);
    
    FA_img = imread(FA_file);
    
    FA_img = FA_img(lims(1,1):lims(1,2),lims(2,1):lims(2,2));
    
    imwrite(FA_img,FA_file);
end

toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
