function register_with_matlab(exp_dir,varargin)
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
% Determine Registration in Gel Images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%write a dummy registration values file, to make sure follow up file count checks work
affine_file = fullfile(base_dir,image_dirs(1).name,'reg_values.txt');
csvwrite(affine_file,diag(ones(3,1)));

[optimizer, metric] = imregconfig('monomodal');

reg_start = tic;
%registration will go in a step-wise fashion, register image 2 to image 1,
%then image 3 to 2, then image 4 to 3, ...
for i_num = 1:(size(image_dirs,1)-1)
    affine_file = fullfile(base_dir,image_dirs(i_num+1).name,'reg_values.txt');
    if (exist(affine_file,'file'))
        continue;
    end
    
    file1 = fullfile(base_dir,image_dirs(i_num).name,filenames.gel);
    file2 = fullfile(base_dir,image_dirs(i_num+1).name,filenames.gel);
    
    img1 = imread(file1);
    img2 = imread(file2);
    
    reg_info = imregtform(img2,img1,'translation',optimizer,metric);
    csvwrite(affine_file,reg_info.T);
    
    if (mod(i_num,round(size(image_dirs,1)/10)) == 0)
        run_so_far = toc(reg_start);
        est_left = (run_so_far/i_num)*(size(image_dirs,1) - i_num);
        fprintf('Done Calculating %d/%d (%d min left)\n',i_num,size(image_dirs,1),round(est_left/60));
    end
end

%load set of translation movements
left_right_trans_set = zeros(size(image_dirs,1),1);
up_down_trans_set = zeros(size(image_dirs,1),1);
for i_num = 2:size(image_dirs,1)
    affine_file = fullfile(base_dir,image_dirs(i_num).name,'reg_values.txt');
    affine_mat = csvread(affine_file);
    left_right_trans_set(i_num) = affine_mat(3,1);
    up_down_trans_set(i_num) = affine_mat(3,2);
end

plot(cumsum(left_right_trans_set),cumsum(up_down_trans_set));
xlabel('Left-Right Translation', 'Fontsize',16)
ylabel('Up-Down Translation', 'Fontsize',16);

output_file = fullfile(base_dir,image_dirs(1).name,'../../puncta_props/translation.pdf');
[output_folder, ~, ~] = fileparts(output_file);
if (not(exist(output_folder,'dir')))
    mkdir(fileparts(output_file));
end
print('-dpdf', output_file);
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply Registration to All Images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i_num = 2:size(image_dirs,1)
    affine_mat = [1,0;0,1; ...
        sum(left_right_trans_set(1:i_num)),sum(up_down_trans_set(1:i_num))];
    
    transform = maketform('affine',affine_mat);
    
    gel_image = imread(fullfile(base_dir,image_dirs(i_num).name,filenames.gel));
    puncta_image = imread(fullfile(base_dir,image_dirs(i_num).name,filenames.puncta));
    
    gel_image = imtransform(gel_image,transform,'XData',[1 size(puncta_image,2)],'YData',[1 size(puncta_image,1)]);
    puncta_image = imtransform(puncta_image,transform,'XData',[1 size(puncta_image,2)],'YData',[1 size(puncta_image,1)]);
    
    imwrite(gel_image,fullfile(base_dir,image_dirs(i_num).name,filenames.gel));
    imwrite(puncta_image,fullfile(base_dir,image_dirs(i_num).name,filenames.puncta));
    
    if (mod(i_num,round(size(image_dirs,1)/10)) == 0)
        fprintf('Done Applying %d/%d\n',i_num,size(image_dirs,1));
    end
end

left_right_history = cumsum(left_right_trans_set);
up_down_history = cumsum(up_down_trans_set);

sample_image = imread(fullfile(base_dir,image_dirs(1).name,filenames.gel));
lims = [1 + ceil(max(up_down_history)), size(sample_image,1) + floor(min(up_down_history)); ...
    1 + ceil(max(left_right_history)), size(sample_image,2) + floor(min(left_right_history))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove the Region Outside the Registered Region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i_num = 1:size(image_dirs,1)
    gel_file = fullfile(base_dir,image_dirs(i_num).name,filenames.gel);
    gel_img = imread(gel_file);
    gel_img = gel_img(lims(1,1):lims(1,2),lims(2,1):lims(2,2));
    imwrite(gel_img,gel_file);
    
    puncta_file = fullfile(base_dir,image_dirs(i_num).name,filenames.puncta);
    puncta_img = imread(puncta_file);
    puncta_img = puncta_img(lims(1,1):lims(1,2),lims(2,1):lims(2,2));
    imwrite(puncta_img,puncta_file);
end

toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
