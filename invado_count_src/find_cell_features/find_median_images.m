function find_median_images(base_dir,varargin)
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('base_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(base_dir,varargin{:});

addpath('matlab_scripts');
filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fields = dir(base_dir);
fields = filter_to_time_series(fields);

total_images = length(dir(fullfile(base_dir,fields(1).name,'individual_pictures'))) - 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Puncta Image Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
images_puncta = cell(0);
for i = 1:total_images
    for j=1:length(fields)
        image_base = fullfile(base_dir,fields(j).name,'individual_pictures');
        image_nums = dir(image_base);
        image_nums = image_nums(3:end);
        image_file = fullfile(image_base,image_nums(i).name, filenames.puncta);
                
        images_puncta{length(images_puncta)+1} = imread(image_file);
    end
        
    if (mod(i,10) == 0)
        disp(['Done reading: ',num2str(i), '/', num2str(total_images)])
    end
end

images = zeros([size(images_puncta{1}),length(images_puncta)]);

for i_num = 1:length(images)
    images(:,:,i_num) = images_puncta{i_num};
end

median_image = median(images,3);

output_base = fullfile(base_dir,fields(1).name,'individual_pictures');
image_nums = dir(output_base);
image_nums = image_nums(3:end);
output_file = fullfile(output_base,image_nums(1).name, filenames.puncta_median);

[output_path] = fileparts(output_file);

if (not(exist(fullfile(output_path),'dir')))
    mkdir(fullfile(output_path));
end

imwrite(uint16(median_image),output_file,'bitdepth',16);

clear images_puncta images;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gel Image Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
images_gel = cell(0);
for i = 1:total_images
    for j=1:length(fields)
        image_base = fullfile(base_dir,fields(j).name,'individual_pictures');
        image_nums = dir(image_base);
        image_nums = image_nums(3:end);
        image_file = fullfile(image_base,image_nums(i).name, filenames.gel);
                
        images_gel{length(images_gel)+1} = imread(image_file);
    end
        
    if (mod(i,10) == 0)
        disp(['Done reading: ',num2str(i), '/', num2str(total_images)])
    end
end

images = zeros([size(images_gel{1}),length(images_gel)]);

for i_num = 1:length(images)
    images(:,:,i_num) = images_gel{i_num};
end

median_image = median(images,3);

output_base = fullfile(base_dir,fields(1).name,'individual_pictures');
image_nums = dir(output_base);
image_nums = image_nums(3:end);
output_file = fullfile(output_base,image_nums(1).name, filenames.gel_median);

[output_path] = fileparts(output_file);

if (not(exist(fullfile(output_path),'dir')))
    mkdir(fullfile(output_path));
end

imwrite(uint16(median_image),output_file,'bitdepth',16);

clear images_gel images;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nuclei Image Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% images_nuc = cell(0);
% for i = 1:total_images
%     for j=1:length(fields)
%         image_base = fullfile(base_dir,fields(j).name,'individual_pictures');
%         image_nums = dir(image_base);
%         image_nums = image_nums(3:end);
%         image_file = fullfile(image_base,image_nums(i).name, filenames.nucleus);
%                 
%         images_nuc{length(images_nuc)+1} = imread(image_file);
%     end
%         
%     if (mod(i,10) == 0)
%         disp(['Done reading: ',num2str(i), '/', num2str(total_images)])
%     end
% end
% 
% images = zeros([size(images_nuc{1}),length(images_nuc)]);
% 
% for i_num = 1:length(images)
%     images(:,:,i_num) = images_nuc{i_num};
% end
% 
% median_image = median(images,3);
% 
% output_base = fullfile(base_dir,fields(1).name,'individual_pictures');
% image_nums = dir(output_base);
% image_nums = image_nums(3:end);
% output_file = fullfile(output_base,image_nums(1).name, filenames.nucleus_median);
% 
% [output_path] = fileparts(output_file);
% 
% if (not(exist(fullfile(output_path),'dir')))
%     mkdir(fullfile(output_path));
% end
% 
% imwrite(uint16(median_image),output_file,'bitdepth',16);
% 
% clear images_nuc images;

toc;