function find_median_images(base_dir,target_file,varargin)
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('base_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('target_file','focal_image.png',@(x)ischar(x)||isarray(x));
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(base_dir,varargin{:});

target_file = i_p.Results.target_file;

addpath('matlab_scripts');
filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fields = dir(base_dir);
fields = filter_to_time_series(fields);

total_images = length(dir(fullfile(base_dir,fields(1).name,'individual_pictures'))) - 2;

images_cell = cell(0);

for i = 1:total_images
    for j=1:length(fields)
        image_base = fullfile(base_dir,fields(j).name,'individual_pictures');
        image_nums = dir(image_base);
        image_nums = image_nums(3:end);
        image_file = fullfile(image_base,image_nums(i).name, target_file);
                
        images_cell{length(images_cell)+1} = imread(image_file);
    end
        
    if (mod(i,10) == 0)
        disp(['Done reading: ',num2str(i)])
    end
end

images = zeros([size(images_cell{1}),length(images_cell)]);

for i_num = 1:length(images)
    images(:,:,i_num) = images_cell{i_num};
    if (mod(i,10) == 0)
        disp(['Done re-ordering: ',num2str(i)])
    end
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
