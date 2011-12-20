function build_raw_montage_file_set(base_dir,target_file,varargin)
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('base_dir',@(x)exist(x,'dir') == 7);
i_p.addRequired('target_file',@(x)ischar(x)||isarray(x));
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(base_dir,target_file,varargin{:});

addpath(genpath('..'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fields = dir(base_dir);
fields = filter_to_time_series(fields);

images_composite = cell(0);
base_image_size = [NaN, NaN];

total_images = length(dir(fullfile(base_dir,fields(1).name,'individual_pictures'))) - 2;
all_pix_vals = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Composite Images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%build up composite images by scan through the image number, so time step, then
%each field
for i = 1:total_images
    images = cell(0);
    
    for j=1:length(fields)
        image_base = fullfile(base_dir,fields(j).name,'individual_pictures');
        image_nums = dir(image_base);
        image_nums = image_nums(3:end);
        image_file = fullfile(image_base,image_nums(i).name, target_file);
        
        %deal with the case where the image does not exist
        if(not(exist(image_file,'file')))
            continue;
        end  
        images{j} = imread(image_file);
        base_image_size = size(images{j});
    end
    
    %scan back through the images collected filling in missing files with
    %Nan blanks
    for j=1:length(fields)
        if(length(images) < j || isempty(images{j}))
            images{j} = NaN*ones(base_image_size);
        end
    end
    
    %reverse every other column of images, to deal with the zig zag pattern
    image_temp = images;
    for j=1:(length(images)/5)
        if (mod(j,2) == 1)
            continue;
        end
        for k=1:5
            image_temp{j*5+k} = images{j*5+(6-k)};
        end
    end
    images = image_temp;
    
	%build up composite image, cat first 5 images in one row
    images_composite{i} = cat(1,images{1:5});
    for j=2:(length(images)/5)
		%then cat each subsequent five images in the next row
        images_composite{i} = cat(2,images_composite{i}, cat(1,images{((j-1)*5+1):(j*5)}));
    end
    
	%finally resize to a reasonable image size
    images_composite{i} = imresize(images_composite{i},[800,NaN]);
    all_pix_vals = [all_pix_vals; images_composite{i}(:)]; %#ok<AGROW>
    
    if (mod(i,10) == 0)
        disp(['Done processing: ',num2str(i), ' of ', target_file])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalize Composite Images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

all_pix_vals = sort(all_pix_vals);
%we only need to do image normalization if the images are not already in
%color (i.e. have more than 2 dimensions)
if (ndims(images_composite{1}) < 3)
    %toss out the top 0.05% of pixels for the image min/max creation
    end_trim_amount = round(length(all_pix_vals)*0.0005);
    images_min_max = double([all_pix_vals(1), all_pix_vals(end-end_trim_amount)]);
    for i=1:length(images_composite)
		%scale and then truncate the composite images
        images_composite{i} = (double(images_composite{i}) - images_min_max(1))/range(images_min_max);
        images_composite{i}(images_composite{i} > images_min_max(2)) = images_min_max(2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output Composite Images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[pathstr, name] = fileparts(target_file); %#ok<ASGLU>
output_dir = fullfile(base_dir,'montage',name);

if (not(exist(fullfile(output_dir),'dir')))
    mkdir(fullfile(output_dir));
end

for i=1:length(images_composite)
    padding_length = num2str(length(num2str(total_images)));
    padded_i_num = sprintf(['%0',padding_length,'d'],i);
    
    imwrite(images_composite{i}, fullfile(output_dir,[padded_i_num,'.png']));
    if (mod(i,10) == 0)
        disp(['Done writing: ',num2str(i),' of ', target_file])
    end
end
toc;
