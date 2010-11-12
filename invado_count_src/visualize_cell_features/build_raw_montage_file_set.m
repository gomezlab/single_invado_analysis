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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fields = dir(base_dir);
%toss out first two entries they are '.' and '..', do same in below
%searches
fields = fields(3:end);
if (strmatch(fields(1).name,'montage'))
    fields = fields(2:end);
end
if (strmatch(fields(1).name,'overall_results'))
    fields = fields(2:end);
end

images_composite = cell(0);
base_image_size = [NaN, NaN];

total_images = length(dir(fullfile(base_dir,fields(1).name,'individual_pictures'))) - 2;

for i = 1:total_images
    images = cell(0);
    
    for j=1:length(fields)
        image_base = fullfile(base_dir,fields(j).name,'individual_pictures');
        image_nums = dir(image_base);
        image_nums = image_nums(3:end);
        image_file = fullfile(image_base,image_nums(i).name, target_file);
        
         %deal with the case where the image does not exist
        if(not(exist(image_file,'file')))
            images{j} = NaN*ones(base_image_size);
            continue;
        end
        
        images{j} = imread(image_file);
        base_image_size = size(images{j});
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
    
    images_composite{i} = cat(1,images{1:5});
    
    for j=2:(length(images)/5)
        images_composite{i} = cat(2,images_composite{i}, cat(1,images{((j-1)*5+1):(j*5)}));
    end
    
    images_composite{i} = imresize(images_composite{i},[800,NaN]);
        
    if (mod(i,10) == 0)
        disp(['Done processing: ',num2str(i)])
    end
end

%we only need to do image normalization if the images are not already in
%color (i.e. have more than 2 dimensions)
if (length(size(images_composite{1})) < 3)
    images_min_max = double([min(min([images_composite{1:end}])),max(max([images_composite{1:end}]))]);
    
    for i=1:length(images_composite)
        images_composite{i} = (double(images_composite{i}) - images_min_max(1))/range(images_min_max);
    end
end


[pathstr, name] = fileparts(target_file); %#ok<ASGLU>
output_dir = fullfile(base_dir,'montage',name);

if (not(exist(fullfile(output_dir),'dir')))
    mkdir(fullfile(output_dir));
end

for i=1:length(images_composite)
    padding_length = num2str(length(num2str(total_images)));
    padded_i_num = sprintf(['%0',padding_length,'d'],i);
    
    imwrite(images_composite{i}, fullfile(output_dir,[padded_i_num,'.png']));
    %     imwrite(gel_composite{i}, fullfile(output_dir,'gel',[padded_i_num,'.png']));
    %     imwrite(cat(3,gel_composite{i},images_composite{i},zeros(size(images_composite{i}))), fullfile(output_dir,'both',[padded_i_num,'.png']));
    
    if (mod(i,10) == 0)
        disp(['Done writing: ',num2str(i)])
    end
end

if (i_p.Results.debug)
    toc;
end