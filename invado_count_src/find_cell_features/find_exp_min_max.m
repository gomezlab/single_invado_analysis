function find_exp_min_max(base_dir,varargin)
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

gel_image_range = [Inf -Inf];
puncta_image_range = [Inf -Inf];

for i = 1:total_images
    for j=1:length(fields)
        image_base = fullfile(base_dir,fields(j).name,'individual_pictures');
        image_nums = dir(image_base);
        image_nums = image_nums(3:end);
                
        gel_image = imread(fullfile(image_base,image_nums(i).name,filenames.gel_filename));
        this_gel_min_max = [min(gel_image(:)), max(gel_image(:))];
        if (this_gel_min_max(1) < gel_image_range(1))
            gel_image_range(1) = this_gel_min_max(1);
        end
        if (this_gel_min_max(2) > gel_image_range(2))
            gel_image_range(2) = this_gel_min_max(2);
        end
        
        puncta_image = imread(fullfile(image_base,image_nums(i).name,filenames.puncta_filename));
        this_puncta_min_max = [min(puncta_image(:)), max(puncta_image(:))];
        if (this_puncta_min_max(1) < puncta_image_range(1))
            puncta_image_range(1) = this_puncta_min_max(1);
        end
        if (this_puncta_min_max(2) > puncta_image_range(2))
            puncta_image_range(2) = this_puncta_min_max(2);
        end

    end
        
    if (mod(i,10) == 0)
        disp(['Done reading: ',num2str(i)])
    end
end

output_base = fullfile(base_dir,fields(1).name,'individual_pictures');
image_nums = dir(output_base);
image_nums = image_nums(3:end);

gel_output_file = fullfile(output_base,image_nums(1).name, filenames.gel_range_file);
puncta_output_file = fullfile(output_base,image_nums(1).name, filenames.puncta_range_file);

[output_path] = fileparts(gel_output_file);

if (not(exist(fullfile(output_path),'dir')))
    mkdir(fullfile(output_path));
end

csvwrite(gel_output_file,gel_image_range);
csvwrite(puncta_output_file,puncta_image_range);