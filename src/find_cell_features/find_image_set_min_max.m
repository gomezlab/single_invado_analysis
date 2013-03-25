function find_image_set_min_max(I_folder, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'FIND_IMAGE_SET_MIN_MAX';

i_p.addRequired('I_folder',@(x)exist(x,'dir') == 7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(I_folder,varargin{:});

addpath('matlab_scripts');
filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base_dir = fullfile(i_p.Results.I_folder,'individual_pictures');
image_dirs = dir(base_dir);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);

gel_range = [Inf -Inf];
puncta_range = [Inf -Inf];

gel_data = [];
puncta_data = [];

for i_num = 1:size(image_dirs,1)
    gel_image = imread(fullfile(base_dir,image_dirs(i_num).name,filenames.gel));
    if (any(size(gel_data)) == 0)
        gel_data = zeros(size(gel_image,1),size(gel_image,2),size(image_dirs,1));
    end
    gel_data(:,:,i_num) = gel_image; %#ok<AGROW>
    gel_range = update_min_max(gel_image(:),gel_range);

    puncta_image = imread(fullfile(base_dir,image_dirs(i_num).name,filenames.puncta));
    if (any(size(puncta_data)) == 0)
        puncta_data = zeros(size(puncta_image,1),size(puncta_image,2),size(image_dirs,1));
    end
    puncta_data(:,:,i_num) = puncta_image; %#ok<AGROW>
    puncta_range = update_min_max(puncta_image(:),puncta_range);
    
    if (mod(i_num,10)==0)
        disp(['Done with ',num2str(i_num),'/',num2str(size(image_dirs,1))])
    end
end

output_file = fullfile(base_dir,image_dirs(i_num).name,filenames.gel_vis_range);
if(not(exist(fileparts(output_file),'dir')))
    mkdir(fileparts(output_file))
end
gel_vis_range = quantile(double(gel_data(:)),[0.001,0.999]);
csvwrite(output_file,gel_vis_range)
csvwrite(fullfile(base_dir,image_dirs(i_num).name,filenames.gel_range),gel_range);

% puncta_vis_range = quantile(double(puncta_data(:)),[0.001,0.999]);
csvwrite(fullfile(base_dir,image_dirs(i_num).name,filenames.puncta_vis_range),puncta_range)
csvwrite(fullfile(base_dir,image_dirs(i_num).name,filenames.puncta_range),puncta_range)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function current_min_max = update_min_max(image_vals,current_min_max)

this_min_max = [min(image_vals(:)),max(image_vals(:))];
if (current_min_max(1) > this_min_max(1))
    current_min_max(1) = this_min_max(1);
end
if (current_min_max(2) < this_min_max(2))
    current_min_max(2) = this_min_max(2);
end

end