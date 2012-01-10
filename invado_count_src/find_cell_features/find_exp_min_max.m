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

gel_image_pix_vals = [];
puncta_image_pix_vals = [];

for i = 1:total_images
    for j=1:length(fields)
        image_base = fullfile(base_dir,fields(j).name,'individual_pictures');
        image_nums = dir(image_base);
        image_nums = image_nums(3:end);
                
        gel_image = imread(fullfile(image_base,image_nums(i).name,filenames.gel));
        gel_image_pix_vals = [gel_image_pix_vals; gel_image(:)]; %#ok<AGROW>
        
        puncta_image = imread(fullfile(image_base,image_nums(i).name,filenames.puncta));
        puncta_image_pix_vals = [puncta_image_pix_vals; puncta_image(:)]; %#ok<AGROW>
    end

    if (mod(i,10) == 0)
        disp(['Done reading: ',num2str(i)]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Range determination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

puncta_image_pix_vals = sort(puncta_image_pix_vals);
gel_image_pix_vals = sort(gel_image_pix_vals);

end_trim_amount = round(size(puncta_image_pix_vals,1)*0.00005);

puncta_image_range = [punta_image_pix_vals(1),punta_image_pix_vals(end); ...
    punta_image_pix_vals(1),punta_image_pix_vals(end-end_trim_amount)];

gel_image_range = [gel_image_pix_vals(1),gel_image_pix_vals(end); ...
    gel_image_pix_vals(1),gel_image_pix_vals(end-end_trim_amount)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_base = fullfile(base_dir,fields(1).name,'individual_pictures');
image_nums = dir(output_base);
image_nums = image_nums(3:end);

gel_output_file = fullfile(output_base,image_nums(1).name, filenames.gel_range);
puncta_output_file = fullfile(output_base,image_nums(1).name, filenames.puncta_range);

[output_path] = fileparts(gel_output_file);

if (not(exist(fullfile(output_path),'dir')))
    mkdir(fullfile(output_path));
end

csvwrite(gel_output_file,gel_image_range);
csvwrite(puncta_output_file,puncta_image_range);


