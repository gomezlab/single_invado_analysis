function flat_field_correct_gel(exp_dir,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);
i_p.addOptional('debug',0,@(x)x == 1 | x == 0);
i_p.parse(exp_dir,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'));
filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base_dir = fullfile(i_p.Results.exp_dir,'individual_pictures');

image_dirs = dir(base_dir);
image_dirs = image_dirs(3:end);

gel_images = [];
for j=1:length(image_dirs)
    gel_image = imread(fullfile(base_dir,image_dirs(j).name,filenames.gel));
    if (any(size(gel_images) == 0))
        gel_images = gel_image;
    else
        gel_images = cat(3, gel_images,gel_image);
    end
end

gel_mean = mean(double(gel_images),3);
gel_mean = gel_mean - mean(gel_mean(:));

[x_coord,y_coord] = meshgrid(1:size(gel_mean,2),1:size(gel_mean,1));
fun = fit([x_coord(:),y_coord(:)],gel_mean(:),'poly22');
flat_fit_linear = feval(fun,[x_coord(:),y_coord(:)]);
plot(fun, [x_coord(:),y_coord(:)], gel_mean(:))

flat_fit = gel_mean;
flat_fit(1:length(flat_fit_linear)) = flat_fit_linear;

for j=1:length(image_dirs)
    gel_image = double(imread(fullfile(base_dir,image_dirs(j).name,filenames.gel)));
    gel_image = gel_image - flat_fit;
    imwrite(uint16(gel_image),fullfile(base_dir,image_dirs(j).name,filenames.gel),'Bitdepth',16);
end
toc;