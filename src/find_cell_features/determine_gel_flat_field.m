function determine_gel_flat_field(exp_dir,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
i_p = inputParser;

% i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);
i_p.addRequired('day_dir',@(x)exist(x,'dir') == 7);
i_p.addOptional('debug',0,@(x)x == 1 | x == 0);
i_p.parse(exp_dir,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'));
filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
field_dirs = dir(i_p.Results.day_dir);
field_dirs = field_dirs(3:end);

gel_images = [];
for i=1:length(field_dirs)
    base_dir = fullfile(i_p.Results.day_dir,field_dirs(i).name,'individual_pictures');
    
    image_dirs = dir(base_dir);
    image_dirs = image_dirs(3:end);
    
    for j=1:length(image_dirs)
        gel_image = imread(fullfile(base_dir,image_dirs(j).name,filenames.gel));
        if (any(size(gel_images) == 0))
            gel_images = gel_image;
        else
            gel_images = cat(3, gel_images,gel_image);
        end
    end
end

% gel_mean = mean(double(gel_images),3);
% gel_mean = gel_mean - mean(gel_mean(:));
% 
% [x_coord,y_coord] = meshgrid(1:size(gel_mean,2),1:size(gel_mean,1));
% fun = fit([x_coord(:),y_coord(:)],gel_mean(:),'poly22');
% flat_fit_linear = feval(fun,[x_coord(:),y_coord(:)]);
% plot(fun, [x_coord(:),y_coord(:)], gel_mean(:))
% 
% flat_fit = gel_mean;
% flat_fit(1:length(flat_fit_linear)) = flat_fit_linear;

gel_median = median(double(gel_images),3);

[x_coord,y_coord] = meshgrid(1:size(gel_median,2),1:size(gel_median,1));
fun = fit([x_coord(:),y_coord(:)],gel_median(:),'poly22');
flat_fit_linear = feval(fun,[x_coord(:),y_coord(:)]);
plot(fun, [x_coord(:),y_coord(:)], gel_median(:))

flat_fit = gel_median;
flat_fit(1:length(flat_fit_linear)) = flat_fit_linear;

imwrite(uint16(flat_fit),fullfile(i_p.Results.day_dir,'flat_field.png'),'Bitdepth',16);

toc;