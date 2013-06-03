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
flat_field_file = fullfile(i_p.Results.exp_dir,'puncta_props','flat_field.png');
if (exist(flat_field_file,'file'))
    flat_field = double(imread(flat_field_file));
else
    return;
end

flat_field = flat_field - mean(flat_field(:));

base_dir = fullfile(i_p.Results.exp_dir,'individual_pictures');

image_dirs = dir(base_dir);
image_dirs = image_dirs(3:end);

for j=1:length(image_dirs)
    gel_image = double(imread(fullfile(base_dir,image_dirs(j).name,filenames.gel)));
    gel_image = gel_image - flat_field;
    imwrite(uint16(gel_image),fullfile(base_dir,image_dirs(j).name,filenames.gel),'Bitdepth',16);
end
toc;