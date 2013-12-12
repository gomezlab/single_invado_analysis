function data_set = read_in_file_set(this_dir,filenames)
data_set = struct;

data_set.this_dir = this_dir;

data_set.gel_vis_range = csvread(fullfile(this_dir, filenames.gel_vis_range));
data_set.gel_image  = double(imread(fullfile(this_dir, filenames.gel)));
data_set.gel_image_norm = (data_set.gel_image - data_set.gel_vis_range(1))/range(data_set.gel_vis_range);

data_set.objects = imread(fullfile(this_dir, filenames.objects));
data_set.objects_perim = imread(fullfile(this_dir, filenames.objects_perim));

data_set.puncta_vis_range = csvread(fullfile(this_dir, filenames.puncta_vis_range));
data_set.puncta_image  = double(imread(fullfile(this_dir, filenames.puncta)));
data_set.puncta_image_norm = (data_set.puncta_image - data_set.puncta_vis_range(1))/range(data_set.puncta_vis_range);

%if the no cells present mask exists, we read in that file, otherwise,
%assume that we are looking at fixed images and the cell mask should be
%used as the no cell region marker
if(exist(fullfile(this_dir, filenames.no_cell_regions), 'file'))
    data_set.no_cells = imread(fullfile(this_dir, filenames.no_cell_regions));
else
    data_set.no_cells = logical(imread(fullfile(this_dir, filenames.cell_mask)));
end

%read in the cell mask file if defined
if(exist(fullfile(this_dir, filenames.cell_mask), 'file'))
    data_set.cell_mask = logical(imread(fullfile(this_dir, filenames.cell_mask)));
end