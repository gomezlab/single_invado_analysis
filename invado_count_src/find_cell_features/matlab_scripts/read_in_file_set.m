function data_set = read_in_file_set(this_dir,filenames)
data_set = struct;

data_set.this_dir = this_dir;

data_set.gel_range = csvread(fullfile(this_dir, filenames.gel_range));

data_set.gel_image  = double(imread(fullfile(this_dir, filenames.gel)));
data_set.gel_image_norm = (data_set.gel_image - data_set.gel_range(1))/range(data_set.gel_range);
data_set.intensity_correction = csvread(fullfile(this_dir, filenames.intensity_correction));
data_set.gel_image_corr = data_set.gel_image*data_set.intensity_correction;

data_set.cell_mask = logical(imread(fullfile(this_dir, filenames.cell_mask)));
data_set.labeled_cells = imread(fullfile(this_dir, filenames.labeled_cell_mask));

data_set.puncta_range = csvread(fullfile(this_dir, filenames.puncta_range));

data_set.puncta_image  = double(imread(fullfile(this_dir, filenames.puncta)));
data_set.puncta_image_norm = (data_set.puncta_image - data_set.puncta_range(1))/range(data_set.puncta_range);

data_set.no_cells = imread(fullfile(this_dir, filenames.no_cells));

