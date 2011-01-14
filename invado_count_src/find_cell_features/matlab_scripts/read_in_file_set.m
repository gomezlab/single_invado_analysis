function data_set = read_in_file_set(this_dir,filenames)
data_set = struct;

%read in and normalize the input gel image
data_set.gel_image  = imread(fullfile(this_dir, filenames.gel_filename));
data_set.gel_no_norm = double(data_set.gel_image);
scale_factor = double(intmax(class(data_set.gel_image)));
data_set.gel_image  = double(data_set.gel_image)/scale_factor;

%read in gel min/max
data_set.gel_range = csvread(fullfile(this_dir, filenames.gel_range_file));

%read in the intensity correction coefficient
data_set.intensity_correction = csvread(fullfile(this_dir, filenames.intensity_correction_filename));

%read in the cell mask file
data_set.cell_mask = logical(imread(fullfile(this_dir, filenames.cell_mask_filename)));

%read in the labeled cells
data_set.labeled_cells = imread(fullfile(this_dir, filenames.labeled_cell_mask_filename));
