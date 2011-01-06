function cfg_struct = add_filenames_to_struct(cfg_struct)
%ADD_FILENAMES_TO_STRUCT    Takes in a struct, to which the filenames for
%                          various files made in the pipeline are added

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg_struct.puncta_labeled_filename = 'puncta_labeled.png';
cfg_struct.puncta_filename = 'focal_image.png';
cfg_struct.puncta_median = '../../../all_field_resources/puncta_median.png';
cfg_struct.puncta_range_file = '../../../all_field_resources/puncta_image_range.csv';

cfg_struct.cell_mask_filename = 'cell_mask.png';
cfg_struct.highlighted_cell_mask_filename = 'highlighted_mask.png';
cfg_struct.labeled_cell_mask_filename = 'cell_mask_labeled.png';

cfg_struct.gel_filename = 'gel.png';
cfg_struct.gel_range_file = '../../../all_field_resources/gel_image_range.csv';

cfg_struct.intensity_correction_filename = 'intensity_correction.csv';