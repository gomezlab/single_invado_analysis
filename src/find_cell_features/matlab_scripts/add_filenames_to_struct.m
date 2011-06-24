function cfg_struct = add_filenames_to_struct(cfg_struct)
%AD_FILENAMES_TO_STRUCT    Takes in a struct, to which the filenames for
%                          various files made in the pipeline are added

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg_struct.affine_matrix = 'affine_matrix.csv';
cfg_struct.registration_diff = 'registration_diffs.csv';
cfg_struct.affine_matrix_cascade = 'affine_matrix_cascade.csv';

cfg_struct.objects = 'puncta_labeled.png';
cfg_struct.objects_perim = 'puncta_labeled_perim.png';
cfg_struct.objects_binary = 'puncta_binary.png';
cfg_struct.objects_highlight = 'puncta_highlight.png';

cfg_struct.puncta = 'puncta.png';
cfg_struct.puncta_unreg = 'unregistered_puncta_image.png';
cfg_struct.puncta_range = '../../puncta_props/puncta_image_range.csv';
cfg_struct.puncta_threshold = '../../puncta_props/puncta_threshold.csv';
cfg_struct.puncta_threshold_plot = '../../puncta_props/puncta_threshold.eps';

cfg_struct.gel = 'gel.png';
cfg_struct.gel_unreg = 'unregistered_gel.png';
cfg_struct.gel_range = '../../puncta_props/gel_image_range.csv';
cfg_struct.bleaching_plot = '../../puncta_props/bleaching_curves.eps';


cfg_struct.binary_shift = 'binary_shift.png';

cfg_struct.cell_mask = 'cell_mask.png';
cfg_struct.cell_mask_highlight = 'highlighted_mask.png';

cfg_struct.intensity_correction = 'intensity_correction.csv';

cfg_struct.no_cell_regions = '../../puncta_props/no_cell_regions.png';