function cfg_struct = add_filenames_to_struct(cfg_struct)
%AD_FILENAMES_TO_STRUCT    Takes in a struct, to which the filenames for
%                          various files made in the pipeline are added

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg_struct.affine_matrix = 'affine_matrix.csv';
cfg_struct.registration_diff = 'registration_diffs.csv';
cfg_struct.affine_matrix_cascade = 'affine_matrix_cascade.csv';

cfg_struct.adhesions_filename = 'puncta_labeled.png';

cfg_struct.puncta = 'focal_image.png';
cfg_struct.puncta_unreg = 'unregistered_focal_image.png';
cfg_struct.puncta_range = '../../adhesion_props/puncta_image_range.csv';

cfg_struct.gel = 'gel.png';
cfg_struct.gel_unreg = 'unregistered_gel.png';
cfg_struct.gel_range = '../../adhesion_props/gel_image_range.csv';

cfg_struct.binary_shift = 'binary_shift.png';

cfg_struct.cell_mask_filename = 'cell_mask.png';
cfg_struct.intensity_correction_filename = 'intensity_correction.csv';