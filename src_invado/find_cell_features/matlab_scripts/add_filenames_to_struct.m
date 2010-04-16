function cfg_struct = add_filenames_to_struct(cfg_struct)
%AD_FILENAMES_TO_STRUCT    Takes in a struct, to which the filenames for
%                          various files made in the pipeline are added

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg_struct.adhesions_filename = 'puncta_labeled.png';
cfg_struct.puncta_filename = 'registered_focal_image.png';
cfg_struct.gel_filename = 'registered_gel.png';
cfg_struct.binary_shift_filename = 'binary_shift.png';
cfg_struct.cell_mask_filename = 'cell_mask.png';
cfg_struct.intensity_correction_filename = 'intensity_correction.csv';
