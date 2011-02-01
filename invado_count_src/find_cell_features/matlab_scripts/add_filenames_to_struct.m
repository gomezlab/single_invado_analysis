function filenames = add_filenames_to_struct(filenames)
%ADD_FILENAMES_TO_STRUCT    Takes in a struct, to which the filenames for
%                          various files made in the pipeline are added

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filenames.puncta_labeled_filename = 'puncta_labeled.png';
filenames.puncta_filename = 'focal_image.png';
filenames.puncta_median = '../../../all_field_resources/puncta_median.png';
filenames.puncta_range_file = '../../../all_field_resources/puncta_image_range.csv';

filenames.cell_mask_filename = 'cell_mask.png';
filenames.highlighted_cell_mask_filename = 'highlighted_mask.png';
filenames.labeled_cell_mask_filename = 'cell_mask_labeled.png';
filenames.labeled_cell_mask_perim = 'cell_mask_perim_labeled.png';

filenames.gel_filename = 'gel.png';
filenames.gel_range_file = '../../../all_field_resources/gel_image_range.csv';

filenames.intensity_correction_filename = 'intensity_correction.csv';

filenames.no_cells = '../../adhesion_props/no_cell_regions.png';

filenames.invader_vis = 'invader_vis.png';

filenames.tracking_raw = '../../adhesion_props/raw_tracking_data.mat';
filenames.tracking = '../../tracking_matrices/tracking_seq.csv';

filenames.cell_props = '../../adhesion_props/cell_props.mat';

filenames.lineage_dir = '../../adhesion_props/lin_time_series/';