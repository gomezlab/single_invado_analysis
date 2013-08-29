function cfg_struct = add_filenames_to_struct(cfg_struct)
%ADD_FILENAMES_TO_STRUCT    Takes in a struct, to which the filenames for
%                          various files made in the pipeline are added

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg_struct.objects = 'puncta_labeled.png';
cfg_struct.objects_perim = 'puncta_labeled_perim.png';
cfg_struct.objects_binary = 'puncta_binary.png';
cfg_struct.objects_highlight = 'puncta_highlight.png';

cfg_struct.objects_above_thresh = 'puncta_above_thresh.png';

cfg_struct.puncta = 'puncta.png';
cfg_struct.puncta_range = '../../puncta_props/image_processing/puncta_range.csv';
cfg_struct.puncta_vis_range = '../../puncta_props/image_processing/puncta_vis_range.csv';
cfg_struct.puncta_threshold = '../../puncta_props/image_processing/puncta_threshold.csv';
cfg_struct.puncta_threshold_plot = '../../puncta_props/image_processing/puncta_threshold.eps';

cfg_struct.cell_mask = 'cell_mask.png';
cfg_struct.cell_mask_highlight = 'highlighted_mask.png';
cfg_struct.no_cell_regions = '../../puncta_props/image_processing/no_cell_regions.png';

cfg_struct.gel = 'gel.png';
cfg_struct.gel_range = '../../puncta_props/image_processing/gel_range.csv';
cfg_struct.gel_vis_range = '../../puncta_props/image_processing/gel_vis_range.csv';
cfg_struct.bleaching_plot = '../../puncta_props/image_processing/bleaching_curves.eps';

cfg_struct.invado_data = '../../puncta_props/invado_analysis/invado_nums.csv';
cfg_struct.not_invado_data = '../../puncta_props/invado_analysis/not_invado_nums.csv';

cfg_struct.tracking = '../../tracking_matrices/tracking_seq.csv';
cfg_struct.tracking_matrix = '../../tracking_matrices/tracking_seq.csv';

cfg_struct.cell_props = '../../puncta_props/image_processing/single_puncta_props.mat';

cfg_struct.lineage_dir = '../../puncta_props/lin_time_series';
