function build_dual_highlight_image(base_dir,varargin)
% BUILD_DUAL_HIGHLIGHT_IMAGE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'build_dual_highlight_image';

i_p.addRequired('base_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(base_dir,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

binary_shift = logical(imread(fullfile(base_dir,'binary_shift.png')));

puncta_binary = logical(imread(fullfile(base_dir,'puncta_binary.png')));

gel_binary = logical(imread(fullfile(base_dir,'degradation_binary.png')));

puncta_image = double(imread(fullfile(base_dir,'registered_focal_image.png')));
puncta_image = puncta_image - min(puncta_image(binary_shift));
puncta_image = puncta_image .* (1/max(puncta_image(binary_shift)));
puncta_image(not(binary_shift)) = 0;

gel_image = double(imread(fullfile(base_dir,'registered_gel.png')));
gel_image = gel_image - min(gel_image(binary_shift));
gel_image = gel_image .* (1/max(gel_image(binary_shift)));
gel_image(not(binary_shift)) = 0;

labeled_all = zeros(size(binary_shift));
labeled_all(puncta_binary) = 1;
labeled_all(gel_binary) = 2;
labeled_all(gel_binary & puncta_binary) = 3;
c_map = [1 0 0; 0 1 0; 1 1 0];

addpath(genpath('..'))

gel_highlight = create_highlighted_image(gel_image,labeled_all,'color_map',c_map);
puncta_highlight = create_highlighted_image(puncta_image,labeled_all,'color_map',c_map);

both = [gel_highlight, 0.5*ones(size(gel_highlight,1),1,3), puncta_highlight];

imwrite(both, fullfile(base_dir,'dual_highlight.png'));

