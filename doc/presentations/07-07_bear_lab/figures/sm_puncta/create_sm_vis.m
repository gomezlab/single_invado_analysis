clear;

load('gel_sm_set.mat');
load('puncta_sm_set.mat');
load('gel_sm_set_no_high.mat');
load('puncta_sm_set_no_high.mat');

addpath('../../../../../src/visualize_cell_features/');

num_images = 10;

sets = cell(0);
sets{1} = all_gel_temp(round(linspace(5,length(all_gel_temp),num_images)));
sets{2} = all_puncta_temp(round(linspace(5,length(all_gel_temp),num_images)));
sets{3} = all_gel_temp_no_high(round(linspace(5,length(all_gel_temp),num_images)));
sets{4} = all_puncta_temp_no_high(round(linspace(5,length(all_gel_temp),num_images)));

montage_sets = cell(0);
for i = 1:length(sets)
    montage_sets{i} = build_montage_image_set(sets{i},'num_cols',num_images);
end

spacer = ones(1,size(montage_sets{1},2),3);

composite_image = [montage_sets{4};spacer;montage_sets{2}; ... 
    spacer;spacer;spacer;montage_sets{1};spacer;montage_sets{3}];

imwrite(composite_image,'montage_composite.png')