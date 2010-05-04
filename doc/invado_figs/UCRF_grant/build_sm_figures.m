load('gel_image_ctrl_11_22_pos01.mat')
load('puncta_image_ctrl_11_22_pos01.mat')

addpath('../../../src_invado/visualize_cell_features/')

gel_images = all_gel{62442};
puncta_images = all_images{62442};

%sampling 8 images
sampled_gel = cell(0);
sampled_puncta = cell(0);
count = 1;
for i = round(linspace(74,118,8))
    sampled_gel{count} = gel_images{i};
    sampled_puncta{count} = puncta_images{i};
    count = count + 1;
end

gel_montage = build_montage_image_set(sampled_gel,'num_cols',length(sampled_gel));
puncta_montage = build_montage_image_set(sampled_puncta,'num_cols',length(sampled_puncta));

spacer = ones(1,size(gel_montage,2),3);

imwrite([puncta_montage;spacer;gel_montage], 'sm_8.png');

%sampling 10 images
sampled_gel = cell(0);
sampled_puncta = cell(0);
count = 1;
for i = round(linspace(74,118,10))
    sampled_gel{count} = gel_images{i};
    sampled_puncta{count} = puncta_images{i};
    count = count + 1;
end

gel_montage = build_montage_image_set(sampled_gel,'num_cols',length(sampled_gel));
puncta_montage = build_montage_image_set(sampled_puncta,'num_cols',length(sampled_puncta));

spacer = ones(1,size(gel_montage,2),3);

imwrite([puncta_montage;spacer;gel_montage], 'sm_10.png');

%sampling 15 images
sampled_gel = cell(0);
sampled_puncta = cell(0);
count = 1;
for i = round(linspace(74,118,15))
    sampled_gel{count} = gel_images{i};
    sampled_puncta{count} = puncta_images{i};
    count = count + 1;
end

gel_montage = build_montage_image_set(sampled_gel,'num_cols',length(sampled_gel));
puncta_montage = build_montage_image_set(sampled_puncta,'num_cols',length(sampled_puncta));

spacer = ones(1,size(gel_montage,2),3);

imwrite([puncta_montage;spacer;gel_montage], 'sm_15.png');

% load('gel_data.mat')
% load('puncta_data.mat')
% 
% addpath('../../../src_invado/visualize_cell_features/')
% 
% gel_images = all_gel{3446};
% puncta_images = all_images{3446};
% 
% %sampling 10 images
% sampled_gel = cell(0);
% sampled_puncta = cell(0);
% count = 1;
% for i = round(linspace(3,49,10))
%     sampled_gel{count} = gel_images{i};
%     sampled_puncta{count} = puncta_images{i};
%     count = count + 1;
% end
% 
% gel_montage = build_montage_image_set(sampled_gel,'num_cols',length(sampled_gel));
% puncta_montage = build_montage_image_set(sampled_puncta,'num_cols',length(sampled_puncta));
% 
% spacer = ones(1,size(gel_montage,2),3);
% 
% imwrite([puncta_montage;spacer;gel_montage], 'sm_10.png');
% 
% %sampling 15 images
% sampled_gel = cell(0);
% sampled_puncta = cell(0);
% count = 1;
% for i = round(linspace(3,49,15))
%     sampled_gel{count} = gel_images{i};
%     sampled_puncta{count} = puncta_images{i};
%     count = count + 1;
% end
% 
% gel_montage = build_montage_image_set(sampled_gel,'num_cols',length(sampled_gel));
% puncta_montage = build_montage_image_set(sampled_puncta,'num_cols',length(sampled_puncta));
% 
% spacer = ones(1,size(gel_montage,2),3);
% 
% imwrite([puncta_montage;spacer;gel_montage], 'sm_15.png');