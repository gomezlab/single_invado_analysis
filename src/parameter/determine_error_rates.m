function determine_error_rates(seg_folder,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_time = tic;

if (not(exist(seg_folder,'dir')))
    error(['Not a folder: ',seg_folder])
end

addpath('../find_cell_features/matlab_scripts/');
addpath('../visualize_cell_features/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Identification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp = dir(seg_folder);
assert(strcmp(temp(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(temp(2).name, '..'), 'Error: expected ".." to be second string in the dir command')

temp = temp(3:end);

to_seg_found = 0;
files = cell(0);
for i=1:length(temp)
    if (regexpi(temp(i).name,'Epi_FITC'))
        source_file = temp(i).name;
    elseif (regexpi(temp(i).name,'to_seg'))
        to_seg_found = 1;
    elseif (regexpi(temp(i).name,'stats'))
    elseif (regexpi(temp(i).name,'consensus.png'))
    elseif (regexpi(temp(i).name,'seg_results'))
    else
        files{length(files)+1} = temp(i).name;
    end
end

if (not(exist('source_file','var')))
    warn('Can not find source file');
end
if (not(to_seg_found))
    warn('Can not find to seg file');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hand Segmentation Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_image = imread(fullfile(seg_folder,files{1}));
hand_seg_images = zeros(size(temp_image,1),size(temp_image,2),length(files));
for i=1:length(files)
    hand_seg_images(:,:,i) = process_hand_segmentation(fullfile(seg_folder,files{i}));
end

consensus = sum(hand_seg_images,3);
consensus = consensus >= 0.5*size(hand_seg_images,3);
consensus_label = bwlabel(consensus, 8);
consensus_props = regionprops(consensus_label,'all');
save(fullfile(seg_folder,'consensus_stats.mat'),'consensus_props');
imwrite(consensus,fullfile(seg_folder,'consensus.png'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Threshold Determination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source Image Filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source_i = double(imread(fullfile(seg_folder,source_file)));
source_i_norm = (source_i - min(source_i(:)))/range(source_i(:));

I_filt = fspecial('disk',11);
blurred_image = imfilter(source_i,I_filt,'same',mean(source_i(:)));
high_passed = source_i - blurred_image;

puncta_mean = mean(double(high_passed(:)));
puncta_std = std(double(high_passed(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seed + Size thresholding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% std_ranges = 0.5:0.5:12;
% min_size_vals = 1:10;
% 
% matches = zeros(length(std_ranges),length(min_size_vals));
% comp_hits = zeros(length(std_ranges),length(min_size_vals));
% % Puncta marking statistics
% for size_index = 1:length(min_size_vals)
%     for t2_index = 1:length(std_ranges)
%         filter_thresh = puncta_mean + puncta_std*std_ranges(t2_index);
%         threshed_image = find_threshed_image(high_passed,filter_thresh);
%         threshed_image_label = bwlabel(threshed_image,4);
%         threshed_image_label = filter_label_mat_by_size(threshed_image_label,min_size_vals(size_index));
%         threshed_image = threshed_image_label > 0;
%         
% %         %segmentation results picture
% %         if (not(exist(fullfile(seg_folder,'seg_results'),'dir')))
% %             mkdir(fullfile(seg_folder,'seg_results'))
% %         end
% %         temp = zeros(size(threshed_image));
% %         temp(threshed_image) = 1;
% %         temp(consensus) = 2;
% %         temp(threshed_image & consensus) = 3;
% %         thresh_highlight = create_highlighted_image(source_i_norm,temp,'color_map',jet(3));
% %         imwrite(thresh_highlight, ... 
% %             fullfile(seg_folder,'seg_results',[num2str(std_ranges(t2_index)),'_',num2str(size_index),'.png']))
%         
%         match_nums = unique(threshed_image.*consensus_label);
%         match_nums = match_nums(match_nums > 0);
%         matches(t2_index,size_index) = length(match_nums);
%         
% %         match_thresh_nums = unique(threshed_image_label.*consensus);
% %         match_thresh_nums = match_thresh_nums(match_thresh_nums > 0);
%         comp_hits(t2_index,size_index) = length(unique(threshed_image_label(:))) - 1 + length(match_nums);
%         
%         1;
%     end
% end
% false_pos = 100*(1 - matches./comp_hits);
% false_neg = 100*(1 - matches./max(consensus_label(:)));
% stats = struct('matches',matches,'comp_hits',comp_hits, ... 
%     'hand_hits',max(consensus_label(:)),'false_pos',false_pos,'false_neg',false_neg); %#ok<NASGU>
% save(fullfile(seg_folder,'stats.mat'),'stats');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expansion Thresholding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t1_vals = 0.5:0.5:5.5;
comp_pixels = zeros(1,length(t1_vals));
consensus_pixels = zeros(1,length(t1_vals));
overlap_pixels = zeros(1,length(t1_vals));
for t1_index = 1:length(t1_vals)
    filter_thresh = puncta_mean + puncta_std*[t1_vals(t1_index), 4];
    threshed_image = find_threshed_image(high_passed,filter_thresh,3);
    threshed_image = imfill(threshed_image,'holes');
    threshed_image_label = bwlabel(threshed_image,8);
    threshed_image_label = filter_label_mat_by_size(threshed_image_label,[-Inf,250]);
    threshed_image = threshed_image_label > 0;
    
    temp = zeros(size(threshed_image)); temp(consensus) = 1; temp(threshed_image) = 2; temp(threshed_image & consensus) = 3;
    
    
    [comp_pixels(t1_index),consensus_pixels(t1_index),overlap_pixels(t1_index)] = ...
        find_pixel_stats(threshed_image,threshed_image_label,consensus,consensus_label);
    1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hand_seg = process_hand_segmentation(file)

hand_seg = double(imread(file));
hand_seg = (hand_seg - min(hand_seg(:)))/range(hand_seg(:));
hand_seg = hand_seg < 0.9;

hand_seg = imfill(hand_seg,'holes');


function threshed_image = find_threshed_image(high_passed_image, filter_thresh, minimum_size)

if (length(filter_thresh) == 1)
    threshed_image = high_passed_image >= filter_thresh;
    threshed_image = remove_edge_objects(threshed_image);
else
    high_threshed_image = high_passed_image >= filter_thresh(2);
    high_threshed_image = remove_edge_objects(high_threshed_image);
    high_threshed_image_label = bwlabel(high_threshed_image,4);
    high_threshed_image_label = filter_label_mat_by_size(high_threshed_image_label,minimum_size);
    high_threshed_image = high_threshed_image_label > 0;
    
    low_threshed_image = high_passed_image >= filter_thresh(1);
    low_thresh_bwlabel = bwlabel(low_threshed_image,4);
    
    overlap_labels = unique(low_thresh_bwlabel.*high_threshed_image);
    overlap_labels = overlap_labels(overlap_labels > 0);
    
    threshed_image = ismember(low_thresh_bwlabel,overlap_labels);
    threshed_image = remove_edge_objects(threshed_image);
end

function [comp_sum,consen_sum,overlap_sum] = find_pixel_stats(comp,comp_label,consensus,consensus_label)

%find the consensus puncta that overlap with the comp puncta,
%throwing out any background overlap
consensus_nums = unique(comp.*consensus_label);
consensus_nums = consensus_nums(consensus_nums > 0);
consensus_overlap = ismember(consensus_label, consensus_nums);
consen_sum = sum(consensus_overlap(:));

%find the comp puncta that overlap with the consensus puncta,
%throwing out any background overlap
comp_nums = unique(comp_label.*consensus);
comp_nums = comp_nums(comp_nums > 0);
comp_overlap = ismember(comp_label, comp_nums);
comp_sum= sum(comp_overlap(:));

dual_overlap = consensus_overlap .* comp_overlap;
overlap_sum = sum(dual_overlap(:));


function filt_image = filter_label_mat_by_size(label_mat,size_mat)

props = regionprops(label_mat,'Area');

if (length(size_mat) == 1)
    filt_image = ismember(label_mat,find([props.Area] >= size_mat(1))).*label_mat;
else
    filt_image = ismember(label_mat,find([props.Area] >= size_mat(1) & [props.Area] <= size_mat(2))).*label_mat;
end