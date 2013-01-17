function build_dual_highlight_image(exp_dir,varargin)
%BUILD_DUAL_HIGHLIGHT_IMAGES    Builds highlighted images showing the
%                               invadopodia (green), non-invadopodia (red),
%                               other puncta (blue) and the cell edge
%                               (purple), placing the image for each frame
%                               of the movie in the individual image
%                               folders

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'BUILD_DUAL_HIGHLIGHT_IMAGES';

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);
% i_p.addRequired('invado_file',@(x)exist(x,'file') == 2);
% i_p.addRequired('non_invado_file',@(x)exist(x,'file') == 2);
i_p.addParamValue('no_scale_bar',0,@(x) islogical(x) || x == 0 || x == 1);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

profile off; profile on;

addpath(genpath('..'))
filenames = add_filenames_to_struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find all the image directories
I_folder = fullfile(exp_dir,'individual_pictures');
image_dirs = dir(I_folder);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tracking matrix reading/filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
invado_data_file = fullfile(I_folder,image_dirs(1).name,filenames.invado_data);
if (exist(invado_data_file,'file'))
    invado_data = csvread(invado_data_file,1,0);
else
    invado_data = zeros(0,3);
end

not_invado_data_file = fullfile(I_folder,image_dirs(1).name,filenames.not_invado_data);
if (exist(not_invado_data_file,'file'))
    not_invado_data = csvread(not_invado_data_file,1,0);
else
    not_invado_data = zeros(0,3);
end

%after loading the tracking sequence filter to only include those puncta
%included in the invadopodia list, remember the list formated so the first
%column contains the lineage number, with the first lineage as 1, so no
%need to translate
tracking_seq = load(fullfile(I_folder,image_dirs(1).name,filenames.tracking_matrix)) + 1;
invado_tracking_seq = tracking_seq(invado_data(:,1),:);
not_invado_tracking_seq = tracking_seq(not_invado_data(:,1),:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image Creation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gel_limits = csvread(fullfile(I_folder,image_dirs(end).name,filenames.gel_range));
puncta_limits = csvread(fullfile(I_folder,image_dirs(end).name,filenames.puncta_range));

output_folder = fullfile(exp_dir,'visualizations','invado_and_not');
if (not(exist(output_folder,'dir')))
    mkdir(output_folder);
end

output_folder_small = fullfile(exp_dir,'visualizations','invado_and_not_small');
if (not(exist(output_folder_small,'dir')))
    mkdir(output_folder_small);
end

output_folder_side = fullfile(exp_dir,'visualizations','side_by_side');
if (not(exist(output_folder_side,'dir')))
    mkdir(output_folder_side);
end

for i = 1:length(image_dirs)
    %check for the presence of adhesions to map onto the last gel image, if
    %not present, skip to the next image folder
    if (not(any(tracking_seq(:,i) > 0)))
        continue;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Gather the adhesion label image and perimeters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        gel_image = double(imread(fullfile(I_folder,image_dirs(i).name,filenames.gel)));
    gel_image = gel_image - gel_limits(1);
    gel_image = gel_image .* (1/gel_limits(2));
    gel_image = cat(3,gel_image,gel_image,gel_image);
    
    puncta_image = double(imread(fullfile(I_folder,image_dirs(i).name,filenames.puncta)));
    puncta_image = puncta_image - puncta_limits(1);
    puncta_image = puncta_image .* (1/puncta_limits(2));
    puncta_image = cat(3,puncta_image,puncta_image,puncta_image);
    
    cell_mask = bwperim(logical(imread(fullfile(I_folder,image_dirs(i).name,filenames.cell_mask))));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Gather the adhesion label perimeters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    puncta_label_perim = imread(fullfile(I_folder,image_dirs(i).name,filenames.objects));
    
    puncta_nums = tracking_seq(:,i);
    puncta_nums = puncta_nums(puncta_nums > 0);
    assert(max(puncta_nums) == length(puncta_nums));
    
    invado_nums = invado_tracking_seq(:,i);
    invado_nums = invado_nums(invado_nums > 0);
    
    not_invado_nums = not_invado_tracking_seq(:,i);
    not_invado_nums = not_invado_nums(not_invado_nums > 0);
    
    puncta_label_perim_invado = ismember(puncta_label_perim,invado_nums);

    puncta_label_perim_not_invado = zeros(size(puncta_label_perim));
    puncta_label_perim_not_invado(ismember(puncta_label_perim,not_invado_nums)) = puncta_label_perim(ismember(puncta_label_perim,not_invado_nums));
    puncta_label_perim_not_invado = im2bw(puncta_label_perim_not_invado,0);
    
    puncta_label_perim_neither = puncta_label_perim;
    puncta_label_perim_neither(puncta_label_perim_invado > 0) = 0;
    puncta_label_perim_neither(puncta_label_perim_not_invado > 0) = 0;
    puncta_label_perim_neither = im2bw(puncta_label_perim_neither,0);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Image Creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    images_to_highlight = {gel_image, puncta_image};
    original_images = images_to_highlight;
    spacer = ones(size(images_to_highlight{1},1),1,3);
    
    for j=1:length(images_to_highlight)
        if (any(images_to_highlight{j}(:) > 1))
            1;
        end
        %cell edge highlighting
        images_to_highlight{j} = create_highlighted_image(images_to_highlight{j},cell_mask,'color_map',[150/255,46/255,166/255]);
        
        %not-invadopodia highlighting
        images_to_highlight{j} = create_highlighted_image(images_to_highlight{j},puncta_label_perim_neither,'color_map',[0,0,1]);
        
        %not-invadopodia highlighting
        images_to_highlight{j} = create_highlighted_image(images_to_highlight{j},puncta_label_perim_not_invado,'color_map',[1,0,0]);
        
        %invadopodia highlighting
        images_to_highlight{j} = create_highlighted_image(images_to_highlight{j},puncta_label_perim_invado,'color_map',[0,1,0]);
    end
    
    output_image = [images_to_highlight{1}, spacer, images_to_highlight{2}];
    imwrite(output_image, fullfile(output_folder,[sprintf('%04d',i),'.png']));
    
    output_image_small = imresize(output_image,0.5);
    imwrite(output_image_small, fullfile(output_folder_small,[sprintf('%04d',i),'.png']));
    
    output_image = [original_images{1}, spacer, original_images{2}];
    output_image = imresize(output_image,0.5);
    imwrite(output_image, fullfile(output_folder_side,[sprintf('%04d',i),'.png']));
    
    if(mod(i,10) == 0), disp(['Done with Image: ',num2str(i),'/',num2str(length(image_dirs))]); end
end

profile off;
if (i_p.Results.debug), profile viewer; end