function build_dual_highlight_image(cfg_file,invado_file,non_invado_file,varargin)
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

i_p.addRequired('cfg_file',@(x)exist(x,'file') == 2);
i_p.addRequired('invado_file',@(x)exist(x,'file') == 2);
i_p.addRequired('non_invado_file',@(x)exist(x,'file') == 2);
i_p.addParamValue('no_scale_bar',0,@(x) islogical(x) || x == 0 || x == 1);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(cfg_file,invado_file,non_invado_file,varargin{:});

if (i_p.Results.debug == 1), profile off; profile on; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process config file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(cfg_file);
while 1
    line = fgetl(fid);
    if ~ischar(line), break; end
    eval(line);
end

addpath(genpath(path_folders));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tracking matrix reading/filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
invado_data = csvread(invado_file);
non_invado_data = csvread(non_invado_file);

%after loading the tracking sequence filter to only include those puncta
%included in the invadopodia list, remember the list formated so the first
%column contains the lineage number, with the first lineage as 1, so no
%need to translate
tracking_seq = load(tracking_seq_file) + 1;
invado_tracking_seq = tracking_seq(invado_data(:,1),:);
non_invado_tracking_seq = tracking_seq(non_invado_data(:,1),:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%final all the image directories
image_dirs = dir(I_folder);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);

gel_limits = csvread(fullfile(I_folder,image_dirs(end).name,'gel_image_range.csv'));
puncta_limits = csvread(fullfile(I_folder,image_dirs(end).name,'gel_image_range.csv'));

for i = 1:length(image_dirs)
    %check for the presence of adhesions to map onto the last gel image, if
    %not present, skip to the next image folder
    if (not(any(tracking_seq(:,i) > 0)))
        continue;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Gather the adhesion label image and perimeters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    binary_shift = logical(imread(fullfile(I_folder,image_dirs(end).name,'binary_shift.png')));
    
    %read in the gel image and normalize to 0-1
    gel_image = double(imread(fullfile(I_folder,image_dirs(i).name,gel_image_filename)));
    gel_image = gel_image - gel_limits(1);
    gel_image = gel_image .* (1/gel_limits(2));
    gel_image(not(binary_shift)) = 0;
    gel_image = cat(3,gel_image,gel_image,gel_image);
    
    puncta_image = double(imread(fullfile(I_folder,image_dirs(i).name,puncta_image_filename)));
    puncta_image = puncta_image - gel_limits(1);
    puncta_image = puncta_image .* (1/puncta_limits(2));
    puncta_image(not(binary_shift)) = 0;
    puncta_image = cat(3,puncta_image,puncta_image,puncta_image);
    
    cell_mask = bwperim(logical(imread(fullfile(I_folder,image_dirs(i).name,edge_filename))));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Gather the adhesion label perimeters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ad_label_perim = imread(fullfile(I_folder,image_dirs(i).name,adhesions_perim_filename));
    
    ad_nums = tracking_seq(:,i);
    ad_nums = ad_nums(ad_nums > 0);
    assert(max(ad_nums) == length(ad_nums));
    
    invado_nums = invado_tracking_seq(:,i);
    invado_nums = invado_nums(invado_nums > 0);
    
    non_invado_nums = non_invado_tracking_seq(:,i);
    non_invado_nums = non_invado_nums(non_invado_nums > 0);
    
    ad_label_perim_invado = zeros(size(ad_label_perim));
    ad_label_perim_invado(ismember(ad_label_perim,invado_nums)) = ad_label_perim(ismember(ad_label_perim,invado_nums));
    ad_label_perim_invado = im2bw(ad_label_perim_invado,0);
    
    ad_label_perim_non_invado = zeros(size(ad_label_perim));
    ad_label_perim_non_invado(ismember(ad_label_perim,non_invado_nums)) = ad_label_perim(ismember(ad_label_perim,non_invado_nums));
    ad_label_perim_non_invado = im2bw(ad_label_perim_non_invado,0);
    
    ad_label_perim_neither = ad_label_perim;
    ad_label_perim_neither(ad_label_perim_invado > 0) = 0;
    ad_label_perim_neither(ad_label_perim_non_invado > 0) = 0;
    ad_label_perim_neither = im2bw(ad_label_perim_neither,0);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Image Creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    images_to_highlight = {gel_image, puncta_image};
    
    for j=1:length(images_to_highlight)
        %cell edge highlighting
        images_to_highlight{j} = create_highlighted_image(images_to_highlight{j},cell_mask,'color_map',[150/255,46/255,166/255]);
        
        %non-invadopodia highlighting
        images_to_highlight{j} = create_highlighted_image(images_to_highlight{j},ad_label_perim_neither,'color_map',[0,0,1]);
        
        %non-invadopodia highlighting
        images_to_highlight{j} = create_highlighted_image(images_to_highlight{j},ad_label_perim_non_invado,'color_map',[1,0,0]);
        
        %invadopodia highlighting
        images_to_highlight{j} = create_highlighted_image(images_to_highlight{j},ad_label_perim_invado,'color_map',[0,1,0]);
    end
    
    spacer = ones(size(images_to_highlight{1},1),1,3);
    
    output_image = [images_to_highlight{1}, spacer, images_to_highlight{2}];
    imwrite(output_image, fullfile(I_folder,image_dirs(i).name,'invado_and_not.png'));
    
    if(i_p.Results.debug), disp(i); end
end

profile off;
if (i_p.Results.debug), profile viewer; end
