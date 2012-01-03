function find_nuclei(exp_dir, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

addpath('matlab_scripts');
addpath('../visualize_cell_features/');
filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fields = dir(exp_dir);
fields = filter_to_time_series(fields);
fields = fields(1);

for j=1:length(fields)
    base_dir = fullfile(exp_dir,fields(j).name,'individual_pictures');
    
    image_dirs = dir(base_dir);
    
    assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
    assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
    assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>
    
    image_dirs = image_dirs(3:end);
    
    puncta_range = csvread(fullfile(base_dir,image_dirs(1).name,filenames.puncta_range));
    
    for i = 1:size(image_dirs)
        puncta_image = double(imread(fullfile(base_dir,image_dirs(i).name,filenames.puncta)));
        nuc_image = double(imread(fullfile(base_dir,image_dirs(i).name,filenames.nucleus)));        
        
        puncta_image_cor = puncta_image .* (mean(nuc_image(:))/mean(puncta_image(:)));
        
        nuc_image_bleed = nuc_image - puncta_image_cor;
        
        nuclei = nuc_image_bleed > 15;
        nuclei_label = bwlabel(nuclei,4);
        
        nuclei_props = regionprops(nuclei_label,'Area'); %#ok<MRPBW>
        
        nuclei = ismember(nuclei_label,find([nuclei_props.Area] > 25));
        nuclei = imfill(nuclei,'holes');
        imwrite(nuclei,fullfile(base_dir,image_dirs(i).name,filenames.nuclei_binary));
        
        puncta_image = (puncta_image - puncta_range(1))./range(puncta_range);
        nuc_highlight = create_highlighted_image(puncta_image,nuclei);
        imwrite(nuc_highlight,fullfile(base_dir,image_dirs(i).name,'nuclei_highlight.png'));
    end
    disp(['Done with ', fields(j).name])
end
toc;