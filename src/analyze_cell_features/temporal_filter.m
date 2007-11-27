function [varargout] = temporal_filter(folder,num,degree,varargin)
% TEMPORAL_FILTER    filters a given binary feature image based on the
%                    presence of overlapping features in a set of images
%
%   Note: This script assumes that the cell image sequence directory
%   structure is present, i.e. each image in the time course and
%   accompanying data is in it's own zero padded directory
%
%   If an output parameter is requested, the number of features
%   removed from the image is returned.
%
%   temporal_filter(DIR,NUM,DEG,OPT) given a directory (DIR) containing
%   time course images and the specific image number of interest (NUM), a
%   given number of time course images (DEG) are used to filter the image
%   of interest, all of the options (OPT) are used as well
%
%   Options:
%       -'out_dir': specifies a directory where the filtered image will be
%        output, defaults to same folder as the input file
%
%       -'out_filename': specifies a filename to use for the filtered
%        output image, defaults to 'time_filtered.png'
%       
%       -'exclude_file': specifies a filename that holds the image numbers
%        that should be excluded from examination
%
%       -'in_filename': specifies a filename to search for in each
%        directory to perform the temporal filtering, defaults to
%        'adhesions.png'
%
%       -'debug': set to 1 to enable debugging output, includes listing of
%        images used in the temporal filter


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'TEMPORAL_FILTER';

i_p.addRequired('folder',@(x)exist(x,'dir')==7);
i_p.addRequired('num', @(x)isnumeric(x) && x>=1);
i_p.addRequired('degree',@(x)isnumeric(x) && x>=1);

i_p.parse(folder,num,degree);

folder_guess = fullfile(folder,num2str(num));
leading_zeros = 1;
while exist(folder_guess,'dir') ~= 7
    leading_zeros = leading_zeros + 1;
    if (leading_zeros > 50)
        error('TEMPORAL_FILTER: can''t find image number folder');
    end
    folder_guess = fullfile(folder,sprintf(['%0',num2str(leading_zeros),'d'],num));
end

i_p.addOptional('out_dir',folder_guess,@(x)exist(x,'dir')==7);
i_p.addOptional('out_filename','time_filtered.png',@ischar);
i_p.addOptional('in_filename','adhesions.png',@ischar);
i_p.addOptional('exclude_filename',0,@(x)exist(x,'file')==2);
i_p.addOptional('debug',0,@(x)x==1);
i_p.parse(folder,num,degree,varargin{:});

%collect all the parsed input values into shorter variable names and read
%in the feature images

experiment_folder = i_p.Results.folder;
num = i_p.Results.num;
degree = i_p.Results.degree;
out_dir = i_p.Results.out_dir;
out_file = i_p.Results.out_filename;
in_file = i_p.Results.in_filename;
debug = i_p.Results.debug;

excluded_frames = 0;
if (i_p.Results.exclude_filename)
    excluded_frames = load(i_p.Results.exclude_filename);
end

full_out_path = fileparts(fullfile(out_dir,out_file));
if (not(exist(full_out_path,'dir')))
    mkdir(full_out_path)
end

base_image_file = fullfile(experiment_folder,num2str(num),in_file);
leading_zeros = 1;
while exist(base_image_file,'file') ~= 2
    leading_zeros = leading_zeros + 1;
    if (leading_zeros > 50)
        error('TEMPORAL_FILTER: can''t find base image file');
    end
    base_image_file = fullfile(experiment_folder,sprintf(['%0',num2str(leading_zeros),'d'],num),in_file);
end

padded_base_num = sprintf(['%0',num2str(leading_zeros),'d'],num);
base_image = imread(base_image_file);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temporal_image_set = ones(size(base_image,1),size(base_image,2),degree);

if (debug), disp(['I_num: ',num2str(num)]), end;

offset = 1;
i_gathered = 0;
while i_gathered < degree
    image_num = offset + num;
    if (find(image_num==excluded_frames))
        if debug, disp(['e_frame: ',num2str(image_num),'  ']), end;        
        offset = offset + 1;
        continue;
    end
    if (exist(fullfile(experiment_folder,sprintf(['%0',num2str(leading_zeros),'d'],image_num),in_file),'file'))
        if (debug), disp(image_num), end;

        i_gathered = i_gathered + 1;
        offset = offset + 1;
        temporal_image_set(:,:,i_gathered) = imread(fullfile(experiment_folder,sprintf('%03d',image_num),in_file));
    else  
        i_gathered = i_gathered + 1;
    end
end

intersected_image = base_image & temporal_image_set(:,:,1);
for i=2:degree
    intersected_image = intersected_image & temporal_image_set(:,:,i);
end

labeled_base = bwlabel(base_image);
out_image = zeros(size(base_image));

features_removed = 0;
for i=1:max(labeled_base(:))
    if (sum(intersected_image(find(labeled_base==i))))
        out_image(find(labeled_base==i)) = 1;
    else
        features_removed = features_removed + 1;
    end
end

imwrite(out_image,fullfile(out_dir,out_file));

if (nargout >= 1)
    varargout(1) = {features_removed};
end

end