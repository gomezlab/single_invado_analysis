function find_puncta_thresh(exp_dir, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('stdev_thresh',4,@(x)isnumeric(x));
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);
i_p.addParamValue('thresh_scan',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

addpath('matlab_scripts');
filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_dir = fullfile(exp_dir,'individual_pictures');
image_dirs = dir(base_dir);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);

all_filt = [];

for i = 1:size(image_dirs)
    image = double(imread(fullfile(base_dir,image_dirs(i).name,filenames.gel)));

    binary_shift = logical(imread(fullfile(base_dir,image_dirs(i).name,filenames.binary_shift)));
    
    only_reg_image = remove_region_outside_registered(image,binary_shift);

    I_filt = fspecial('disk',11);
    blurred_image = imfilter(only_reg_image,I_filt,'same',mean(only_reg_image(:)));
    high_passed_image = only_reg_image - blurred_image;
    
    all_filt = [all_filt high_passed_image(:)']; %#ok<AGROW>
    disp(['Done with ',num2str(i),'/',num2str(size(image_dirs,1))]);
end
all_filt = double(all_filt);

threshold = mean(all_filt(:)) + i_p.Results.stdev_thresh*std(all_filt(:));

csvwrite(fullfile(base_dir,image_dirs(1).name,filenames.puncta_threshold),threshold)

hist(all_filt(:),100);
xlabel('High Pass Filtered Intensity','FontSize',16,'FontName','Helvetica');
ylabel('Pixel Count','FontSize',16,'FontName','Helvetica');
y_limits = ylim();

for i=1:length(threshold)
    this_thresh = threshold(i);
    line([this_thresh,this_thresh],[0,y_limits(2)],'Color','red', ... 
        'LineStyle','--','LineWidth',3);
end

set(gca, 'FontName','Helvetica','FontSize',16,'Box','off');
set(gcf, 'PaperPositionMode', 'auto');
print('-depsc2', fullfile(base_dir,image_dirs(1).name,filenames.puncta_threshold_plot));
close;

toc;

if (i_p.Results.thresh_scan)
    props_folder = fileparts(fullfile(base_dir,image_dirs(1).name,filenames.puncta_threshold));
    props_folder = fullfile(props_folder,'thresh_scan');
    if (not(exist(props_folder,'dir')))
        mkdir(props_folder);
    end
    
    mean_val = mean(all_filt(:));
    stdev_val = std(all_filt(:));
    
    for low_thresh = 0.5:0.1:5
        for high_thresh = 0.5:0.1:5
            thresh_vals = mean_val + [low_thresh high_thresh]*stdev_val;
            output_filename = [num2str(low_thresh),'_',num2str(high_thresh),'.csv'];
            csvwrite(fullfile(props_folder,output_filename),thresh_vals);
        end
    end
end