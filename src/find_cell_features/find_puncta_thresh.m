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
single_image_thresh = [];

for i_num = 1:size(image_dirs)
    image = double(imread(fullfile(base_dir,image_dirs(i_num).name,filenames.puncta)));
    
    I_filt = fspecial('disk',11);
    blurred_image = imfilter(image,I_filt,'same',mean(image(:)));
    high_passed_image = image - blurred_image;
    
    single_image_thresh = [single_image_thresh ...
        mean(high_passed_image(:)) + i_p.Results.stdev_thresh*std(high_passed_image(:))]; %#ok<AGROW>
    all_filt = [all_filt high_passed_image(:)']; %#ok<AGROW>
    
    if (mod(i_num,10)==0)
        disp(['Done with ',num2str(i_num),'/',num2str(size(image_dirs,1))])
    end
end
all_filt = double(all_filt);
plot(single_image_thresh)
xlabel('Image Number');
ylabel('Puncta Threshold');
print('-depsc2', fullfile(base_dir,'..','puncta_props','single_image_thresholds.eps'));
close;

mean_all_filt = mean(all_filt(:));
stdev_all_filt = std(all_filt(:));

csvwrite(fullfile(base_dir,image_dirs(1).name,filenames.puncta_threshold), ...
    [mean_all_filt,stdev_all_filt]);

hist(all_filt(:),100);
xlabel('High Pass Filtered Intensity','FontSize',16,'FontName','Helvetica');
ylabel('Pixel Count','FontSize',16,'FontName','Helvetica');
y_limits = ylim();

for i=1:5
	line([i*stdev_all_filt,i*stdev_all_filt],[0,y_limits(2)],'Color','red', ...
		'LineStyle','--','LineWidth',3);
end

set(gca, 'FontName','Helvetica','FontSize',16,'Box','off');
set(gcf, 'PaperPositionMode', 'auto');
print('-depsc2', fullfile(base_dir,image_dirs(1).name,filenames.puncta_threshold_plot));
close;

toc;

if (i_p.Results.thresh_scan)
    props_folder = fileparts(fullfile(base_dir,image_dirs(1).name,filenames.puncta_threshold));
    props_folder = fullfile(props_folder,'thresh_scan_seed');
    if (not(exist(props_folder,'dir')))
        mkdir(props_folder);
    end
    
    mean_val = mean(all_filt(:));
    stdev_val = std(all_filt(:));
    
    for std_count = 0.5:0.1:20
        thresh_val = mean_val + std_count*stdev_val;
        output_filename = sprintf('%0.1f.csv',std_count);
        csvwrite(fullfile(props_folder,output_filename),thresh_val);
    end
    
    props_folder = fileparts(fullfile(base_dir,image_dirs(1).name,filenames.puncta_threshold));
    props_folder = fullfile(props_folder,'thresh_scan_expand');
    if (not(exist(props_folder,'dir')))
        mkdir(props_folder);
    end
    
    for std_count = 0.1:0.1:4.9
        thresh_val = mean_val + stdev_val*[std_count,5];
        output_filename = sprintf('%0.1f.csv',std_count);
        csvwrite(fullfile(props_folder,output_filename),thresh_val);
    end
end
