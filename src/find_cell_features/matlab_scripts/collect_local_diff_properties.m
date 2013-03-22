function diffs = collect_local_diff_properties(data_struct,this_ad,gel_min_val,varargin)
% COLLECT_LOCAL_DIFF_PROPERTIES    find the local difference properties of
%                                  a provided set of images and a given
%                                  puncta

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;
i_p.FunctionName = 'COLLECT_LOCAL_DIFF_PROPERTIES';

i_p.addRequired('data_struct',@isstruct);
i_p.addRequired('this_ad',@islogical);
i_p.addRequired('gel_min_val',@(x)isnumeric(x));

i_p.addOptional('debug',0,@(x)x == 1 || x == 0);

i_p.parse(data_struct,this_ad,gel_min_val,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diffs = struct();

background_region = logical(imdilate(this_ad,strel('disk',5,0)));
%we don't want to include any areas that have been identified as other
%objects
background_region = and(background_region,not(data_struct.objects));

diffs.gel_intensity_puncta = mean(data_struct.gel_image(this_ad)) - gel_min_val;

diffs.Global_gel_diff = mean(data_struct.gel_image(this_ad)) - mean(data_struct.gel_image(:));

if (sum(sum(background_region)) == 0)
    diffs.Local_gel_diff = NaN;
    diffs.Local_gel_diff_percent = NaN;
else
    diffs.background_mean = mean(data_struct.gel_image(background_region)) - gel_min_val;
    diffs.Local_gel_diff = diffs.background_mean - diffs.gel_intensity_puncta;
    diffs.Local_gel_diff_percent = 100*(diffs.Local_gel_diff/diffs.background_mean);
end