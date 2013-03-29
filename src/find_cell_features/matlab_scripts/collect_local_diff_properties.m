function diffs = collect_local_diff_properties(data_struct,this_object,varargin)
% COLLECT_LOCAL_DIFF_PROPERTIES    find the local difference properties of
%                                  a provided set of images and a given
%                                  puncta

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;
i_p.FunctionName = 'COLLECT_LOCAL_DIFF_PROPERTIES';

i_p.addRequired('data_struct',@isstruct);
i_p.addRequired('this_object',@islogical);

i_p.addOptional('debug',0,@(x)x == 1 || x == 0);

i_p.parse(data_struct,this_object,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diffs = struct();

background_region = logical(imdilate(this_object,strel('disk',5,0)));
%we don't want to include any areas that have been identified as other
%objects, need to include this_object to deal with pre-birth diff calcs
background_region = background_region & not(data_struct.objects) & not(this_object);

diffs.gel_intensity_puncta = mean(data_struct.gel_image(this_object));

if (sum(sum(background_region)) == 0)
    diffs.gel_background_intensity = NaN;
    diffs.Local_gel_diff = NaN;
else
    diffs.gel_background_intensity = mean(data_struct.gel_image(background_region));
    diffs.Local_gel_diff = diffs.gel_background_intensity - diffs.gel_intensity_puncta;
end