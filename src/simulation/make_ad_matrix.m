function adhesion_matrix = make_ad_matrix(radius, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'make_ad_matrix';

i_p.addRequired('radius',@(x)x >= 0);
i_p.addParamValue('mean_intensity',0.4237,@isnumeric);
i_p.addParamValue('outer_intensity',0.3067,@isnumeric);

i_p.parse(radius, varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (i_p.Results.radius == 0)    
    adhesion_matrix = [i_p.Results.mean_intensity];
    return;
end

adhesion_matrix = zeros(i_p.Results.radius*2+1);
center = ceil(size(adhesion_matrix,1)/2);
adhesion_matrix(center,center) = 1;
% adhesion_matrix(bwdist(adhesion_matrix) > i_p.Results.radius) = NaN;
adhesion_matrix(bwdist(adhesion_matrix) <= i_p.Results.radius) = 1;

ad_matrix_temp = adhesion_matrix;
pixel_counts = [];
while(any(ad_matrix_temp(:)))
    temp_erode = imerode(ad_matrix_temp,strel('disk', 1));
    pixel_counts = [pixel_counts, sum(ad_matrix_temp(:)) - sum(temp_erode(:))]; %#ok<AGROW>
    ad_matrix_temp = temp_erode;
end

assert(sum(pixel_counts) == sum(adhesion_matrix(:)))

intensity_eq = '';

for i=1:length(pixel_counts)
    frac = ['(', num2str(pixel_counts(i)), '/', num2str(sum(adhesion_matrix(:))), ')'];
    intensity = ['(', num2str(i_p.Results.outer_intensity), ' + (', num2str(i), ' - 1) * x)'];
    intensity_eq = [intensity_eq, frac, '*', intensity, ' + ']; %#ok<AGROW>
end

intensity_eq = [intensity_eq, '(-1*', num2str(i_p.Results.mean_intensity), ')'];
step_size = solve(intensity_eq);

ad_matrix_temp = adhesion_matrix;
for i=1:length(pixel_counts)
    temp_erode = imerode(ad_matrix_temp,strel('disk', 1));
    
    this_layer_intensity = i_p.Results.outer_intensity + (i-1)*step_size;
    
    assert(sum(sum(ad_matrix_temp & not(temp_erode))) == pixel_counts(i))
    
    adhesion_matrix(ad_matrix_temp & not(temp_erode)) = this_layer_intensity;
    ad_matrix_temp = temp_erode;
end