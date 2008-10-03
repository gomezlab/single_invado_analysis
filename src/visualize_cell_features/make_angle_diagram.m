exp_name = 'time_series_08';
tracking_mat = load(fullfile('../../results/focal_adhesions',exp_name,'tracking_matrices/tracking_seq.csv'))+1;
data = load(fullfile('../../results/focal_adhesions',exp_name,'adhesion_props/lin_time_series/Angle_to_center.csv'));
% orient_data = load(fullfile('../../results/focal_adhesions',exp_name,'adhesion_props/lin_time_series/Orientation.csv'));
area_data = load(fullfile('../../results/focal_adhesions',exp_name,'adhesion_props/lin_time_series/Area.csv'));
assert(all(size(tracking_mat) == size(data)),'Error: tracking and data matrices are not the same size')

ad_lifetimes = zeros(size(tracking_mat,1),1);
for i = 1:size(tracking_mat,1), ad_lifetimes(i) = sum(tracking_mat(i,:) >= 1); end

ad_lifetime_remain = zeros(size(tracking_mat));
for i = 1:size(tracking_mat,1)
    tracking_row = tracking_mat(i,:);
    for j = 1:size(tracking_mat,2)
        ad_lifetime_remain(i,j) = sum(tracking_row(j:end) >= 1);
    end
end
%
% rows_to_include = zeros(1,size(tracking_mat,2));
% for i = 1:size(tracking_mat,1)
%     if (sum(tracking_mat(i,:) >= 1) >= 20), rows_to_include(i) = 1; end
% end
% ad_lifetime_remain = ad_lifetime_remain(rows_to_include,:);
% tracking_mat = tracking_mat(rows_to_include,:);
% data = data(rows_to_include,:);

area_sums = zeros(360,size(tracking_mat,2));
for i = 1:size(tracking_mat,1)
    tracking_row = tracking_mat(i,:);
    for j = 1:size(tracking_mat,2)
        %         surrounding_entries = [NaN, tracking_row(j), NaN];
        %         try surrounding_entries(1) = tracking_row(j - 1); end %#ok<TRYNC>
        %         try surrounding_entries(3) = tracking_row(j + 1); end %#ok<TRYNC>
        %         if (surrounding_entries(2) >= 1 && surrounding_entries(1) < 1)
        %             image(floor(data(i,j)*(360/(2*pi)))+1,j,1:3) = image(floor(data(i,j)*(360/(2*pi)))+1,j,1:3) + cat(3,0,0.1,0);
        %         end

        if (not(isnan(data(i,j))))
            angle = round(data(i,j)*(360/(2*pi)));
            if (angle == 0), angle = 360; end
            assert(angle <= 360 & angle > 0,'%d %d',angle, data(i,j))
%             if (all(image(angle,j,1:3) == 1))
%                 image(angle,j,1:3) = cat(3,0,0,0);
%             end
            area_sums(angle,j) = area_sums(angle,j) + area_data(i,j);
        end
    end
end
area_sums = log(area_sums);

image = zeros(360, size(tracking_mat,2), 3);
cmap = hot(1000);
for i = 1:size(image,1)
    for j = 1:size(image,2)
        if (area_sums(i,j) == -Inf), continue; end
        
        scaled_area = (area_sums(i,j) - min(min(area_sums(area_sums ~= -Inf))))/range(area_sums(area_sums ~= -Inf));
        if (scaled_area == 0), scaled_area = 1e-10; end
        assert(scaled_area >= 0 && scaled_area <= 1, '%d %d', scaled_area, area_sums(i,j));
        
        cmap_row = ceil(scaled_area*size(cmap,1));
        assert(cmap_row > 0 && cmap_row <= size(cmap,1),'Cmap row is out of range: %d', cmap_row);
        
        image(i,j,1:3) = cat(3,cmap(cmap_row,1),cmap(cmap_row,2),cmap(cmap_row,3));
    end
end

imwrite(image,['as_',exp_name,'.png'])