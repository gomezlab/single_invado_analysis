exp_name = 'time_series_23';
tracking_mat = load(fullfile('../../results/focal_adhesions',exp_name,'tracking_matrices/tracking_seq.csv'))+1;
data = load(fullfile('../../results/focal_adhesions',exp_name,'adhesion_props/lin_time_series/Average_adhesion_signal.csv'));
assert(all(size(tracking_mat) == size(data)),'Error: tracking and data matrices are not the same size')

filtered_tracking_mat = zeros(1,size(tracking_mat,2));
filtered_data = zeros(1,size(tracking_mat,2));
for i = 1:size(tracking_mat,1)
    tracking_row = tracking_mat(i,:);
    data_row = data(i,:);
    if (sum(tracking_row >= 1) >= 20)
        filtered_tracking_mat(size(filtered_tracking_mat,1) + 1,:) = tracking_row;
        filtered_data(size(filtered_data,1) + 1,:) = data_row;
    end
end
filtered_tracking_mat = filtered_tracking_mat(2:end,:);
filtered_data = filtered_data(2:end,:);
assert(all(size(filtered_tracking_mat) == size(filtered_data)),'Error: filtered tracking and data matricies are not the same size');

max_data = max(max(filtered_data));
min_data = min(min(filtered_data));

tracking_mat = filtered_tracking_mat;
image = ones([size(tracking_mat),3]);
for i = 1:size(tracking_mat,1)
    tracking_row = tracking_mat(i,:);
    for j = 1:size(tracking_mat,2)
        surrounding_entries = [NaN, tracking_row(j), NaN];
        try surrounding_entries(1) = tracking_row(j - 1); end %#ok<TRYNC>
        try surrounding_entries(3) = tracking_row(j + 1); end %#ok<TRYNC>

        %Birth
        if (surrounding_entries(2) >= 1 && surrounding_entries(1) < 1)
            image(i,j,1:3) = [0,1,0];
        elseif ((surrounding_entries(1) >= 1 && surrounding_entries(3) >= 1) || (surrounding_entries(2) >= 1 && j == 1) || (surrounding_entries(2) >= 1 && j == size(tracking_mat,2)))
            scaled_data_value = (filtered_data(i,j) - min_data)/(max_data - min_data);
            assert(scaled_data_value <= 1 && scaled_data_value >= 0)
            
            image(i,j,1:3) = scaled_data_value*[1,1,1];
            assert(not(isnan(filtered_data(i,j))),'Error: filtered data point (%d, %d) is NaN',i,j);
        elseif (surrounding_entries(2) >= 1 && surrounding_entries(3) < 1)
            image(i,j,1:3) = [1,0,0];
        end
    end
end

imwrite(image,['lineage_sum_',exp_name,'.png'])