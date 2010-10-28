coord_data = read.csv('stage_corners.STG',header=F,skip=4);

#coord_data = coord_data[5:dim(coord_data)[[1]],];

stage_positions = unique(as.numeric(coord_data$V8));

output_data = c()

for (i in stage_positions) {
	this_stage_data = subset(coord_data, V8 == i);

	this_stage_corners = subset(this_stage_data, V4 != 1111);
	this_stage_not_corners = subset(this_stage_data, V4 == 1111);

	z_focus_model = lm(V4 ~ V2 + V3, data=this_stage_corners)
	
	z_focus_predict = predict(z_focus_model, this_stage_not_corners);
	
	z_focus_predict = as.numeric(sprintf("%0.2f",z_focus_predict))
	
	this_stage_data$V4[this_stage_data$V4 == 1111] = z_focus_predict
	this_stage_data$V6[this_stage_data$V6 == 1111] = z_focus_predict	
	output_data = rbind(output_data, this_stage_data);
}

write.table(output_data,'stage_pos_predicted.csv',row.names=F,col.names=F,sep=',')