#default corner filename, used mostly for debugging from the R command line
corner_file = 'determined_corners.STG'

## Command line argument processing
args = commandArgs(TRUE);
if (length(args) != 0) {
    debug = FALSE;
    
	#split out the arguments from the passed in parameters and assign variables 
	#in the current scope
    for (this_arg in commandArgs()) {
        split_arg = strsplit(this_arg,"=",fixed=TRUE)
        if (length(split_arg[[1]]) == 1) {
            assign(split_arg[[1]][1], TRUE);
        } else {
            assign(split_arg[[1]][1], split_arg[[1]][2]);
        }
    }
}

#the first four lines are header data that we aren't interested in
coord_data = read.csv(corner_file,header=F,skip=4);

stage_positions = unique(as.numeric(coord_data$V8));

output_data = c()

for (i in stage_positions) {
	this_stage_data = subset(coord_data, V8 == i);
    
    #1111 is the special numeric code used to identify non-corner data
	this_stage_corners = subset(this_stage_data, V4 != 1111);
    stopifnot(dim(this_stage_corners)[[1]] == 4)

	this_stage_not_corners = subset(this_stage_data, V4 == 1111);

    # V2 - X 
    # V3 - Y 
    # V4 - Z
	z_focus_model = lm(V4 ~ V2 + V3, data=this_stage_corners)
	
	z_focus_predict = predict(z_focus_model, this_stage_not_corners);
	
	z_focus_predict = as.numeric(sprintf("%0.2f",z_focus_predict))
	
	this_stage_data$V4[this_stage_data$V4 == 1111] = z_focus_predict
	this_stage_data$V6[this_stage_data$V6 == 1111] = z_focus_predict	
	output_data = rbind(output_data, this_stage_data);
}

write.table(output_data,'stage_pos_predicted.csv',row.names=F,col.names=F,sep=',')
