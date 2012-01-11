################################################################################
# invado_analysis_lib.R: the R functions associated with processing the
# invaodopodia data
################################################################################

library(Hmisc);

get_non_nan_data <- function(files) {
    data = c();
    
    for (this_file in files) {
        temp = read.csv(this_file,header=F)
        temp = temp[!is.na(temp)];
        data = c(data,temp);

        if (any(temp < -500)) {
            print(this_file)
        }
    }
    
    return(data)
}


################################################################################
# Main Program
################################################################################

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

	print(data_dir);
    if (exists('data_dir')) {
        exp_props = gather_invado_properties(data_dir,
            results.file = file.path('models','puncta_props_corr.Rdata'));
        
        pdf(file.path(data_dir,'p_vals.pdf'))
        hist(exp_props$p_value);
        graphics.off()

        data_types_to_include = c('overall_filt','p_value','mean_local_diff', 'pre_diff_p_value');
        
        filter_sets = build_filter_sets(exp_props);
        
        invado_lineage_data = subset(exp_props, filter_sets$invado_filter, select = data_types_to_include);
        local_diff_invado_lineage_data = subset(exp_props, filter_sets$local_diff_filter, 
            select = data_types_to_include);
        
        not_invado_lineage_data = subset(exp_props, filter_sets$not_invado_filter, 
            select = data_types_to_include);
        
        write.table(invado_lineage_data, file.path(data_dir, 'invado_data.csv'), 
            row.names=F, col.names=F, sep=',')
        write.table(local_diff_invado_lineage_data, file.path(data_dir, 'local_invado_data.csv'), 
            row.names=F, col.names=F, sep=',')
        
        write.table(not_invado_lineage_data, file.path(data_dir, 'not_invado_data.csv'), 
            row.names=F, col.names=F, sep=',')
    }
}
