################################################################################
# invaod_analysis_lib.R: the R functions associated with processing the
# invaodopodia data
################################################################################

gather_invado_properties <- function(results_dirs, build_degrade_plots = FALSE, 
    conf.level = 0.95, degrade_file = "Local_gel_diff.csv", results.file = NA) {
    
    all_props = list();

    for (this_exp_dir in results_dirs) {
        if (! file.exists(file.path(this_exp_dir,'single_lin.csv'))) {
            warning(paste("Unable to find", file.path(this_exp_dir,'single_lin.csv'),
                  "moving on to the next file."))
            next;
        }
        
        #Reading in raw data
        print(paste('Working on:',this_exp_dir));
        lineage_data = read.table(file.path(this_exp_dir,'single_lin.csv'), 
            sep=",",header=T);
        degrade_data = read.table(file.path(this_exp_dir,'lin_time_series', degrade_file), 
            sep=",",header=F);
        stopifnot(dim(lineage_data)[[1]] == dim(degrade_data)[[1]])
        
        #Building the filter sets
        longev_filter = ! is.na(lineage_data$longevity) & lineage_data$longevity >= 5;
        no_split_birth_filt = ! is.na(lineage_data$split_birth_status) & ! lineage_data$split_birth_status;
        death_filt = ! is.na(lineage_data$death_status) & lineage_data$death_status;

        overall_filt = longev_filter & no_split_birth_filt & death_filt;
        
        all_props$overall_filt = c(all_props$overall_filt,which(overall_filt));
        all_props$experiment = c(all_props$experiment,rep(basename(dirname(this_exp_dir)),length(which(overall_filt))));
        props_to_include = c("longevity","largest_area","birth_i_num");
        for (i in props_to_include) {
            all_props[[i]] = c(all_props[[i]], lineage_data[[i]][overall_filt]);
        }

        time_to_conf_neg = c();
        time_to_neg = c();
        for (i in which(overall_filt)) {
            only_data = na.omit(as.numeric(degrade_data[i,]));
            
            all_props$mean_vals = c(all_props$mean_vals, mean(only_data))
            all_props$max_local_diff = c(all_props$max_local_diff, max(only_data));
            all_props$min_local_diff = c(all_props$min_local_diff, min(only_data));
            
            test_results = t.test(only_data,conf.level=conf.level);
            
            all_props$low_conf_int = c(all_props$low_conf_int, test_results$conf.int[1]);
            all_props$high_conf_int = c(all_props$high_conf_int, test_results$conf.int[2]);
            all_props$p_values = c(all_props$p_values, test_results$p.value);

            #search for the minimum number of points (at least five) needed to
            #verify that the local gel difference is negative
            if (test_results$conf.int[2] < 0) {
                time_to_neg = c(time_to_neg, which(only_data < 0)[[1]])

                min_found = 0;
                for (i in 5:length(only_data)) {
                    short_conf_int = t.test(only_data[1:i],conf.level=conf.level)$conf.int;
                    if (short_conf_int[2] < 0 && ! min_found) {
                        min_found = 1;
                        time_to_conf_neg = c(time_to_conf_neg, i);
                    }
                }
            } else {
                time_to_conf_neg = c(time_to_conf_neg, NA);
                time_to_neg = c(time_to_neg, NA)
            }
        }
        all_props$time_to_conf_neg = c(all_props$time_to_conf_neg, time_to_conf_neg);
        all_props$time_to_neg = c(all_props$time_to_neg, time_to_neg);
        
        if (! is.na(results.file)) {
            this_file = file.path(this_exp_dir,results.file);
            if (! file.exists(dirname(this_file))) {
                dir.create(dirname(this_file),recursive=TRUE);
            }
            all_props_frame = as.data.frame(all_props);
            save(all_props_frame,file = this_file);
        }
    }
    return(as.data.frame(all_props));
}

################################################################################
# Main Program
################################################################################

args = commandArgs(TRUE);
if (length(args) != 0) {
    debug = FALSE;
    min_length = 10;
    
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
            results.file = file.path('models','puncta_props.Rdata'));
        exp_props = gather_invado_properties(data_dir, 
            results.file = file.path('models','puncta_props_corr.Rdata'), degrade_file="Local_gel_diff_corr.csv");
        
        invado_lineage_nums = exp_props$overall_filt[exp_props$high_conf_int < 0];
        invado_p_values = exp_props$p_values[exp_props$high_conf_int < 0];
        invado_mean_vals = exp_props$mean_vals[exp_props$high_conf_int < 0];
        non_invado_lineage_nums = exp_props$overall_filt[exp_props$high_conf_int >= 0];
        non_invado_p_values = exp_props$p_values[exp_props$high_conf_int >= 0];
        non_invado_mean_vals = exp_props$mean_vals[exp_props$high_conf_int >= 0];
        
        write.table(cbind(invado_lineage_nums,invado_p_values,invado_mean_vals), file.path(data_dir, 'invado_nums.csv'), row.names=F, col.names=F, sep=',')
        write.table(cbind(non_invado_lineage_nums,non_invado_p_values, non_invado_mean_vals), file.path(data_dir, 'non_invado_nums.csv'), row.names=F, col.names=F, sep=',')
    }
}
