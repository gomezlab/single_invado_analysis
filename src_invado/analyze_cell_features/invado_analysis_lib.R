################################################################################
# invaod_analysis_lib.R: the R functions associated with processing the
# invaodopodia data
################################################################################

library(Hmisc);

gather_invado_properties <- function(results_dirs, build_degrade_plots = FALSE, 
    conf.level = 0.95, degrade_file = "Local_gel_diff_corr.csv", results.file = NA, build_plots=TRUE) {
    
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
        pre_diff_data = read.table(file.path(this_exp_dir,'lin_time_series', 'Pre_birth_diff.csv'), 
            sep=",",header=F);
        area_data = read.table(file.path(this_exp_dir,'lin_time_series', 'Area.csv'), 
            sep=",",header=F);
        stopifnot(dim(lineage_data)[[1]] == dim(degrade_data)[[1]])
        
        #Building the filter sets
        
        #exclude puncta that don't live for at least 5 time steps
        longev_filter = ! is.na(lineage_data$longevity) & lineage_data$longevity >= 5;
        
        #exclude puncta that underwent split births
        no_split_birth_filt = ! is.na(lineage_data$split_birth_status) & ! lineage_data$split_birth_status;

        #exclude puncta that disappeared due to merge events (coded as 0 in the
        #lineage data file)
        death_filt = ! is.na(lineage_data$death_status) & lineage_data$death_status;

        overall_filt = longev_filter & no_split_birth_filt & death_filt;
        
        all_props$overall_filt = c(all_props$overall_filt,which(overall_filt));
        all_props$experiment = c(all_props$experiment,rep(basename(dirname(this_exp_dir)),length(which(overall_filt))));
        props_to_include = c("longevity","largest_area","birth_i_num");
        for (i in props_to_include) {
            all_props[[i]] = c(all_props[[i]], lineage_data[[i]][overall_filt]);
        }

        #analyzing each of the puncta in the filtered set for invadopodia filtering processes
        if (build_plots) {
            pdf(file.path(this_exp_dir,'local_degrade_plots.pdf'));
        }
        for (i in which(overall_filt)) {
            only_data = na.omit(as.numeric(degrade_data[i,]));
            
            all_props$mean_vals = c(all_props$mean_vals, mean(only_data))
            all_props$max_local_diff = c(all_props$max_local_diff, max(only_data));
            all_props$min_local_diff = c(all_props$min_local_diff, min(only_data));
            
            test_results = t.test(only_data,conf.level=conf.level);
            
            all_props$low_conf_int = c(all_props$low_conf_int, test_results$conf.int[1]);
            all_props$high_conf_int = c(all_props$high_conf_int, test_results$conf.int[2]);
            all_props$p_value = c(all_props$p_value, test_results$p.value);
            
            only_pre_diff_data = na.omit(as.numeric(pre_diff_data[i,]));
            
            only_pre_diff_test = t.test(only_pre_diff_data,conf.level=conf.level);

            local_pre_test = t.test(only_data - only_pre_diff_data,conf.level=conf.level);
            all_props$local_pre_p_value = c(all_props$local_pre_p_value, local_pre_test$p.value);
            
            only_area_data = na.omit(as.numeric(area_data[i,]));

            if (build_plots) {
                #plotting the results of the tests
                time_points = seq(from=0,by=5,along.with=only_data);
                
                par(bty='n', mar=c(5,4,4,5))
                matplot(time_points, cbind(only_data, only_pre_diff_data), typ='l', xlab='Time (min)', ylab='Local Diff', main=i)
                legend('topleft',c('Local Diff','First Image Local Diff'), fill=c('black','red'))
                segments(0,0,max(time_points),0, col='green', lty=4)
                errbar(max(time_points), test_results$estimate, test_results$conf.int[2], test_results$conf.int[1], add=T)
                errbar(max(time_points), only_pre_diff_test$estimate, only_pre_diff_test$conf.int[2], only_pre_diff_test$conf.int[1], add=T, col='red')
                
                #Adding the areas to the same plot
                # plot_props = par('usr');
                # 
                # scale_factor = (plot_props[4]-plot_props[3])/max(only_area_data)
                # lines(time_points,only_area_data*scale_factor+plot_props[3],col='purple')
                # 
                # tick_labels = (axTicks(2) - plot_props[3])/scale_factor;

                # axis(4,at=axTicks(2),labels=sprintf('%.3f',(axTicks(2) - plot_props[3])/scale_factor))
            }
        }
        if (build_plots) {
            graphics.off();
        }
    
        all_props = as.data.frame(all_props);

        if (! is.na(results.file)) {
            this_file = file.path(this_exp_dir,results.file);
            if (! file.exists(dirname(this_file))) {
                dir.create(dirname(this_file),recursive=TRUE);
            }
            save(all_props,file = this_file);
        }
    }
    return(all_props);
}

build_filter_sets <- function(raw_data_set) {
    filter_sets = list();
    filter_sets$invado_filter = raw_data_set$high_conf_int < 0 & raw_data_set$local_pre_p_value < 0.05;
    filter_sets$non_invado_filter = ! filter_sets$invado_filter;

    return(filter_sets);
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
            results.file = file.path('models','puncta_props_corr.Rdata'));
        
        data_types_to_include = c('overall_filt','p_value','mean_vals');
        
        filter_sets = build_filter_sets(exp_props);

        invado_lineage_data = subset(exp_props, filter_sets$invado_filter, select = data_types_to_include);
        
        non_invado_lineage_data = subset(exp_props, filter_sets$non_invado_filter, select = data_types_to_include);
        
        write.table(invado_lineage_data, file.path(data_dir, 'invado_data.csv'), row.names=F, col.names=F, sep=',')
        write.table(non_invado_lineage_data, file.path(data_dir, 'non_invado_data.csv'), row.names=F, col.names=F, sep=',')
    }
}
