################################################################################
# invaod_analysis_lib.R: the R functions associated with processing the
# invaodopodia data
################################################################################

gather_invado_properties <- function(results_dirs, build_degrade_plots = FALSE, conf.level = 0.95) {
    all_props = list();

    for (this_exp_dir in results_dirs) {
        if (! file.exists(file.path(this_exp_dir,'single_lin.csv'))) {
            warning(paste("Unable to find", file.path(this_exp_dir,'single_lin.csv'),
                  "moving on to the next file."))
            next;
        }
        
        print(paste('Working on:',this_exp_dir));
        lineage_data = read.table(file.path(this_exp_dir,'single_lin.csv'), 
            sep=",",header=T);
        degrade_data = read.table(file.path(this_exp_dir,'lin_time_series/Local_gel_diff.csv'), 
            sep=",",header=F);
        stopifnot(dim(lineage_data)[[1]] == dim(degrade_data)[[1]])

        longev_filter = ! is.na(lineage_data$longevity) & lineage_data$longevity >= 5;
        all_props$longev_filter = c(all_props$longev_filter,which(longev_filter));
        all_props$longevity = c(all_props$longevity, lineage_data$longevity[longev_filter]);

        final_conf_ints = matrix(ncol=2);
        for (i in which(longev_filter)) {
            only_data = na.omit(as.numeric(degrade_data[i,]));
            conf_int = t.test(only_data,conf.level=conf.level)$conf.int;
            final_conf_ints = rbind(final_conf_ints, conf_int[1:2]); 
        }
        final_conf_ints = final_conf_ints[2:dim(final_conf_ints)[[1]],];
        
        all_props$low_conf_int = c(all_props$low_conf_int, final_conf_ints[,1]);
        all_props$high_conf_int = c(all_props$high_conf_int, final_conf_ints[,2]);

        if (build_degrade_plots) {
            min_diff = min(degrade_data[which(longev_filter),], na.rm=T);
            max_diff = max(degrade_data[which(longev_filter),], na.rm=T);
            
            filtered_longev = lineage_data$longevity[longev_filter];

            pdf(paste('degrade_mats_', this_exp, '.pdf',sep=''))
            for (i in which(longev_filter)) {
                
                only_data = na.omit(as.numeric(degrade_data[i,]));

                plot(only_data, ylim=c(min_diff, max_diff), xlab='Time (5 min)', ylab='Local Gel Diff', main=i)
                lines(c(0,130),c(0,0), col='red')
                
                for (j in 1:length(only_data)) {
                    if (j > 1) {
                        # tryCatch(conf_int = t.test(only_data[1:j])$conf.int,error=next);
                        conf_int = t.test(only_data[1:j])$conf.int;
                        errbar(j,mean(only_data[1:j]), conf_int[2], conf_int[1], col='green', add=T)
                    } else {
                        points(j, mean(only_data[1:j]), col='green',pch=19);
                    }
                }
            }
            graphics.off()
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
        exp_props = gather_invado_properties(data_dir);
        invado_lineage_nums = exp_props$longev_filter[exp_props$high_conf_int < 0];
        non_invado_lineage_nums = exp_props$longev_filter[exp_props$high_conf_int >= 0];


        # if (model_type == 'average') {
        #     average_model = gather_bilinear_models_from_dirs(data_dir,
        #             data_file='Average_adhesion_signal.csv', min_length = min_length,
        #             results.file=file.path('..','models','intensity.Rdata'), debug=debug)
        #         write_assembly_disassembly_periods(average_model[[1]],file.path(data_dir,'..'))	
        # }
        # if (model_type == 'cell_background') {
        #     temp = gather_bilinear_models_from_dirs(data_dir, 
        #             data_file='CB_corrected_signal.csv', min_length = min_length,
        #             results.file=file.path('..','models','CB_corrected.Rdata'), debug=debug)
        # }
        # if (model_type == 'local_background') {
        #     temp = gather_bilinear_models_from_dirs(data_dir, 
        #             data_file='Background_corrected_signal.csv', min_length = min_length,
        #             results.file=file.path('..','models','local_corrected.Rdata'), debug=debug)
        # }
        # if (model_type == 'area') {
        #     temp = gather_bilinear_models_from_dirs(data_dir, 
        #             data_file='Area.csv', min_length = min_length, log.trans = FALSE,
        #             results.file=file.path('..','models','area.Rdata'), debug=debug)
        # }
        # if (model_type == 'box_intensity') {
        #     temp = gather_bilinear_models_from_dirs(data_dir, 
        #             data_file='Box_intensity.csv', min_length = min_length,
        #             results.file=file.path('..','models','box.Rdata'), debug=TRUE)
        # }
        # if (model_type == 'background_correlation_model') {
        #     intensity_data <- 
        #         read.table(file.path(data_dir,'Background_corrected_signal.csv'), 
        #                 sep=',', header=FALSE);

        #     centroid_x <- read.table(file.path(data_dir,'Centroid_x.csv'), sep=',', header=FALSE);
        #     centroid_y <- read.table(file.path(data_dir,'Centroid_y.csv'), sep=',', header=FALSE);

        #     results = correlate_signal_vs_dist(intensity_data, centroid_x, centroid_y)
        #         output_file = file.path(data_dir,'..','models','background_corr.Rdata')
        #         if (! file.exists(dirname(output_file))) {
        #             dir.create(dirname(output_file),recursive=TRUE)
        #         }
        #     save(results,file = output_file);
        # }
    }
}

