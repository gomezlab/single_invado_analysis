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
        no_split_birth_filt = ! is.na(lineage_data$split_birth_status) & ! lineage_data$split_birth_status;
        death_filt = ! is.na(lineage_data$death_status) & lineage_data$death_status;
        eccentricity_filt = ! is.na(lineage_data$average_eccentricity) & lineage_data$average_eccentricity > 3;

        overall_filt = longev_filter & no_split_birth_filt & death_filt & eccentricity_filt;

        all_props$overall_filt = c(all_props$overall_filt,which(overall_filt));
        all_props$longevity = c(all_props$longevity, lineage_data$longevity[overall_filt]);

        final_conf_ints = matrix(ncol=2);
        for (i in which(overall_filt)) {
            only_data = na.omit(as.numeric(degrade_data[i,]));
            conf_int = t.test(only_data,conf.level=conf.level)$conf.int;
            final_conf_ints = rbind(final_conf_ints, conf_int[1:2]); 
        }
        final_conf_ints = final_conf_ints[2:dim(final_conf_ints)[[1]],];
        
        all_props$low_conf_int = c(all_props$low_conf_int, final_conf_ints[,1]);
        all_props$high_conf_int = c(all_props$high_conf_int, final_conf_ints[,2]);

        if (build_degrade_plots) {
            min_diff = min(degrade_data[which(overall_filt),], na.rm=T);
            max_diff = max(degrade_data[which(overall_filt),], na.rm=T);
            
            filtered_longev = lineage_data$longevity[overall_filt];

            pdf(paste('degrade_mats_', this_exp, '.pdf',sep=''))
            for (i in which(overall_filt)) {
                
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
        invado_lineage_nums = exp_props$overall_filt[exp_props$high_conf_int < 0];
        non_invado_lineage_nums = exp_props$overall_filt[exp_props$high_conf_int >= 0];
        
        write.table(t(invado_lineage_nums), file.path(data_dir, 'invado_nums.csv'), row.names=F, col.names=F, sep=',')
        write.table(t(non_invado_lineage_nums), file.path(data_dir, 'non_invado_nums.csv'), row.names=F, col.names=F, sep=',')
    }
}

