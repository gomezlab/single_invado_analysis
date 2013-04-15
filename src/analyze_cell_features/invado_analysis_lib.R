################################################################################
# invado_analysis_lib.R: the R functions associated with processing the
# invaodopodia data
################################################################################

library(Hmisc);

gather_invado_properties <- function(results_dirs, build_degrade_plots = FALSE,
                                     conf.level = 0.95, results.file = NA,
                                     build_plots=TRUE, debug=FALSE) {
    
    for (this_exp_dir in results_dirs) {
        all_props = list();
        
        ########################################################################
        #Reading in raw data
        ########################################################################
        print(paste('Working on:',this_exp_dir));
        data_folder = file.path(this_exp_dir,'lin_time_series');
        area_data = read.table(file.path(data_folder, 'Area.csv'), 
            sep=",",header=F);
        local_diff_data = read.table(file.path(data_folder,'Local_gel_diff.csv'), 
            sep=",",header=F);
        local_diff_corrected_data = read.table(file.path(data_folder,'Local_diff_corrected.csv'), 
            sep=",",header=F);
        pre_diff_data = read.table(file.path(data_folder,'Pre_birth_diff.csv'), 
            sep=",",header=F);
        edge_dist_data = read.table(file.path(data_folder,'Centroid_dist_from_edge.csv'), 
            sep=",",header=F);

        ########################################################################
        #Building the filter sets
        ########################################################################
        longevity_uncertain = rowSums(! is.na(area_data))

        #We only want to consider puncta that live for at least 5 time steps
        longev_filter = ! is.na(longevity_uncertain) & longevity_uncertain >= 5;
        
        overall_filt = longev_filter;
        
        if (sum(overall_filt) == 0) {
            print("None of the lineages passed the filtering step, returning from function with no output.");
            return(all_props);
        }
        
        all_props$birth_observed = is.na(area_data[overall_filt,1]);
        all_props$death_observed = is.na(area_data[overall_filt,dim(area_data)[2]]);
        
        longevity = longevity_uncertain;
        longevity[! (all_props$birth_observed & all_props$death_observed)] = NA;

        all_props$lineage_nums = which(overall_filt)
        all_props$experiment = rep(this_exp_dir,sum(overall_filt));
        all_props$longevity_uncertain = longevity_uncertain[overall_filt]
        all_props$longevity = longevity[overall_filt]

        all_props$mean_area = rowMeans(area_data[overall_filt,],na.rm=T);
        all_props$mean_edge_dist = rowMeans(edge_dist_data[overall_filt,],na.rm=T);
        
        if (build_plots) {
            pdf(file.path(this_exp_dir,'local_degrade_plots.pdf'));
        }
        #analyzing each of the puncta in the filtered set to identify invadopodia
        for (lin_num in which(overall_filt)) {
            ###############################################
            # Local Diff
            ###############################################
            local_diff = na.omit(as.numeric(local_diff_data[lin_num,]));

            all_props$mean_local_diff = c(all_props$mean_local_diff, mean(local_diff))
            all_props$max_local_diff = c(all_props$max_local_diff, max(local_diff));
            all_props$min_local_diff = c(all_props$min_local_diff, min(local_diff));
            
            stat_tests = list();

            stat_tests$local_diff = tryCatch(
                t.test(local_diff,conf.level=conf.level), 
                error = t.test.error);
            
            all_props$low_conf_int = c(all_props$low_conf_int, stat_tests$local_diff$conf.int[1]);
            all_props$high_conf_int = c(all_props$high_conf_int, stat_tests$local_diff$conf.int[2]);
            all_props$p_value = c(all_props$p_value, stat_tests$local_diff$p.value);
            
            ###############################################
            # Pre-birth
            ###############################################
            pre_diff = na.omit(as.numeric(pre_diff_data[lin_num,]));
            
            stat_tests$pre_diff = tryCatch(
                t.test(pre_diff,conf.level=conf.level), 
                error = t.test.error);
            
            all_props$mean_pre_diff = c(all_props$mean_pre_diff,mean(pre_diff));
            all_props$mean_pre_diff_p_value = c(all_props$mean_pre_diff_p_value,
                                                stat_tests$pre_diff$p.value);

            ###############################################
            # Local Diff Corrected
            ###############################################
            local_diff_corrected = na.omit(as.numeric(local_diff_corrected_data[lin_num,]));
            
            stat_tests$local_diff_corrected = tryCatch(
                t.test(local_diff_corrected,conf.level=conf.level), 
                error = t.test.error);

            all_props$local_diff_corrected_p_value = c(all_props$local_diff_corrected_p_value, 
                stat_tests$local_diff_corrected$p.value);
            all_props$mean_local_diff_corrected = c(all_props$mean_local_diff_corrected, 
                as.numeric(stat_tests$local_diff_corrected$estimate));
            
            ###############################################
            # Local Diff Corrected vs Pre Diff
            ###############################################
            stat_tests$local_diff_pre_diff = tryCatch(
                t.test(pre_diff,local_diff_corrected), 
                error = t.test.error);
            
            all_props$local_pre_p_value = c(all_props$local_pre_p_value, 
                stat_tests$local_diff_pre_diff$p.value);

            if (length(all_props$mean_local_diff_corrected) != length(all_props$mean_local_diff)) {
                # browser()
            }

            if (build_plots) {
                all_three_sets = cbind(local_diff, pre_diff, local_diff_corrected);
                # browser();
                build_single_invado_plot(all_three_sets,stat_tests,lin_num);
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

build_single_invado_plot <- function(data_sets,stat_tests, lin_num) {
    time_points = seq(from=0,by=5,along.with=data_sets[,1]);
    
    par(bty='n', mar=c(4,4,2,0))
    matplot(time_points, data_sets, 
        typ='l', lty=c(1,2,4), xlab='Time (min)', ylab='Difference Metric', main=lin_num, 
        lwd=2, xlim=c(0,max(time_points)*1.05))
    
    plot_limits = par("usr");
    
    legend('topleft',c('Local Diff','Pre-birth Local Diff', 'Local Diff - Pre-birth Diff' ), 
        fill=c('black','red', 'green'))
    segments(0,0,max(time_points),0,lty=4)

    errbar(max(time_points)*1.01, stat_tests$local_diff$estimate, 
        stat_tests$local_diff$conf.int[2], stat_tests$local_diff$conf.int[1], add=T)
    errbar(max(time_points)*1.03, stat_tests$pre_diff$estimate, 
        stat_tests$pre_diff$conf.int[2], stat_tests$pre_diff$conf.int[1], add=T, col='red')
    errbar(max(time_points)*1.05, stat_tests$local_diff_corrected$estimate, 
        stat_tests$local_diff_corrected$conf.int[2], stat_tests$local_diff_corrected$conf.int[1], 
        add=T, col='green')

    #Adding the areas to the same plot
    # plot_props = par('usr');
    # 
    # scale_factor = (plot_props[4]-plot_props[3])/max(only_area_data)
    # lines(time_points,only_area_data*scale_factor+plot_props[3],col='purple')
    # 
    # tick_labels = (axTicks(2) - plot_props[3])/scale_factor;

    # axis(4,at=axTicks(2),labels=sprintf('%.3f',(axTicks(2) - plot_props[3])/scale_factor))
}

t.test.error <- function(e) {
    list(conf.int = c(Inf, -Inf), p.value = 1)
}

build_filter_sets <- function(raw_data_set, conf.level = 0.99,min_mean_local_diff_corrected = NA) {
    filter_sets = list();

    filter_sets$local_diff_filter = raw_data_set$mean_local_diff > 0 &
        raw_data_set$p_value < (1 - conf.level);
    
    filter_sets$local_pre_p_value = raw_data_set$local_pre_p_value < (1 - conf.level);

    filter_sets$pre_diff_filter = ! raw_data_set$birth_observed |
        (raw_data_set$mean_local_diff_corrected > 0 & 
         raw_data_set$local_diff_corrected_p_value < (1 - conf.level));

    filter_sets$invado_filter = filter_sets$local_diff_filter & filter_sets$pre_diff_filter & 
        filter_sets$local_pre_p_value;

    if (!is.na(min_mean_local_diff_corrected)) {
        filter_sets$min_local_diff_corrected = 
            raw_data_set$mean_local_diff_corrected > min_mean_local_diff_corrected;
        filter_sets$invado_filter = filter_sets$invado_filter & filter_sets$min_local_diff_corrected;
    }

    filter_sets$not_invado_filter = ! filter_sets$invado_filter;

    return(filter_sets);
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
        
        if (dim(exp_props)[1] == 0) {
            print("Didn't find any lineages to analyze, probably nothing long-lived enough.");
            return;
        }
		
		data_types_to_include = c('lineage_nums', 'p_value', 'mean_local_diff',
								  'local_diff_corrected_p_value',
								  'mean_local_diff_corrected');
        
        filter_sets = build_filter_sets(exp_props);
        
		invado_lineage_data = subset(exp_props, filter_sets$invado_filter,
									 select = data_types_to_include);
        
        not_invado_lineage_data = subset(exp_props, filter_sets$not_invado_filter, 
            select = data_types_to_include);
        
        write.table(invado_lineage_data, file.path(data_dir, 'invado_data.csv'), 
            row.names=F, col.names=T, sep=',')
        write.table(not_invado_lineage_data, file.path(data_dir, 'not_invado_data.csv'), 
            row.names=F, col.names=T, sep=',')
    }
}
