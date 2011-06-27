################################################################################
# invado_analysis_lib.R: the R functions associated with processing the
# invaodopodia data
################################################################################

library(Hmisc);

gather_invado_properties <- function(results_dirs, build_degrade_plots = FALSE, 
    conf.level = 0.95, degrade_file = "Local_gel_diff.csv", results.file = NA, 
    build_plots=TRUE, debug=FALSE) {
    
    for (this_exp_dir in results_dirs) {
        all_props = list();
        
        if (! file.exists(file.path(this_exp_dir,'single_lin.csv'))) {
            warning(paste("Unable to find", file.path(this_exp_dir,'single_lin.csv'),
                  "moving on to the next file."))
            next;
        }
        
        #Reading in raw data
        print(paste('Working on:',this_exp_dir));
        lineage_data = read.table(file.path(this_exp_dir,'single_lin.csv'), 
            sep=",",header=T);
        local_diff_data = read.table(file.path(this_exp_dir,'lin_time_series', degrade_file), 
            sep=",",header=F);
        pre_diff_data = read.table(file.path(this_exp_dir,'lin_time_series', 'Pre_birth_diff.csv'), 
            sep=",",header=F);
        area_data = read.table(file.path(this_exp_dir,'lin_time_series', 'Area.csv'), 
            sep=",",header=F);
        edge_dist_data = read.table(file.path(this_exp_dir,'lin_time_series', 'Centroid_dist_from_edge.csv'), 
            sep=",",header=F);
        stopifnot(dim(lineage_data)[[1]] == dim(local_diff_data)[[1]])
        
        #Building the filter sets
        
        #exclude puncta that don't live for at least 5 time steps
        longev_filter = ! is.na(lineage_data$longevity) & lineage_data$longevity >= 5;
        
        #exclude puncta that underwent split births
        no_split_birth_filt = ! is.na(lineage_data$split_birth_status) & ! lineage_data$split_birth_status;

        #exclude puncta that disappeared due to merge events (coded as 0 in the
        #lineage data file)
        death_filt = ! is.na(lineage_data$death_status) & lineage_data$death_status;

        overall_filt = longev_filter & no_split_birth_filt & death_filt;

        all_props$lineage_nums = c(all_props$lineage_nums,which(overall_filt));
        all_props$experiment = c(all_props$experiment,
            rep(basename(dirname(this_exp_dir)),length(which(overall_filt))));
        props_to_include = c("longevity","largest_area","birth_i_num","mean_area");
        for (i in props_to_include) {
            all_props[[i]] = c(all_props[[i]], lineage_data[[i]][overall_filt]);
        }

        if (build_plots) {
            pdf(file.path(this_exp_dir,'local_degrade_plots.pdf'));
        }

        #analyzing each of the puncta in the filtered set for invadopodia filtering processes
        for (lin_num in which(overall_filt)) {
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
            
            #Pre-birth local difference data
            pre_diff = na.omit(as.numeric(pre_diff_data[lin_num,]));
            
            stat_tests$pre_diff = tryCatch(
                t.test(pre_diff,conf.level=conf.level), 
                error = t.test.error);

            stat_tests$pre_local_diff = tryCatch(
                t.test(local_diff - pre_diff,conf.level=conf.level), 
                error = t.test.error);
            all_props$pre_diff_p_value = c(all_props$pre_diff_p_value, stat_tests$pre_local_diff$p.value);
            all_props$mean_pre_diff = c(all_props$mean_pre_diff, as.numeric(stat_tests$pre_local_diff$estimate));

            only_area_data = na.omit(as.numeric(area_data[lin_num,]));
            only_edge_dist_data = na.omit(as.numeric(edge_dist_data[lin_num,]));

            all_props$mean_edge_dist = c(all_props$mean_edge_dist, mean(only_edge_dist_data));
            
            # if (lin_num == 515) {
            #     all_three_sets = cbind(local_diff, pre_diff, local_diff - pre_diff);
            #     # browser();
            #     svg('single_invado_ctrl_02_23_pos04_515.svg')
            #     build_presentation_single_invado_plot(all_three_sets,stat_tests,lin_num);
            #     graphics.off()
            # }

            if (build_plots) {
                all_three_sets = cbind(local_diff, pre_diff, local_diff - pre_diff);
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
    errbar(max(time_points)*1.05, stat_tests$pre_local_diff$estimate, 
        stat_tests$pre_local_diff$conf.int[2], stat_tests$pre_local_diff$conf.int[1], 
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

build_presentation_single_invado_plot <- function(data_sets,stat_tests, lin_num) {
    time_points = seq(from=0,by=5,along.with=data_sets[,1]);

    par(bty='n',mar=c(3.25,2.83,0.5,0),cex=2,mgp=c(2.1,1,0))
    matplot(time_points, data_sets, 
        typ='l', lty=c(1,2,4), xlab='Time since Invadopodia Birth (min)', ylab='Fluorescence Difference', 
        lwd=5, xlim=c(0,max(time_points)*1.05), ylim=c(1.05*min(data_sets), 3.25*max(data_sets)))
    
    plot_limits = par("usr");

    #place the legend
    x_pos = (plot_limits[2] - plot_limits[1])*0.075 + plot_limits[1];
    y_pos = (plot_limits[4] - plot_limits[3])*1 + plot_limits[3];
    legend(x_pos,y_pos,c('Local','Pre-birth', 'Local - Pre-birth' ), 
        col=c('black','red', 'green'), lty=c(1,2,4), lwd=3)
    
    #draw the zero line
    segments(0,0,max(time_points),0,lty=4)
    
    #put in the error bars
    par(cex = 0.5, lwd=2);
    errbar(max(time_points)*1.01, stat_tests$local_diff$estimate, 
        stat_tests$local_diff$conf.int[2], stat_tests$local_diff$conf.int[1], add=T,lwd=3)
    errbar(max(time_points)*1.03, stat_tests$pre_diff$estimate, 
        stat_tests$pre_diff$conf.int[2], stat_tests$pre_diff$conf.int[1], add=T, col='red',lwd=3)
    errbar(max(time_points)*1.05, stat_tests$pre_local_diff$estimate, 
        stat_tests$pre_local_diff$conf.int[2], stat_tests$pre_local_diff$conf.int[1], 
        add=T, col='green',lwd=3)
}
 
t.test.error <- function(e) {
    list(conf.int = c(Inf, -Inf), p.value = 1)
}

build_filter_sets <- function(raw_data_set, conf.level = 0.95) {
    filter_sets = list();

    # filter_sets$local_diff_filter = raw_data_set$high_conf_int < 0;
    filter_sets$pre_diff_filter = raw_data_set$mean_pre_diff < 0 & raw_data_set$pre_diff_p_value < (1 - conf.level);
    filter_sets$local_diff_filter = raw_data_set$p_value < (1 - conf.level);
    # filter_sets$pre_diff_filter = raw_data_set$pre_diff_p_value < 1 - conf.level;
    
    filter_sets$invado_filter = filter_sets$local_diff_filter & filter_sets$pre_diff_filter;
    filter_sets$not_invado_filter = ! filter_sets$invado_filter;

    return(filter_sets);
}

push <- function(data_set, new_data) {
    if (is.list(data_set)) {
        curr_length = length(data_set);
        new_data_pos = curr_length + 1;
        data_set[[new_data_pos]] = new_data;
    } else {
        data_set = c(data_set, new_data);
    }
    return(data_set);
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
