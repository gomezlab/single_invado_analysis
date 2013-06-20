###############################################################################
# Data Loading
###############################################################################

load_results_data_frame <- function(dirs,file,variable_name, debug=FALSE) {
    #simple function to read in all the R data files in each of the specified
    #dirs, searching for the provided variable name after every load, note that
    #all the data is concatenated into a single data frame
	
    results = data.frame()
	for (i in 1:length(dirs)) {
		this_file = file.path(dirs[i],file)
        if (debug) {
            print(paste("Loading:", file.path(dirs[i],file)))
        }
		if (file.exists(this_file)) {
			load(file.path(dirs[i],file))
            browser()
			results = rbind(results,get(variable_name))
		}
	}
	results
}

###############################################################################
# Data Summary
###############################################################################

gather_barplot_properties <- function(data_sets, bootstrap.rep = 50000) {
    plot_props = list()
    if (is.numeric(data_sets)) {
        plot_props$mean = mean(data_sets);
        temp_conf_data = determine_mean_conf_int(data_sets, bootstrap.rep);

        plot_props$yminus = temp_conf_data[1];
        plot_props$yplus = temp_conf_data[2];
    } else {
        for (i in 1:length(data_sets)) {
            plot_props$mean = c(plot_props$mean, mean(data_sets[[i]]));

            temp_conf_data = determine_mean_conf_int(data_sets[[i]], bootstrap.rep);

            plot_props$yminus = c(plot_props$yminus, temp_conf_data[1]);
            plot_props$yplus = c(plot_props$yplus, temp_conf_data[2]);
        }
    }
    return(plot_props);
}

determine_mean_conf_int <- function(data, bootstrap.rep = 50000) {
	require(boot);
	boot_samp = boot(data, function(data,indexes) mean(data[indexes],na.rm=T), bootstrap.rep);
		
    conf_int = boot.ci(boot_samp, type="perc", conf=0.95)$percent[4:5]

    return(conf_int)
}
