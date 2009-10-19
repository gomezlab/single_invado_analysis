rm(list = ls())
source('FA_analysis_lib.R')
library(lattice)
library(geneplotter)
library(Hmisc)
debug = TRUE;

################################################################################
#Result loading
################################################################################
raw_data <- list()

exp_dirs <- Sys.glob('../../results/FA_RhoG/*')
exp_dirs <- exp_dirs[file_test('-d',exp_dirs)]

for (exp_type_dir in exp_dirs) {
    exp_name = basename(exp_type_dir)
    model_folder_set = Sys.glob(file.path(exp_type_dir,'*/*/models'))
    raw_data[[exp_name]] = load_results(model_folder_set,file.path('intensity.Rdata'), debug=TRUE)
}

print('Done Loading Data')
########################################
#Result filtering
########################################
processed = list();

for (exp_name in names(raw_data)) {
        if (debug) {
            print(paste("Filtering", exp_name));
        }
        processed$no_filt[[exp_name]] = filter_results(raw_data[[exp_name]], 
            min_R_sq = -Inf, max_p_val = Inf);
        
        processed$only_signif[[exp_name]] = filter_results(raw_data[[exp_name]], 
            min_R_sq = -Inf, max_p_val = 0.05);
        
        processed$high_Rsq[[exp_name]] = filter_results(raw_data[[exp_name]]);
}

#rm(raw_data)
gc()

print('Done Filtering Data')
out_folder = '../../doc/publication/figures'
dir.create(out_folder,recursive=TRUE, showWarnings=FALSE);
stop();
################################################################################
#Plotting
################################################################################

dir.create(dirname(file.path(out_folder, 'RhoG', 'kinetics_comparison.svg')), 
    recursive=TRUE, showWarnings=FALSE);

svg(file.path(out_folder, 'RhoG', 'kinetics_comparison.svg'), height=7*(1/2));
layout(rbind(c(1,2)))
par(mar=c(4,4,1,0), bty='n');
boxplot_with_points(list(processed$only_signif$KD$assembly$slope, 
                         processed$only_signif$control$assembly$slope), 
                    names=c("KD", "Control"), inc.points=FALSE,
                    ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')))
boxplot_with_points(list(processed$only_signif$KD$disassembly$slope, 
                         processed$only_signif$control$disassembly$slope), 
                    names=c("KD", "Control"), inc.points=FALSE,
                    ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')))
graphics.off()
