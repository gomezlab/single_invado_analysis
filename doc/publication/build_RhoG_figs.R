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
par(mar=c(4,4.5,1,0), bty='n');
boxplot_with_points(list(processed$only_signif$KD$assembly$slope, 
                         processed$only_signif$control$assembly$slope), 
                    names=c("KD", "Control"), inc.points=FALSE,
                    ylab=expression(paste('Assembly Rate (',min^-1,')',sep='')))
boxplot_with_points(list(processed$only_signif$KD$disassembly$slope, 
                         processed$only_signif$control$disassembly$slope), 
                    names=c("KD", "Control"), inc.points=FALSE,
                    ylab=expression(paste('Disassembly Rate (',min^-1,')',sep='')))
graphics.off()

for (i in unique(processed$only_signif$KD$assembly$exp_num)) {
    for (j in unique(processed$only_signif$KD$assembly$exp_num)) {
        if (j == i) {
            next;
        }
        if (i > j) {
            next;
        }
        filtered_data = subset(processed$only_signif$KD$assembly, exp_num != i | exp_num != j);
        filtered_dis = subset(processed$only_signif$KD$disassembly, exp_num != i | exp_num != j);
        
        print(paste(i,j))

        p_vals = determine_median_p_value(filtered_data$slope, processed$only_signif$control$assembly$slope);
        print(p_vals$ratio_conf)
        p_vals = determine_median_p_value(filtered_dis$slope, processed$only_signif$control$disassembly$slope);
        print(p_vals$ratio_conf)


    }
}

svg(file.path(out_folder, 'RhoG', 'single_exp.svg'));
layout(rbind(c(1,2),c(3,4)))
par(mar=c(4,4.5,1,0), bty='n');
kd_exp_assem = boxplot(slope ~ exp_num, data=processed$only_signif$KD$assembly, 
    notch=T, varwidth=T,ylim=c(0,0.12), xlab='KD Experiments',
    ylab=expression(paste('Assembly Rate (',min^-1,')')), pch=19,cex=0.25)
boxplot(slope ~ exp_num, data=processed$only_signif$control$assembly, 
    notch=T, varwidth=T, ylim=c(0,0.12), xlab='Control Experiments',
    ylab=expression(paste('Assembly Rate (',min^-1,')')), pch=19,cex=0.25)


kd_exp_disassem = boxplot(slope ~ exp_num, data=processed$only_signif$KD$disassembly, 
    notch=T, varwidth=T, ylim=c(0,0.12), xlab='KD Experiments',
    ylab=expression(paste('Disassembly Rate (',min^-1,')')), pch=19,cex=0.25)
boxplot(slope ~ exp_num, data=processed$only_signif$control$disassembly, 
    notch=T, varwidth=T, ylim=c(0,0.12),xlab='Control Experiments',
    ylab=expression(paste('Disassembly Rate (',min^-1,')')), pch=19,cex=0.25)

graphics.off()
