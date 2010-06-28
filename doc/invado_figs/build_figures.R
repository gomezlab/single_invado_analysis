rm(list = ls())
source('../../src/analyze_cell_features/invado_analysis_lib.R')
source('../../fa_src/analyze_cell_features/FA_analysis_lib.R')
library(lattice)
library(geneplotter)
library(Hmisc)
library(boot)
source('errbar.s')

debug=TRUE;

################################################################################
# Data Reading/Initial Processing
################################################################################

dirs_to_load = list();
dirs_to_load$control = Sys.glob('../../results/Invadopodia/control/*/adhesion_props/');
dirs_to_load$coro1C_kd = Sys.glob('../../results/Invadopodia/Coro1C/*/adhesion_props/');
dirs_to_load$coro1B_kd = Sys.glob('../../results/Invadopodia/Coro1B/*/adhesion_props/');
dirs_to_load$CotL_kd = Sys.glob('../../results/Invadopodia/CotL/*/adhesion_props/');
dirs_to_load$Cort_kd = Sys.glob('../../results/Invadopodia/Cort/*/adhesion_props/');
dirs_to_load$BB94 = Sys.glob('../../results/Invadopodia/BB94/*/adhesion_props/');
# dirs_to_load$Exo70 = Sys.glob('../../results/Invadopodia/Exo70/*/adhesion_props/');
dirs_to_load$FAK = Sys.glob('../../results/Invadopodia/FAK/*/adhesion_props/');
dirs_to_load$LKB = Sys.glob('../../results/Invadopodia/LKB/*/adhesion_props/');

raw_data = list();
for (exp_type in names(dirs_to_load)) {
    raw_data[[exp_type]] = load_results_data_frame(dirs_to_load[[exp_type]], 
        file.path('models','puncta_props_corr.Rdata'), 'all_props',debug=F);
}

data_sets = list();
diff_conf = list();
combined_data_sets = list()
for (exp_name in names(raw_data)) {
    this_data_set = raw_data[[exp_name]];

    #Filter set of 0.95
    filter_sets = build_filter_sets(this_data_set);
    data_sets[[exp_name]]$invado = subset(this_data_set, filter_sets$invado_filter);
    data_sets[[exp_name]]$not_invado = subset(this_data_set, filter_sets$not_invado_filter);
    
    combined_data_sets$longevity = c(combined_data_sets$longevity, data_sets[[exp_name]]$invado$longevity);
    combined_data_sets$mean_area = c(combined_data_sets$mean_area, data_sets[[exp_name]]$invado$mean_area);
    combined_data_sets$mean_vals = c(combined_data_sets$mean_vals, data_sets[[exp_name]]$invado$mean_vals);
    combined_data_sets$exp_name = c(combined_data_sets$exp_name,
        rep(exp_name,length(data_sets[[exp_name]]$invado$longevity)));

    for (this_conf in as.character(seq(0.95,0.05,by=-0.05))) {
        # dealing with a strange bug in the as.character sequence above where the "0.1"
        # entry becomes "0.0999999999999999"
        if (this_conf == "0.0999999999999999") {
            this_conf = "0.1";
        }
        if (! is.list(diff_conf[[this_conf]])) {
            diff_conf[[this_conf]] = list();
        }
        filter_sets = build_filter_sets(this_data_set, conf.level=as.numeric(this_conf));
        diff_conf[[this_conf]][[exp_name]]$invado = subset(this_data_set, filter_sets$invado_filter);
        diff_conf[[this_conf]][[exp_name]]$not_invado = subset(this_data_set, filter_sets$not_invado_filter);
    }
}

combined_data_sets = as.data.frame(combined_data_sets);

# data_sets$day_1_control = subset(data_sets$invado_control, regexpr('11_16',data_sets$invado_control$experiment) == 1)
# data_sets$day_2_control = subset(data_sets$invado_control, regexpr('11_22',data_sets$invado_control$experiment) == 1)
# data_sets$day_3_control = subset(data_sets$invado_control, regexpr('02_23',data_sets$invado_control$experiment) == 1)
# 
# data_sets$day_1_coro1C = subset(data_sets$invado_coro1C_kd, regexpr('11_14',data_sets$invado_coro1C_kd$experiment) == 1)
# data_sets$day_2_coro1C = subset(data_sets$invado_coro1C_kd, regexpr('11_15',data_sets$invado_coro1C_kd$experiment) == 1)
# data_sets$day_3_coro1C = subset(data_sets$invado_coro1C_kd, regexpr('11_22',data_sets$invado_coro1C_kd$experiment) == 1)
# 
ts_props = list()

for (exp_type in names(dirs_to_load)) {
    ts_props[[exp_type]]$bleaching_curve = load_data_files(dirs_to_load[[exp_type]], c('bleaching_curves.csv'), debug=F, inc_exp_names=T);
}

################################################################################
# Processing
################################################################################

matrix_sets = list(
    mean = matrix(ncol = length(diff_conf[["0.95"]]), nrow = length(diff_conf)),
    yplus = matrix(ncol = length(diff_conf[["0.95"]]), nrow = length(diff_conf)),
    yminus = matrix(ncol = length(diff_conf[["0.95"]]), nrow = length(diff_conf))
);

area_comparison_sets = list(
    invado = matrix_sets,
    not_invado = matrix_sets
);

longevity_comparison_sets = area_comparison_sets;
mean_diff_sets = area_comparison_sets;

n_counts = list(
    invado = matrix_sets$mean,
    not_invado = matrix_sets$mean,
    percentage = matrix_sets$mean
)

exp_field_counts = list();

cutoff_sequence = c();
exp_name_sequence = c();

cutoff_count = 1;
for (cutoff in names(diff_conf)) {
    cutoff_sequence = c(cutoff_sequence, as.numeric(cutoff));

    print(paste('Starting on cutoff value:', cutoff))

    exp_count = 1;
    for (exp_name in names(diff_conf[[cutoff]])) {
        if (cutoff_count == 1) {
            exp_name_sequence = c(exp_name_sequence, exp_name);
        } else {
            stopifnot(exp_name_sequence[exp_count] == exp_name);
        }

        for (invado_class in c('invado','not_invado')) {
            #Area comparisons
            temp_props = gather_barplot_properties(diff_conf[[cutoff]][[exp_name]][[invado_class]]$mean_area, 
                bootstrap.rep = 100);
            area_comparison_sets[[invado_class]]$mean[cutoff_count, exp_count] = temp_props$mean;
            area_comparison_sets[[invado_class]]$yplus[cutoff_count, exp_count] = temp_props$yplus;
            area_comparison_sets[[invado_class]]$yminus[cutoff_count, exp_count] = temp_props$yminus;
            
            #Longevity comparisons
            temp_props = gather_barplot_properties(diff_conf[[cutoff]][[exp_name]][[invado_class]]$longevity*5,
                bootstrap.rep = 100);
            longevity_comparison_sets[[invado_class]]$mean[cutoff_count, exp_count] = temp_props$mean;
            longevity_comparison_sets[[invado_class]]$yplus[cutoff_count, exp_count] = temp_props$yplus;
            longevity_comparison_sets[[invado_class]]$yminus[cutoff_count, exp_count] = temp_props$yminus;
            
            #Mean Local Difference Comparisons
            temp_props = gather_barplot_properties(diff_conf[[cutoff]][[exp_name]][[invado_class]]$mean_vals,
                bootstrap.rep = 100);
            mean_diff_sets[[invado_class]]$mean[cutoff_count, exp_count] = temp_props$mean;
            mean_diff_sets[[invado_class]]$yplus[cutoff_count, exp_count] = temp_props$yplus;
            mean_diff_sets[[invado_class]]$yminus[cutoff_count, exp_count] = temp_props$yminus;

            n_counts[[invado_class]][cutoff_count, exp_count] = 
                length(diff_conf[[cutoff]][[exp_name]][[invado_class]]$mean_area);

            exp_field_counts[[exp_name]] = length(dirs_to_load[[exp_name]]);
        }
        exp_count = exp_count + 1;
    }
    cutoff_count = cutoff_count + 1;
}

for (invado_class in c('invado','not_invado')) {
    for (type in names(longevity_comparison_sets[[invado_class]])) {
        colnames(longevity_comparison_sets[[invado_class]][[type]]) <- exp_name_sequence;
        rownames(longevity_comparison_sets[[invado_class]][[type]]) <- cutoff_sequence;
        
        colnames(area_comparison_sets[[invado_class]][[type]]) <- exp_name_sequence;
        rownames(area_comparison_sets[[invado_class]][[type]]) <- cutoff_sequence;
    }
}

n_counts$percentage = n_counts$invado/(n_counts$not_invado+n_counts$invado);
stop('Done with Loading/Processing')

################################################################################
# Plotting
################################################################################

####################
# Area Comparisons
####################

dir.create('area',recursive=TRUE, showWarnings=FALSE);

area_names = c()
for (i in 1:length(colnames(area_comparison_sets$invado$mean))) {
    area_names = c(area_names, paste(colnames(area_comparison_sets$invado$mean)[[i]], " (",n_counts$invado[1,i], ")", sep=''))
}

svg(file.path('area','area_barplots.svg'), width=14)
par(mar=c(2,4,0.5,0))
x_pos = barplot(area_comparison_sets$invado$mean[1,], names=area_names,
        ylab='Mean Puncta Area (\u03BCm\u00B2)', ylim=c(0,max(area_comparison_sets$invado$yplus[1,])))
errbar(t(x_pos),area_comparison_sets$invado$mean[1,],
       area_comparison_sets$invado$yplus[1,], area_comparison_sets$invado$yminus[1,],
       add=TRUE,cex=1E-10,lwd=1.5)
# segments(-1,area_comparison_sets$invado$yminus[1,1],20,area_comparison_sets$invado$yminus[1,1])
# segments(-1,area_comparison_sets$invado$yplus[1,1],20,area_comparison_sets$invado$yplus[1,1])
graphics.off()

####################
# Longevity Comparisons
####################
dir.create('longevity',recursive=TRUE, showWarnings=FALSE);

svg(file.path('longevity','longevity_barplots.svg'), width=14)
par(mar=c(2,4,0.5,0))
x_pos = barplot(longevity_comparison_sets$invado$mean[1,], names=area_names,
        ylab='Puncta Longevity (min)', ylim=c(0,max(longevity_comparison_sets$invado$yplus[1,])))
errbar(t(x_pos),longevity_comparison_sets$invado$mean[1,],
       longevity_comparison_sets$invado$yplus[1,], longevity_comparison_sets$invado$yminus[1,],
       add=TRUE,cex=1E-10,lwd=1.5)
# segments(-1,longevity_comparison_sets$invado$yminus[1,1],20,longevity_comparison_sets$invado$yminus[1,1])
# segments(-1,longevity_comparison_sets$invado$yplus[1,1],20,longevity_comparison_sets$invado$yplus[1,1])
graphics.off()


####################
# Mean Local Diff Comparisons
####################

dir.create('mean_local',recursive=TRUE, showWarnings=FALSE);
svg(file.path('mean_local','local_diff_means.svg'), width=14)
par(mar=c(2,4,0.5,0))
x_pos = barplot(mean_diff_sets$invado$mean[1,], names=area_names,
       ylab='Mean Degradation Level', ylim=c(min(mean_diff_sets$invado$yminus[1,]), 0))
errbar(t(x_pos),mean_diff_sets$invado$mean[1,],
       mean_diff_sets$invado$yplus[1,], mean_diff_sets$invado$yminus[1,],
       add=TRUE,cex=1E-10,lwd=1.5)
# segments(-1,mean_diff_sets$invado$yminus[1,1],20,mean_diff_sets$invado$yminus[1,1])
# segments(-1,mean_diff_sets$invado$yplus[1,1],20,mean_diff_sets$invado$yplus[1,1])
graphics.off();

all_mean_data = list();
all_mean_data$area = area_comparison_sets$invado$mean[1,];
all_mean_data$longevity = longevity_comparison_sets$invado$mean[1,];
all_mean_data$mean_diff = mean_diff_sets$invado$mean[1,];
splom(area * longevity * mean_diff, data=all_mean_data)

####################
# Area p-value Comparisons - line plot
####################

svg(file.path('area','area_p_value.svg'))
layout(matrix(1:6,nrow=3,ncol=2))
par(mar=c(4.5,4,0,0), bty='n')

for (i in 1:dim(area_comparison_sets$invado$mean)[2]) {
    invado_data = list();
    invado_data$mean = area_comparison_sets$invado$mean[,i];
    invado_data$yplus = area_comparison_sets$invado$yplus[,i];
    invado_data$yminus = area_comparison_sets$invado$yminus[,i];
    
    not_invado_data = list();
    not_invado_data$mean = area_comparison_sets$invado$mean[,1];
    not_invado_data$yplus = area_comparison_sets$invado$yplus[,1];
    not_invado_data$yminus = area_comparison_sets$invado$yminus[,1];

    min_max = range(invado_data, not_invado_data);

#Plotting the Invadopodia at each cutoff value
    plot(cutoff_sequence, invado_data$mean, type='l', ylim=min_max, #xlim=c(0,1), 
        col='red', xlab='P-value Cutoff',ylab='Mean Puncta Area (\u03BCm\u00B2)', lwd=3);
    errbar(t(cutoff_sequence),invado_data$mean, 
        invado_data$yplus, invado_data$yminus,
        add=TRUE,cex=1E-10, col='red', lwd=3)

#Plotting the Non-invadopodia at each cutoff value
    lines(cutoff_sequence, not_invado_data$mean, col='blue', lwd=3);
    errbar(t(cutoff_sequence),not_invado_data$mean, 
        not_invado_data$yplus, not_invado_data$yminus,
        add=TRUE,cex=1E-10, col='blue', lwd=3)
    
    plot_dims = par("usr");
    if (i == 1) {
        legend('topleft',c('Invadopodia','Not Invadopodia'), fill=c('red','blue'))
        text(0.05*(plot_dims[2] - plot_dims[1])+plot_dims[1],
            (plot_dims[4]-plot_dims[3])/2+plot_dims[3],
            exp_name_sequence[i],adj=c(0,0.5));
    } else {
        text(0.5*(plot_dims[2] - plot_dims[1])+plot_dims[1],
            plot_dims[4],
            exp_name_sequence[i],adj=c(0.5,1));
    }
}

graphics.off()

####################
# Longevity p-value Comparisons - line plot
####################

svg(file.path('longevity','longevity_p_value.svg'))
layout(matrix(1:6,nrow=3,ncol=2))
par(mar=c(4.5,4,0,0), bty='n')

for (i in 1:dim(longevity_comparison_sets$invado$mean)[2]) {
    invado_data = list();
    invado_data$mean = longevity_comparison_sets$invado$mean[,i];
    invado_data$yplus = longevity_comparison_sets$invado$yplus[,i];
    invado_data$yminus = longevity_comparison_sets$invado$yminus[,i];
    
    not_invado_data = list();
    not_invado_data$mean = longevity_comparison_sets$invado$mean[,1];
    not_invado_data$yplus = longevity_comparison_sets$invado$yplus[,1];
    not_invado_data$yminus = longevity_comparison_sets$invado$yminus[,1];

    min_max = range(invado_data, not_invado_data);

#Plotting the Invadopodia at each cutoff value
    plot(cutoff_sequence, invado_data$mean, type='l', ylim=min_max, #xlim=c(0,1), 
        col='red', xlab='P-value Cutoff',ylab='Longevity (min)', lwd=3);
    errbar(t(cutoff_sequence),invado_data$mean, 
        invado_data$yplus, invado_data$yminus,
        add=TRUE,cex=1E-10, col='red', lwd=3)

#Plotting the Non-invadopodia at each cutoff value
    lines(cutoff_sequence, not_invado_data$mean, col='blue', lwd=3);
    errbar(t(cutoff_sequence),not_invado_data$mean, 
        not_invado_data$yplus, not_invado_data$yminus,
        add=TRUE,cex=1E-10, col='blue', lwd=3)
    
    plot_dims = par("usr");
    if (i == 1) {
        legend('topleft',c('Invadopodia','Not Invadopodia'), fill=c('red','blue'))
        text(0.05*(plot_dims[2] - plot_dims[1])+plot_dims[1],
            (plot_dims[4]-plot_dims[3])/2+plot_dims[3],
            exp_name_sequence[i],adj=c(0,0.5));
    } else {
        text(0.5*(plot_dims[2] - plot_dims[1])+plot_dims[1],
            plot_dims[4],
            exp_name_sequence[i],adj=c(0.5,1));
    }
}

graphics.off()


####################
# Longevity/Area Comparisons
####################
dir.create('longevity_size_comp',recursive=TRUE, showWarnings=FALSE);

data_names = c('Control','Coro1B KD', 'Coro1C KD', 'Coactosin KD')

for (i in 1:length(longev_sets_to_use)) {
    svg(file.path('longevity_size_comp',paste(data_names[[i]], '.svg', sep='')))
    par(bty='n');
    plot(longev_sets_to_use[[i]], area_sets_to_use[[i]], main=data_names[[i]],
        xlab='Longevity (min)', ylab='Mean Puncta Area (\u03BCm\u00B2)')
    graphics.off()
}

svg(file.path('longevity_size_comp','all.svg'))
par(bty='n');
layout(rbind(c(1,2),c(3,4)))
for (i in 1:length(longev_sets_to_use)) {
    plot(longev_sets_to_use[[i]], area_sets_to_use[[i]], main=data_names[[i]],
        xlab='Longevity (min)', ylab='Mean Puncta Area (\u03BCm\u00B2)')
}
graphics.off()
