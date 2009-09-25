rm(list = ls());
source('FA_analysis_lib.R');
library(lattice);
library(geneplotter);
library(Hmisc);

################################################################################
#Result loading
################################################################################
raw_data <- list()

#Rapamycin results
exp_dirs_rap <- Sys.glob('../../results/rap_src/*/adhesion_props/models/');
exp_dirs_rap <- exp_dirs_rap[file_test('-d',exp_dirs_rap)];
raw_data$rap = load_results(exp_dirs_rap,file.path('intensity.Rdata'));

exp_dirs_rap <- Sys.glob('../../results/rap_src_control/*/adhesion_props/models/');
exp_dirs_rap <- exp_dirs_rap[file_test('-d',exp_dirs_rap)];
raw_data$rap_ctrl = load_results(exp_dirs_rap,file.path('intensity.Rdata'));

exp_dirs_rap <- Sys.glob('../../results/rap_src_pax/*/adhesion_props/models/');
exp_dirs_rap <- exp_dirs_rap[file_test('-d',exp_dirs_rap)];
raw_data$rap_pax = load_results(exp_dirs_rap,file.path('intensity.Rdata'));

print('Done Loading Data')

################################################################################
#Processing
################################################################################

results_rap_onlysignif = filter_results(raw_data$rap, min_R_sq = -Inf, max_p_val = 0.05);
results_rap_ctrl_onlysignif = filter_results(raw_data$rap_ctrl, min_R_sq = -Inf, max_p_val = 0.05);
results_rap_pax_onlysignif = filter_results(raw_data$rap_pax, min_R_sq = -Inf, max_p_val = 0.05);

results_rap = filter_results(raw_data$rap);
results_rap_ctrl = filter_results(raw_data$rap_ctrl);
results_rap_pax = filter_results(raw_data$rap_pax);

#Data splitting

#Rap Src FAK image splitting
rap_src_pre = list();
rap_src_post = list();

assembly_pre_filt = results_rap_onlysignif$assembly$exp_num %% 2 == 0;
assembly_post_filt = results_rap_onlysignif$assembly$exp_num %% 2 == 1;
rap_src_pre$assembly = results_rap_onlysignif$assembly[assembly_pre_filt,];
rap_src_post$assembly = results_rap_onlysignif$assembly[assembly_post_filt,];

disassembly_pre_filt = results_rap_onlysignif$disassembly$exp_num %% 2 == 0;
disassembly_post_filt = results_rap_onlysignif$disassembly$exp_num %% 2 == 1;
rap_src_pre$disassembly = results_rap_onlysignif$disassembly[disassembly_pre_filt,];
rap_src_post$disassembly = results_rap_onlysignif$disassembly[disassembly_post_filt,];

joint_pre_filt = results_rap_onlysignif$joint$exp_num %% 2 == 0;
joint_post_filt = results_rap_onlysignif$joint$exp_num %% 2 == 1;
rap_src_pre$joint = results_rap_onlysignif$joint[joint_pre_filt,];
rap_src_post$joint = results_rap_onlysignif$joint[joint_post_filt,];

rap_src_stage_lengths = gather_stage_lengths(rap_src_pre,rap_src_post);

#Rap Src FAK Control image splitting
rap_src_pre_ctrl = list();
rap_src_post_ctrl = list();

assembly_pre_filt = results_rap_ctrl_onlysignif$assembly$exp_num %% 2 == 0;
assembly_post_filt = results_rap_ctrl_onlysignif$assembly$exp_num %% 2 == 1;
rap_src_pre_ctrl$assembly = results_rap_ctrl_onlysignif$assembly[assembly_pre_filt,];
rap_src_post_ctrl$assembly = results_rap_ctrl_onlysignif$assembly[assembly_post_filt,];

disassembly_pre_filt = results_rap_ctrl_onlysignif$disassembly$exp_num %% 2 == 0;
disassembly_post_filt = results_rap_ctrl_onlysignif$disassembly$exp_num %% 2 == 1;
rap_src_pre_ctrl$disassembly = results_rap_ctrl_onlysignif$disassembly[disassembly_pre_filt,];
rap_src_post_ctrl$disassembly = results_rap_ctrl_onlysignif$disassembly[disassembly_post_filt,];

joint_pre_filt = results_rap_ctrl_onlysignif$joint$exp_num %% 2 == 0;
joint_post_filt = results_rap_ctrl_onlysignif$joint$exp_num %% 2 == 1;
rap_src_pre_ctrl$joint = results_rap_ctrl_onlysignif$joint[joint_pre_filt,];
rap_src_post_ctrl$joint = results_rap_ctrl_onlysignif$joint[joint_post_filt,];

rap_src_ctrl_stage_lengths = gather_stage_lengths(rap_src_pre_ctrl,rap_src_post_ctrl);

#Paxillin Imaging splitting
rap_pax_pre = list();
rap_pax_post = list();

assembly_pre_filt = results_rap_pax_onlysignif$assembly$exp_num %% 2 == 0;
assembly_post_filt = results_rap_pax_onlysignif$assembly$exp_num %% 2 == 1;
rap_pax_pre$assembly = results_rap_pax_onlysignif$assembly[assembly_pre_filt,];
rap_pax_post$assembly = results_rap_pax_onlysignif$assembly[assembly_post_filt,];

disassembly_pre_filt = results_rap_pax_onlysignif$disassembly$exp_num %% 2 == 0;
disassembly_post_filt = results_rap_pax_onlysignif$disassembly$exp_num %% 2 == 1;
rap_pax_pre$disassembly = results_rap_pax_onlysignif$disassembly[disassembly_pre_filt,];
rap_pax_post$disassembly = results_rap_pax_onlysignif$disassembly[disassembly_post_filt,];

joint_pre_filt = results_rap_pax_onlysignif$joint$exp_num %% 2 == 0;
joint_post_filt = results_rap_pax_onlysignif$joint$exp_num %% 2 == 1;
rap_pax_pre$joint = results_rap_pax_onlysignif$joint[joint_pre_filt,];
rap_pax_post$joint = results_rap_pax_onlysignif$joint[joint_post_filt,];

rap_pax_stage_lengths = gather_stage_lengths(rap_pax_pre,rap_pax_post);

pax_lifetimes = determine_mean_p_value(rap_pax_pre$joint$total_lifetime, rap_pax_post$joint$total_lifetime)
fak_lifetimes = determine_mean_p_value(rap_src_pre$joint$total_lifetime, rap_src_post$joint$total_lifetime)

################################################################################
#Plotting
################################################################################
out_folder = '../../doc/publication/figures'
dir.create(out_folder,recursive=TRUE, showWarnings=FALSE);

########################################
#RAP-SRC Plotting
########################################

#Lifetime Plots
dir.create(dirname(file.path(out_folder,'rapr_src','rapr_src_rates_nofilt.pdf')), 
    recursive=TRUE, showWarnings=FALSE);
svg(file.path(out_folder,'rapr_src','rapr_FAK_lifetimes.svg'))
par(bty='n', mar=c(3,4,0,0))
temp = plot_stage_length_data(rap_src_stage_lengths, type='side_by_side', names=c('Before', 'After'));
graphics.off()

svg(file.path(out_folder,'rapr_src','rapr_pax_lifetimes.svg'))
par(bty='n', mar=c(3,4,0,0))
temp = plot_stage_length_data(rap_pax_stage_lengths, type='side_by_side', names=c('Before', 'After'));
graphics.off()

svg(file.path(out_folder,'rapr_src','rapr_overall_lifetimes.svg'))
plot_max = max(pax_lifetimes$conf_ints, fak_lifetimes$conf_ints) + 4;
par(bty='n', mar=c(2,4,0,0))
barplot(cbind(fak_lifetimes$median_vals, pax_lifetimes$median_vals), 
        beside=TRUE, names=c('FAK', 'Paxillin'), legend=c('Before', 'After'), 
        xlim=c(0,7.5), ylim=c(0,plot_max), ylab='Lifetime (min)')
errbar(c(1.5,2.5,4.5,5.5),c(fak_lifetimes$median_vals, pax_lifetimes$median_vals), 
       c(fak_lifetimes$conf_ints[,2],pax_lifetimes$conf_ints[,2]),
       c(fak_lifetimes$conf_ints[,1],pax_lifetimes$conf_ints[,1]), add=TRUE, cex=1E-20, xlab='', ylab='')

pax_max = max(pax_lifetimes$conf_ints) + 0.75;
upper_left = c(4.5, pax_max+1.5);
lower_right = c(5.5, pax_max);
plot_signif_bracket(upper_left, lower_right, over_text='*')
graphics.off()

#Unfiltered Rate Plotting
pdf(file.path(out_folder,'rapr_src','rapr_FAK_rates_nofilt.pdf'))
layout(rbind(c(1,2),c(3,4)))
par(bty='n', mar=c(2,4,1,0))

boxplot_with_points(list(rap_src_pre$assembly$slope, 
                         rap_src_post$assembly$slope),
                    names=c('Before', 'After'), main='Assembly', notch=T,
                    ylab=expression(paste('Rate (',min^-1,')',sep='')))

boxplot_with_points(list(rap_src_pre$disassembly$slope, 
                         rap_src_post$disassembly$slope),
                    names=c('Before', 'After'), main='Disassembly', notch=T,
                    ylab=expression(paste('Rate (',min^-1,')',sep='')))

par(bty='n', mar=c(2,4,1.5,0))

boxplot_with_points(list(rap_src_pre_ctrl$assembly$slope, 
                         rap_src_post_ctrl$assembly$slope),
                    names=c('Before', 'After'), main='Assembly Control', notch=T,
                    ylab=expression(paste('Rate (',min^-1,')',sep='')))

boxplot_with_points(list(rap_src_pre_ctrl$disassembly$slope, 
                         rap_src_post_ctrl$disassembly$slope),
                    names=c('Before', 'After'), main='Disassembly Control', notch=T,
                    ylab=expression(paste('Rate (',min^-1,')',sep='')))
graphics.off()

pdf(file.path(out_folder,'rapr_src','rapr_pax_rates_nofilt.pdf'), height=(7/2))
layout(rbind(c(1,2)))
par(bty='n', mar=c(2,4,1,0))

boxplot_with_points(list(rap_pax_pre$assembly$slope, 
                         rap_pax_post$assembly$slope),
                    names=c('Before', 'After'), main='Assembly', notch=T,
                    ylab=expression(paste('Rate (',min^-1,')',sep='')), point_cex=0.05)

boxplot_with_points(list(rap_pax_pre$disassembly$slope, 
                         rap_pax_post$disassembly$slope),
                    names=c('Before', 'After'), main='Disassembly', notch=T,
                    ylab=expression(paste('Rate (',min^-1,')',sep='')), point_cex=0.05)
graphics.off()


