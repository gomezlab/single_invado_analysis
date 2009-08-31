rm(list = ls())
source('FA_analysis_lib.R')
library(lattice)
library(geneplotter)
library(Hmisc)

################################################################################
#Result loading
################################################################################
raw_data <- list()
single_props <- list()

#Rapamycin results
exp_dirs_rap <- Sys.glob('../../results/rap_src/*/adhesion_props/models/')
exp_dirs_rap <- exp_dirs_rap[file_test('-d',exp_dirs_rap)]
raw_data$results_rap = load_results(exp_dirs_rap,file.path('intensity.Rdata'))
raw_data$corr_results_rap = load_results(exp_dirs_rap,file.path('local_corrected.Rdata'))

exp_dirs_rap <- Sys.glob('../../results/rap_src_control/*/adhesion_props/models/')
exp_dirs_rap <- exp_dirs_rap[file_test('-d',exp_dirs_rap)]
raw_data$results_rap_ctrl = load_results(exp_dirs_rap,file.path('intensity.Rdata'))
raw_data$corr_results_rap_ctrl = load_results(exp_dirs_rap,file.path('local_corrected.Rdata'))

print('Done Loading Data')

################################################################################
#Plotting
################################################################################
out_folder = '../../doc/publication/figures'
dir.create(out_folder,recursive=TRUE, showWarnings=FALSE);

########################################
#RAP-SRC Plotting
########################################

#Data splitting
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

#Lifetime Plot
dir.create(dirname(file.path(out_folder,'rapr_src','rapr_src_rates_nofilt.pdf')), 
    recursive=TRUE, showWarnings=FALSE);
svg(file.path(out_folder,'rapr_src','rapr_src_lifetimes.svg'))
par(bty='n', mar=c(3,4,0,0))
temp = plot_stage_length_data(rap_src_stage_lengths, type='side_by_side', names=c('Before', 'After'));
graphics.off()

#Unfiltered Rate Plotting
pdf(file.path(out_folder,'rapr_src','rapr_src_rates_nofilt.pdf'))
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
