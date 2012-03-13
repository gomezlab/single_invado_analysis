rm(list = ls())
source('../../../invado_count_src/analyze_cell_features/invado_analysis_lib.R')
source('misc_plotting_properties.R')
debug=TRUE;

################################################################################
# Data Reading/Initial Processing
################################################################################

degrade_files = list();

base_dir = '../../../results/Invado_count_old/'
degrade_files$control = Sys.glob(file.path(base_dir,'*/*ntrol/overall*/degrade_percentage_70.csv'))
degrade_files$BB94 = Sys.glob(file.path(base_dir,'*/BB94/overall*/degrade_percentage_70.csv'))
degrade_files$DMSO = Sys.glob(file.path(base_dir,'*/DMSO/overall*/degrade_percentage_70.csv'))
degrade_files$FAK = Sys.glob(file.path(base_dir,'*/FAK/overall*/degrade_percentage_70.csv'))
degrade_files$PP2 = Sys.glob(file.path(base_dir,'*/PP2/overall*/degrade_percentage_70.csv'))
degrade_files$Noc = Sys.glob(file.path(base_dir,'*/Noc/overall*/degrade_percentage_70.csv'))
degrade_files$PurvA = Sys.glob(file.path(base_dir,'*/PurvA/overall*/degrade_percentage_70.csv'))

degrade_percents = list()

for (exp_type in names(degrade_files)) {
    for (i in 1:length(degrade_files[[exp_type]])) {
        temp = as.numeric(read.csv(degrade_files[[exp_type]][i],header=F));
        degrade_percents[[exp_type]] = rbind(degrade_percents[[exp_type]],temp*100);
    }
}

degrade_conf = list()
for (exp_type in names(degrade_files)) {
    degrade_conf[[exp_type]]$upper = colConfUpper(degrade_percents[[exp_type]])
    degrade_conf[[exp_type]]$lower = colConfLower(degrade_percents[[exp_type]])
}

time = seq(1.25,along.with=degrade_percents$control[1,],by=0.5)

colors = get_exp_colors()

################################################################################
# Plotting
################################################################################

###########################################################
# Percentage of cells that have degraded
###########################################################

########################################
# Single Exp
########################################

svg(file.path('single_exp','control.svg'),width=3.5,height=3.5);
par(bty='n',mar=c(2.7,2.5,0.3,0), mgp=c(1.6,0.5,0),xpd=T)
source('../../../invado_count_src/analyze_cell_features/invado_analysis_lib.R')
plot_single_exp_mat(time,degrade_percents$control,degrade_conf$control$upper,
    degrade_conf$control$lower, type='l',ylim=c(0,100),
    col=colors$control,xlab='Time Since Plating (hrs)',
    ylab='Fraction of Cells with Invadopodia',xlim=c(0,max(time)))
graphics.off()

svg(file.path('single_exp','DMSO.svg'),width=3.5,height=3.5);
par(bty='n',mar=c(2.7,2.5,0.3,0), mgp=c(1.6,0.5,0),xpd=T)
plot_single_exp_mat(time,degrade_percents$DMSO,degrade_conf$DMSO$upper,
    degrade_conf$DMSO$lower, type='l',ylim=c(0,100),
    col=colors$DMSO,xlab='Time Since Plating (hrs)',
    ylab='Fraction of Cells with Invadopodia',xlim=c(0,max(time)))
graphics.off()

svg(file.path('single_exp','BB94.svg'),width=3.5,height=3.5);
par(bty='n',mar=c(2.7,2.5,0.3,0), mgp=c(1.6,0.5,0),xpd=T)
plot_single_exp_mat(time,degrade_percents$BB94,degrade_conf$BB94$upper,
    degrade_conf$BB94$lower, type='l',ylim=c(0,100),
    col=colors$BB94,xlab='Time Since Plating (hrs)',
    ylab='Fraction of Cells with Invadopodia',xlim=c(0,max(time)))
graphics.off()

svg(file.path('single_exp','FAK.svg'),width=3.5,height=3.5);
par(bty='n',mar=c(2.7,2.5,0.3,0), mgp=c(1.6,0.5,0),xpd=T)
plot_single_exp_mat(time,degrade_percents$FAK,degrade_conf$FAK$upper,
    degrade_conf$FAK$lower,type='l',ylim=c(0,100),
    col=colors$FAK,xlab='Time Since Plating (hrs)',
    ylab='Fraction of Cells with Invadopodia',xlim=c(0,max(time)))
graphics.off()

svg(file.path('single_exp','PP2.svg'),width=3.5,height=3.5);
par(bty='n',mar=c(2.7,2.5,0.3,0), mgp=c(1.6,0.5,0),xpd=T)
plot_single_exp_mat(time,degrade_percents$PP2,degrade_conf$PP2$upper,
    degrade_conf$PP2$lower,type='l',ylim=c(0,100),
    col=colors$PP2,xlab='Time Since Plating (hrs)',
    ylab='Fraction of Cells with Invadopodia',xlim=c(0,max(time)))
graphics.off()

svg(file.path('single_exp','Noc.svg'),width=3.5,height=3.5);
par(bty='n',mar=c(2.7,2.5,0.3,0), mgp=c(1.6,0.5,0),xpd=t)
plot_single_exp_mat(time,degrade_percents$Noc,degrade_conf$Noc$upper,
    degrade_conf$Noc$lower,type='l',ylim=c(0,100),
    col=colors$Noc,xlab='Time since plating (hrs)',
    ylab='Fraction of Cells with Invadopodia',xlim=c(0,max(time)))
graphics.off()

svg(file.path('single_exp','PurvA.svg'),width=3.5,height=3.5);
par(bty='n',mar=c(2.7,2.5,0.3,0), mgp=c(1.6,0.5,0),xpd=t)
plot_single_exp_mat(time,degrade_percents$PurvA,degrade_conf$PurvA$upper,
    degrade_conf$PurvA$lower,type='l',col=colors$PurvA,
    ylim=c(0,100),ylab='Fraction of Cells with Invadopodia',
    xlim=c(0,max(time)),xlab='Time since Plating (hrs)')
graphics.off()

########################################
# Multiple Exp
########################################

source('../../../invado_count_src/analyze_cell_features/invado_analysis_lib.R')

svg(file.path('degrade_curves.svg'),width=3,height=3);
par(bty='n',mar=c(2.7,2.5,0.3,0), mgp=c(1.6,0.5,0),xpd=t)
plot_multiple_exp_mat(time,degrade_percents,degrade_conf,
    cols=colors,exp_names=c("BB94","DMSO"),
    ylim=c(0,55),ylab='Fraction of Cells with Invadopodia',
    xlim=c(0,max(time)),xlab='Time since Plating (hrs)')
legend('topleft',c('DMSO','BB94'),
    fill=c(colors$DMSO,colors$BB94),bty='n')
graphics.off()
