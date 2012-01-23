rm(list = ls())
source('../../../invado_count_src/analyze_cell_features/invado_analysis_lib.R')
debug=TRUE;

################################################################################
# Data Reading/Initial Processing
################################################################################

degrade_files = list();

base_dir = '../../../results/Invado_count/'
degrade_files$control = Sys.glob(file.path(base_dir,'*/*ntrol/overall*/degrade*'))
degrade_files$BB94 = Sys.glob(file.path(base_dir,'*/BB94/overall*/degrade*'))
degrade_files$DMSO = Sys.glob(file.path(base_dir,'*/DMSO/overall*/degrade*'))
degrade_files$FAK = Sys.glob(file.path(base_dir,'*/FAK/overall*/degrade*'))
degrade_files$PP2 = Sys.glob(file.path(base_dir,'*/PP2/overall*/degrade*'))
degrade_files$Noc = Sys.glob(file.path(base_dir,'*/Noc/overall*/degrade*'))
degrade_files$PurvA = Sys.glob(file.path(base_dir,'*/PurvA/overall*/degrade*'))

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

colors = list(control = rgb(0,0,0), control_light = rgb(0,0,0,0.5),
    DMSO = rgb(1,0,0), DMSO_light = rgb(1,0,0,0.5),
    BB94 = rgb(47/255,143/255,51/255), BB94_light = rgb(24/255,82/255,27/255,0.5),
    PP2 = rgb(77/255,0/255,104/255), PP2_light = rgb(77/255,0/255,104/255,0.5),
    Noc = rgb(255/255,128/255,0/255), Noc_light = rgb(255/255,128/255,0/255,0.5),
    PurvA = rgb(0/255,255/255,241/255), PurvA_light = rgb(0/255,255/255,241/255,0.5),
    FAK = rgb(14/255,19/255,149/255), FAK_light = rgb(14/255,19/255,149/255,0.5)
    )

################################################################################
# Plotting
################################################################################

svg(file.path('single_exp','control.svg'),width=3.5,height=3.5);
par(bty='n',mar=c(2.7,2.5,1.3,0), mgp=c(1.6,0.5,0),xpd=T)
plot_single_exp_mat(time,degrade_percents$control,degrade_conf$control$upper,
    degrade_conf$control$lower, type='l',ylim=c(0,100),
    col=colors$control,xlab='Time Since Plating (hrs)',
    ylab='Percentage of Degradative Cells',xlim=c(0,max(time)))
graphics.off()

svg(file.path('single_exp','DMSO.svg'),width=3.5,height=3.5);
par(bty='n',mar=c(2.7,2.5,1.3,0), mgp=c(1.6,0.5,0),xpd=T)
plot_single_exp_mat(time,degrade_percents$DMSO,degrade_conf$DMSO$upper,
    degrade_conf$DMSO$lower, type='l',ylim=c(0,100),
    col=colors$DMSO,xlab='Time Since Plating (hrs)',
    ylab='Percentage of Degradative Cells',xlim=c(0,max(time)))
graphics.off()

svg(file.path('single_exp','BB94.svg'),width=3.5,height=3.5);
par(bty='n',mar=c(2.7,2.5,1.3,0), mgp=c(1.6,0.5,0),xpd=T)
plot_single_exp_mat(time,degrade_percents$BB94,degrade_conf$BB94$upper,
    degrade_conf$BB94$lower, type='l',ylim=c(0,100),
    col=colors$BB94,xlab='Time Since Plating (hrs)',
    ylab='Percentage of Degradative Cells',xlim=c(0,max(time)))
graphics.off()

svg(file.path('single_exp','FAK.svg'),width=3.5,height=3.5);
par(bty='n',mar=c(2.7,2.5,1.3,0), mgp=c(1.6,0.5,0),xpd=T)
plot_single_exp_mat(time,degrade_percents$FAK,degrade_conf$FAK$upper,
    degrade_conf$FAK$lower,type='l',ylim=c(0,100),
    col=colors$FAK,xlab='Time Since Plating (hrs)',
    ylab='Percentage of Degradative Cells',xlim=c(0,max(time)))
graphics.off()

svg(file.path('single_exp','PP2.svg'),width=3.5,height=3.5);
par(bty='n',mar=c(2.7,2.5,1.3,0), mgp=c(1.6,0.5,0),xpd=T)
plot_single_exp_mat(time,degrade_percents$PP2,degrade_conf$PP2$upper,
    degrade_conf$PP2$lower,type='l',ylim=c(0,100),
    col=colors$PP2,xlab='Time Since Plating (hrs)',
    ylab='Percentage of Degradative Cells',xlim=c(0,max(time)))
graphics.off()

svg(file.path('single_exp','noc.svg'),width=3.5,height=3.5);
par(bty='n',mar=c(2.7,2.5,1.3,0), mgp=c(1.6,0.5,0),xpd=t)
plot_single_exp_mat(time,degrade_percents$noc,degrade_conf$noc$upper,
    degrade_conf$noc$lower,type='l',ylim=c(0,100),
    col=colors$noc,xlab='time since plating (hrs)',
    ylab='percentage of degradative cells',xlim=c(0,max(time)))
graphics.off()

svg(file.path('single_exp','PurvA.svg'),width=3.5,height=3.5);
par(bty='n',mar=c(2.7,2.5,1.3,0), mgp=c(1.6,0.5,0),xpd=t)
plot_single_exp_mat(time,degrade_percents$PurvA,degrade_conf$PurvA$upper,
    degrade_conf$PurvA$lower,type='l',ylim=c(0,100),
    col=colors$PurvA,xlab='time since plating (hrs)',
    ylab='percentage of degradative cells',xlim=c(0,max(time)))
graphics.off()
