this_exp_dir = '../../results/Invadopodia/control/11_22_pos01/adhesion_props/';

degrade_data = read.table(file.path(this_exp_dir,'lin_time_series', 'Local_gel_diff_corr.csv'), 
    sep=",",header=F);
pre_diff_data = read.table(file.path(this_exp_dir,'lin_time_series', 'Pre_birth_diff_corr.csv'), 
    sep=",",header=F);

all_props = list();

i = 62442

only_data = na.omit(as.numeric(degrade_data[i,]));

all_props$mean_vals = c(all_props$mean_vals, mean(only_data))
all_props$max_local_diff = c(all_props$max_local_diff, max(only_data));
all_props$min_local_diff = c(all_props$min_local_diff, min(only_data));

test_results = tryCatch(t.test(only_data,conf.level=0.95), error = t.test.error);

all_props$high_conf_int = c(all_props$high_conf_int, test_results$conf.int[2]);
all_props$p_value = c(all_props$p_value, test_results$p.value);

#Pre-birth local difference data
only_pre_diff_data = na.omit(as.numeric(pre_diff_data[i,]));

only_pre_diff_test = tryCatch(t.test(only_pre_diff_data,conf.level=0.95), error = t.test.error);

pre_test = tryCatch(t.test(only_data - only_pre_diff_data,conf.level=0.95), error = t.test.error);
all_props$pre_p_value = c(all_props$pre_p_value, pre_test$p.value);
all_props$pre_high_conf_int = c(all_props$pre_high_conf_int, pre_test$conf.int[2]);

#Plotting the data
time_points = seq(from=0,by=5,along.with=only_data);

svg(file.path('UCRF_grant','sm_data_plot.svg'),height=8)

par(bty='n', mar=c(4,4,1,0))
layout(rbind(c(1),c(2)), heights=c(0.2,1))

#place holder for small multiple diagram
plot.new()
mtext('A',adj=-0.13,side=3,line=-0.5,cex=2)

#plot all the lines and error bars in the second plot
matplot(time_points, cbind(only_data, only_pre_diff_data, only_data - only_pre_diff_data), typ='l', xlab='Time (min)', ylab='Difference Metric', lwd=3)

plot_limits = par("usr");
legend(plot_limits[1],plot_limits[4]*0.95,c('Local Diff','Pre-birth Local Diff', 'Local Diff - Pre-birth Diff' ), fill=c('black','red', 'green'))
segments(0,0,max(time_points),0, col='purple', lty=4)

plot_limits = par("usr");

errbar(max(time_points)+(plot_limits[2]-max(time_points))*0.4, test_results$estimate, test_results$conf.int[2], test_results$conf.int[1], add=T)
errbar(max(time_points)+(plot_limits[2]-max(time_points))*0.4, only_pre_diff_test$estimate, 
    only_pre_diff_test$conf.int[2], only_pre_diff_test$conf.int[1], add=T, col='red')
errbar(max(time_points)+(plot_limits[2]-max(time_points))*0.4, pre_test$estimate, 
    pre_test$conf.int[2], pre_test$conf.int[1], add=T, col='green')
mtext('B',adj=-0.13,side=3,line=-2,cex=2)

graphics.off()
