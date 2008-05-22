plot_ad_seq <- function (results,index,dir) {
	
	ad_seq = results$exp_data[index,]
	ad_seq = t(ad_seq[!(is.nan(ad_seq))]);
	if (! file.exists(dir)) {
		dir.create(dir)
	}
	
	pdf(paste(dir,index,'.pdf',sep=''));
	par(mfrow=c(1,2),bty='n', mar=c(5,4,1,1))
			
	if (is.finite(results$early_slope[index])) {
		early_ad_seq = ad_seq[1:results$early_offset[index]];
#		print(this_ad_seq)
		early_ad_seq = early_ad_seq/early_ad_seq[1];
#		print(log(this_ad_seq))
		
		plot(1:results$early_offset[index],early_ad_seq,xlab='Time (minutes)',ylab='Normalized Paxillin Signal')
		
		x = c(0,results$early_offset[index]);
		y = c(results$early_slope[index]*x[1] + results$early_inter[index],
		   	  results$early_slope[index]*x[2] + results$early_inter[index])
		
		lines(x,y,col='red',lwd=2)
	}
	
	if (is.finite(results$late_slope[index])) {
		this_ad_seq = ad_seq[(length(ad_seq) - results$late_offset[index]) : length(ad_seq)];
		this_ad_seq = this_ad_seq/this_ad_seq[1];
		
		plot((length(ad_seq) - results$late_offset[index]) : length(ad_seq),
			 this_ad_seq, xlab='Time (minutes)', ylab='Normalized Paxillin Signal')
		
		x = c(length(ad_seq) - results$late_offset[index],length(ad_seq));
		y = c(results$late_slope[index]*x[1] + results$late_inter[index],
		   	  results$late_slope[index]*x[2] + results$late_inter[index])
		
		lines(x,y,col='red',lwd=2)
	}
		
	dev.off();
}

#plot_ad_seq(norm_all_data[[1]],498,'./')

gather_linear_regions <- function (dirs, min_length=10, data_file='Average_adhesion_signal.csv', col_lims = NaN, normed = 1, log_lin = 0) {
	
	results = list()
	
	for (k in 1:length(dirs)) {
		
		if (is.na(dirs[[k]])) {
			next
		}
		
		print(dirs[[k]])
		
		ad_sig <- read.table(paste(dirs[[k]],data_file,sep=''),header = FALSE, sep  = ',');
		ad_props <- read.table(paste(dirs[[k]],'../single_lin.csv',sep=''), header = TRUE, sep=',');
		
		if (is.numeric(col_lims) && length(col_lims) == 2) {
			ad_sig = ad_sig[,col_lims[1]:col_lims[2]];
		}
		
		rows <- dim(ad_sig)[[1]]
		cols <- dim(ad_sig)[[2]]

		m_results <- list(early_offset = array(min_length, dim = c(rows)),
						  early_R_sq = array(0, dim = c(rows)),
						  early_inter = array(NaN, dim = c(rows)), 
						  early_slope = array(NaN, dim = c(rows)),
						  late_offset = array(min_length, dim = c(rows)),
						  late_R_sq = array(0, dim = c(rows)),
						  late_inter = array(NaN, dim = c(rows)), 
						  late_slope = array(NaN, dim = c(rows)),
						  exp_data = ad_sig,
						  exp_props = ad_props
						 );
	
		for (i in 1:rows) {
			temp = ad_sig[i,];
			temp = temp[!(is.nan(temp))];
			
			if (dim(temp)[[2]] == 0) {
				m_results$seq_length[[i]] = 0;
				next;
			}	
			
			this_ad_sig = data.frame(y = t(temp[1,]), x = 1:dim(temp)[[2]]);
			names(this_ad_sig) <- c('y','x');
		
			m_results$seq_length[[i]] = length(temp);
		
			#Skip over adhesions which don't live long enough
			if (m_results$seq_length[[i]] < min_length * 2) {
				next
			}
			#Skip over adhesions where we don't see the entire life cycle
			if (is.finite(ad_sig[i,1])) {
				next
			}
	
			cur_offset = m_results$early_offset[[i]];
			for (j in cur_offset:dim(this_ad_sig)[[1]]) {
				early_subset = this_ad_sig[1:cur_offset,]
				if (normed) {
					early_subset$y = early_subset$y/early_subset$y[1]
				}
				if (log_lin) {
					early_subset$y = log(early_subset$y)
				}
				
				model <- lm(y ~ x, data = early_subset)
		
				summary <- summary(model);
			
				if (summary$r.squared > m_results$early_R_sq[[i]]) {
					m_results$early_R_sq[[i]] = summary$r.squared;
					m_results$early_length[[i]] = dim(early_subset)[[1]];
					m_results$early_model[[i]] = model;
					m_results$early_offset[[i]] = cur_offset;
					m_results$early_inter[[i]] = coef(model)[[1]];
					m_results$early_slope[[i]] = coef(model)[[2]];
				}
				cur_offset = cur_offset + 1
			}
			
			cur_offset = m_results$late_offset[[i]]
			for (j in cur_offset:dim(this_ad_sig)[[1]]) {
				late_subset = this_ad_sig[(dim(this_ad_sig)[[1]]-cur_offset):dim(this_ad_sig)[[1]],]
				if (normed) {
#					late_subset$y = late_subset$y/late_subset$y[[length(late_subset$y)]]
					late_subset$y = late_subset$y[[1]]/late_subset$y
				}
				if (log_lin) {
					late_subset$y = log(late_subset$y)
				}
	
				model <- lm(y ~ x, data = late_subset)
				
				summary <- summary(model);
			
				if (summary$r.squared > m_results$late_R_sq[[i]]) {
					m_results$late_R_sq[[i]] = summary$r.squared;
					m_results$late_length[[i]] = dim(late_subset)[[1]];
					m_results$late_model[[i]] = model;
					m_results$late_offset[[i]] = cur_offset;
					m_results$late_inter[[i]] = coef(model)[[1]];
					m_results$late_slope[[i]] = coef(model)[[2]];
				}
				cur_offset = cur_offset + 1
			}
		}
		results[[k]] <- m_results
	}
	results
}

plot_lin_reg_set <- function(results,dir,file='linear_regions.pdf', hist_file=NA) {
	
	early_slope = c()
	early_error = c()
	early_n = c()
	late_slope  = c()
	late_error  = c()
	late_n = c()
	for (i in 1:10) {
		slope_ests = exp_set_slope_estimate(results,r_cutoff = ((i - 1)/10));
		early_slope[i] = slope_ests$early
		early_n[i] = slope_ests$early_n
		early_error[i] = slope_ests$early_sem
		late_slope[i] = slope_ests$late
		late_n[i] = slope_ests$late_n
		late_error[i] = slope_ests$late_sem
	}
	
	library(Hmisc)
    pdf(paste(dir,file,sep=''),width=8.5,height=8.5,pointsize=14);
    par(mfrow=c(2,2),bty='n', mar=c(5,4,1,1))
	
	#Plot 1 - Slope versus R squared (early)
#    plot(results$early_slope[results$early_R_sq != 0],
#         results$early_R_sq[results$early_R_sq != 0],
#         xlab='Formation Slope', ylab='R Squared'
#        );

	slopes = results$early_slope 
	R_sqs = results$early_R_sq
	init_filter = R_sqs != 0;
	all_x = slopes[init_filter]
	all_y = R_sqs[init_filter]
	cols = rainbow(10);
	c_l = length(cols);

	dist_list = results$exp_props$starting_edge_dist
	dist_range = c(min(dist_list[init_filter]),max(dist_list[init_filter]));
	
	indexes = init_filter & dist_list <= (1/c_l)*(dist_range[2]-dist_range[1])
	x = slopes[indexes]
	y = R_sqs[indexes]
	plot(x, y, col=cols[1], 
		 xlim=c(min(all_x),max(all_x)), ylim=c(min(all_y),max(all_y)), 
		 cex = 0.5, xlab='Formation Slope (/min)', ylab='R Squared')
	
	for (i in 2:c_l) {
		indexes = init_filter & dist_list <= (i/c_l)*(dist_range[2]-dist_range[1])
		indexes = indexes & dist_list > ((i-1)/c_l)*(dist_range[2]-dist_range[1])
		x = slopes[indexes]
		y = R_sqs[indexes]
		points(x,y,col=cols[i], cex = 0.5)
	}
	
	x_range = c(0.14,0.15);
	y_range = c(0.1,0.6);
	
	for (i in 1:c_l) {
		rect(x_range[1], ((i-1)/c_l)*(y_range[2] - y_range[1]) + y_range[1],
			 x_range[2], (i/c_l)*(y_range[2] - y_range[1]) + y_range[1],
			 col=cols[i], border = NA)
	}
	text(x_range[2],y_range[1],sprintf('%.1f',dist_range[1]),pos=4)
	text(x_range[2],y_range[2],sprintf('%.1f',dist_range[2]),pos=4)

	#Plot 2 - R squared versus Slope (late)
#    plot(results$late_slope[results$late_R_sq != 0 & results$exp_props$death_status],
#		  results$late_R_sq[results$late_R_sq != 0 & results$exp_props$death_status],
#         xlab='Decay Slope (/min)', ylab='R Squared', cex = 0.5
#        )

	slopes = results$late_slope
	R_sqs = results$late_R_sq
	init_filter = R_sqs != 0  & results$exp_props$death_status
	all_x = slopes[init_filter]
	all_y = R_sqs[init_filter]
	cols = rainbow(10);
	c_l = length(cols);

	dist_list = results$exp_props$starting_edge_dist
	dist_range = c(min(dist_list[init_filter]),max(dist_list[init_filter]));	
	indexes = init_filter & dist_list <= (1/c_l)*(dist_range[2]-dist_range[1])
	x = slopes[indexes]
	y = R_sqs[indexes]
	plot(x, y, col=cols[1], 
		 xlim=c(min(all_x),max(all_x)), ylim=c(min(all_y),max(all_y)), 
		 cex = 0.5, xlab='Decay Slope (/min)', ylab='R Squared')
	
	for (i in 2:c_l) {
		indexes = init_filter & dist_list <= (i/c_l)*(dist_range[2]-dist_range[1])
		indexes = indexes & dist_list > ((i-1)/c_l)*(dist_range[2]-dist_range[1])
		x = slopes[indexes]
		y = R_sqs[indexes]
		points(x,y,col=cols[i], cex = 0.5)
	}
	
	x_range = c(max(all_x)*0.8,max(all_x)*0.85);
	y_range = c(0.1,0.6);
	
	for (i in 1:c_l) {
		rect(x_range[1], ((i-1)/c_l)*(y_range[2] - y_range[1]) + y_range[1],
			 x_range[2], (i/c_l)*(y_range[2] - y_range[1]) + y_range[1],
			 col=cols[i], border = NA)
	}
	text(x_range[2],y_range[1],sprintf('%.1f',dist_range[1]),pos=4)
	text(x_range[2],y_range[2],sprintf('%.1f',dist_range[2]),pos=4)
    
	#Plot 3
    errbar(seq(0,0.9,by=0.1), early_slope, early_slope - early_error, early_slope + early_error, xlab = 'R squared cutoff', ylab='Accumulation Rate (/min)')
    for (i in 1:10) {
    	x_points = seq(0,0.9,by=0.1)
    	text(x_points[i],early_slope[i]+early_error[i],early_n[i])
    }
	
	#Plot 4
    errbar(seq(0,0.9,by=0.1), late_slope, late_slope - late_error, late_slope + late_error, xlab = 'R squared cutoff', ylab='Decay Rate (/min)')
    for (i in 1:10) {
    	x_points = seq(0,0.9,by=0.1)
    	text(x_points[i],late_slope[i]+late_error[i],late_n[i])
    }
	
	if (! is.na(hist_file)) {
	    pdf(paste(dir,hist_file,sep=''),width=8.5,height=8.5,pointsize=14);
    	par(mfrow=c(2,2),bty='n', mar=c(5,4,1,1))
        
	    #Plot 1
    	hist(results$early_length, xlab = 'Linear Sequence Length (min)', main='')
    
	    #Plot 2
    	hist(results$late_length, xlab = 'Linear Sequence Length (min)', main='')
    
	    #Plot 3
    	hist(results$early_slope[results$early_R_sq > 0.9], xlab = 'Slope (/min)', main='')

	    #Plot 4
    	hist(results$late_slope[results$late_R_sq > 0.9 & results$exp_props$death_status], xlab = 'Slope (/min)', main='')
    	dev.off()
    }
    
    dev.off()
}

exp_set_slope_estimate <- function(results,r_cutoff=0.9) {
	early_slopes <- c()
	late_slopes <- c()
	
	if (! is.null(names(results))) {
		results = list(results);
	}	
	for (i in 1:length(results)) {
		res = results[[i]];
#		early_slopes <- c(early_slopes, res$early_slope[res$early_R_sq > r_cutoff & res$early_slope > 0]);
#		late_slopes  <- c(late_slopes, res$late_slope[res$late_R_sq > r_cutoff & res$late_slope < 0]);
		early_slopes <- c(early_slopes, res$early_slope[res$early_R_sq > r_cutoff]);
		late_slopes  <- c(late_slopes, res$late_slope[res$late_R_sq > r_cutoff & res$exp_props$death_status]);
	}
#	print(length(early_slopes))
#	print(length(late_slopes))
	slopes <- list(early = mean(early_slopes),
				   early_n = length(early_slopes),
				   early_cv = sd(early_slopes)/mean(early_slopes),
				   early_sem = sd(early_slopes)/sqrt(length(early_slopes)),
				   early_sd = sd(early_slopes),
				   late = mean(late_slopes),	
				   late_n = length(late_slopes),
				   late_cv = sd(late_slopes)/mean(late_slopes),
				   late_sem = sd(late_slopes)/sqrt(length(late_slopes)),
				   late_sd = sd(late_slopes)
				  )
}

####################################################################################################
# Main Program
####################################################################################################

#Paxillin
prefix = '../../results/focal_adhesions/';

dirs = c('time_series_01/adhesion_props/lin_time_series/','time_series_04/adhesion_props/lin_time_series/',
         'time_series_05/adhesion_props/lin_time_series/','time_series_06/adhesion_props/lin_time_series/',
         'time_series_07/adhesion_props/lin_time_series/','time_series_08/adhesion_props/lin_time_series/',
         'time_series_09/adhesion_props/lin_time_series/','time_series_10/adhesion_props/lin_time_series/',
         'time_series_11/adhesion_props/lin_time_series/','time_series_12/adhesion_props/lin_time_series/',
         'time_series_13/adhesion_props/lin_time_series/','time_series_14/adhesion_props/lin_time_series/',
         'time_series_15/adhesion_props/lin_time_series/','time_series_16/adhesion_props/lin_time_series/',
         'time_series_17/adhesion_props/lin_time_series/','time_series_18/adhesion_props/lin_time_series/',
         'time_series_19/adhesion_props/lin_time_series/','time_series_20/adhesion_props/lin_time_series/',
         'time_series_21/adhesion_props/lin_time_series/','time_series_22/adhesion_props/lin_time_series/',
         'time_series_23/adhesion_props/lin_time_series/');

pre_dirs <- c()
for (i in 1:length(dirs)) {
#for (i in 1:1) {
	pre_dirs[i] = paste(prefix,dirs[i],sep='')
}

#print('No Norm')
#all_data = gather_linear_regions(pre_dirs,normed = 0)
print('Norm')
norm = gather_linear_regions(pre_dirs)
print('Norm Log lin')
log = gather_linear_regions(pre_dirs,log_lin = 1)

#for (i in 1:length(dirs)) {
#for (i in 1:1) {
#	if (is.na(pre_dirs[i])) {
#		next
#	}	
#	plot_lin_reg_set(norm[[i]],paste(pre_dirs[i],'../plots/',sep=''),hist_file='lin_hist.pdf')
#}

#all_est <- exp_set_slope_estimate(all_data[[1]])
#norm_est <- exp_set_slope_estimate(norm_all_data[[1]])
#log_est <- exp_set_slope_estimate(log_all_data[[1]])
#
#
##FAK
#prefix = '../../results/FAK/';
#
#dirs = c('time_series_01/adhesion_props/lin_time_series/');
#
#pre_dirs <- c()
#for (i in 1:length(dirs)) {
##for (i in 2:2) {	
#	pre_dirs[i] = paste(prefix,dirs[i],sep='')
#}

#pre_fak = gather_linear_regions(pre_dirs,col_lims = c(1,33))
#post_fak = gather_linear_regions(pre_dirs,col_lims = c(34,90))

#for (i in 1:length(dirs)) {
#	if (is.na(pre_dirs[i])) {
#		next
#	}	
#	plot_lin_reg_set(norm_all_data[[i]],paste(pre_dirs[i],'../plots/',sep=''))
#}

#rm(i)