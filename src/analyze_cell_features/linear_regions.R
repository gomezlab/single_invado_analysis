################################################################################
#linear_regions.R: various functions used to find and plot the linear regions
#  and associated data from the focal adhesion identification/analysis programs
################################################################################

########################################
#Data fitting functions
########################################

gather_bilinear_models_from_dirs <- function (dirs, min_length=10, 
	data_file='Average_adhesion_signal.csv', col_lims = NA, 
	normed = TRUE, log.trans = TRUE, boot.samp = NA, results.file = NA,
	save.exp_data = TRUE, debug = FALSE) {
	
	results = list()
	
	for (i in 1:length(dirs)) {
		if (is.na(dirs[[i]])) {
			next
		}
		
		if (debug) {
			print(dirs[[i]])
		}
			
		exp_data <- as.matrix(read.table(file.path(dirs[[i]],data_file),header = FALSE, sep  = ','));

		exp_props <- read.table(file.path(dirs[[i]],'../single_lin.csv'), header = TRUE, sep=',');
		
		#process the col_lim parameter passed in if values were passed in
		this_col_lim = NA 
		if (! is.na(as.matrix(col_lims)[1,1])) {
			if(dim(as.matrix(col_lims))[[2]] == 1) {
				this_col_lim = c(col_lims[i],dim(exp_data)[[2]])
			} else {
				this_col_lim = col_lims[i,]
			}
		}
		
		results[[i]] <- gather_bilinear_models(exp_data, exp_props, 
							min_length = min_length, col_lims = this_col_lim, 
							normed = normed, log.trans = log.trans, 
							boot.samp = boot.samp, save.exp_data = save.exp_data, debug=debug)
        
        results[[i]]$exp_dir = dirs[[i]]
        if (! is.na(results.file)) {
            this_result = results[[i]]
            save(this_result,file = file.path(dirs[[i]],results.file))
        }
	}
	
	results
}

gather_bilinear_models <- function(data_set, props, 
	min_length = 10, col_lims = NA, normed = TRUE, 
	log.trans = TRUE, boot.samp = NA, save.exp_data = TRUE, debug = FALSE) {
		
	if (is.numeric(col_lims) && length(col_lims) == 2) {
		data_set = data_set[,col_lims[1]:col_lims[2]];
	}

	rows <- dim(data_set)[[1]]
	cols <- dim(data_set)[[2]]
		  
	results <- list()
	results$exp_props = props
	
	if (save.exp_data) {
		results$exp_data = data_set;
	}
	results$stable_data_set = list();
	for (i in 1:rows) {
		if (i > 500) {
			#next
		}
		
		this_data_set = as.vector(data_set[i,])
		numeric_data_set = this_data_set[! is.nan(this_data_set)]
		if (length(numeric_data_set) == 0) {
			next
		}

		if (i %% 100 == 0 & debug) {
			print(i)
		}
		
		#Skip over adhesions which don't live long enough
		if (length(numeric_data_set) < (min_length * 2)) {
			next
		}
		
		these_exp_props = results$exp_props[i,]
		
		temp_results = find_optimum_bilinear_fit(this_data_set, these_exp_props, normed = normed, 
			min_length = min_length, log.trans = log.trans)

		results$early$offset[i] = temp_results$early$offset
		results$early$R_sq[i]   = temp_results$early$R_sq
		results$early$inter[i]  = temp_results$early$inter
		results$early$slope[i]  = temp_results$early$slope
		results$early$fold_change[i]  = temp_results$early$fold_change
		results$early$residual[i]  = temp_results$early$residual
		
		results$late$offset[i]  = temp_results$late$offset
		results$late$R_sq[i]    = temp_results$late$R_sq
		results$late$inter[i]   = temp_results$late$inter
		results$late$slope[i]   = temp_results$late$slope
		results$late$fold_change[i]  = temp_results$late$fold_change
		results$late$residual[i]   = temp_results$late$residual
		
		#if either of the offset values are NA, then we weren't able 
		#to get a fit for one side of the data, don't calculate a 
		#stable lifetime because we don't know how long the adhesion 
		#was around before the movie starts or after it ends
		if (! is.na(results$late$offset[i]) && ! is.na(results$early$offset[i])) {
			results$stable_data_set[[i]] = numeric_data_set[results$early$offset[i]:(length(numeric_data_set) - results$late$offset[i])]
			
			results$stable_lifetime[i] = length(results$stable_data_set[[i]])
			results$stable_mean[i] = mean(results$stable_data_set[[i]])
			results$stable_variance[i] = var(results$stable_data_set[[i]])
			if (length(results$stable_data_set[[i]]) == 1) {
				results$stable_variance[i] = 0;
			}
		}
	}

	results <- pad_results_to_row_length(results, rows)

	if (is.numeric(boot.samp)) {
		results$sim_results <- gather_linear_regions.boot(results, min_length = min_length, 
			col_lims = col_lims, normed = normed, log.trans = log.trans, boot.samp = boot.samp)
	}
	
	results
}

pad_results_to_row_length <- function(results, desired_length) {
	for (i in 1:length(results$early)) {
		for (j in (length(results$early[[i]]) + 1):desired_length) {
			results$early[[i]][j] = NA
		}
	}
	for (i in 1:length(results$late)) {
		for (j in (length(results$late[[i]]) + 1):desired_length) {
			results$late[[i]][j] = NA
		}
	}
	for (i in (length(results$stable_lifetime) + 1):desired_length) {
		results$stable_lifetime[i] = NA
		results$stable_data_set[i] = NA
		results$stable_variance[i] = NA
		results$stable_mean[i] = NA		
	}
	results
}

find_optimum_bilinear_fit <- function(initial_data_set, exp_props, normed = TRUE, min_length = 10, log.trans = TRUE) {

	results = list(initial_data_set = initial_data_set)
	resid = list(early = list(), late = list())	
	results$filt_init = initial_data_set[! is.nan(initial_data_set)]
	if (length(results$filt_init) == 0) {
		return(c(NA, NA))
	}
	this_data_set = data.frame(y = results$filt_init, x = 1:length(results$filt_init))

	#Search the beginning of the sequence for a linear fit
	early_slope_calculated = FALSE;
	if (is.nan(initial_data_set[1]) & ! exp_props$split_birth_status) {
		early_slope_calculated = TRUE;
		for (j in min_length:dim(this_data_set)[[1]]) {
			early_subset = this_data_set[1:j,]
			if (normed) {
				early_subset$y = early_subset$y/early_subset$y[1]
			}
			if (log.trans) {
				early_subset$y = log(early_subset$y)
			}
				
			model <- lm(y ~ x, data = early_subset)
			summary <- summary(model);
			
			results$early$R_sq[j] = summary$adj.r.squared
			#dealing with a degerate case, where lm produce NaN for the R squared 
			#value when the data set is a flat line, see:
			#	>data <- data.frame(x = c(1,2,3), y = c(1,1,1))
			#	>summary(lm(y ~ x, data=data))
			if (is.nan(results$early$R_sq[j])) {
				results$early$R_sq[j] = 1
			}
			
			results$early$length[j] = dim(early_subset)[[1]]
			results$early$offset[j] = j
			results$early$inter[j] = coef(model)[[1]]
			results$early$slope[j] = coef(model)[[2]]

			if (log.trans) {
				results$early$fold_change[j] = max(early_subset$y)
			} else {				
				results$early$fold_change[j] = max(early_subset$y)/min(early_subset$y)
			}
			resid$early[[j]] = as.numeric(resid(model))
		}
	} else {
		results$early$R_sq[1] = 0
		
		results$early$length[1] = NA
		results$early$offset[1] = NA
		results$early$inter[1] = NA
		results$early$slope[1] = NA
		results$early$fold_change[1] = NA

		resid$early[[1]] = NA
	}
	
	#Search the end of the sequence for a linear fit
	late_slope_calculated = FALSE;
	if (is.nan(initial_data_set[length(initial_data_set)]) & exp_props$death_status) {
		late_slope_calculated = TRUE;
		for (j in min_length:dim(this_data_set)[[1]]) {
			late_subset = this_data_set[(dim(this_data_set)[[1]]-j):dim(this_data_set)[[1]],]
			if (normed) {
				late_subset$y = late_subset$y[1]/late_subset$y
			}
			if (log.trans) {
				late_subset$y = log(late_subset$y)
			}
	
			model <- lm(y ~ x, data = late_subset)
			summary <- summary(model);
			
			results$late$R_sq[j] = summary$adj.r.squared
			#dealing with a degerate case, where lm produce NaN for the R squared 
			#value when the data set is a flat line, see:
			#	>data <- data.frame(x = c(1,2,3), y = c(1,1,1))
			#	>summary(lm(y ~ x, data=data))
			if (is.nan(results$late$R_sq[j])) {
				results$late$R_sq[j] = 1
			}
			
			results$late$length[j] = dim(late_subset)[[1]]
			results$late$offset[j] = j
			results$late$inter[j] = coef(model)[[1]]
			results$late$slope[j] = coef(model)[[2]]

			if (log.trans) {
				results$late$fold_change[j] = max(late_subset$y)
			} else {						
				results$late$fold_change[j] = max(late_subset$y)/min(late_subset$y)
			}
			resid$late[[j]] = as.numeric(resid(model))
		}
	} else {
		results$late$R_sq[1] = 0
				
		results$late$length[1] = NA
		results$late$offset[1] = NA
		results$late$inter[1] = NA
		results$late$slope[1] = NA
		results$late$fold_change[1] = NA

		resid$late[[1]] = NA	
	}

	best_indexes = find_best_offset_combination(results, min_length = min_length)
	
	#With the R squared matrix calculated reset the r_sq componenets to NA, if needed since 
	#there were no fits calculated for them
	if (! early_slope_calculated) {
		results$early$R_sq[1] = NA
	}	
	if (! late_slope_calculated) {
		results$late$R_sq[1] = NA
	}
	
	best_results = list()
				
	best_results$early = as.data.frame(results$early)[best_indexes[1],]
	best_results$early$residual = resid$early[best_indexes[1]]
	best_results$late = as.data.frame(results$late)[best_indexes[2],]
	best_results$late$residual = resid$late[best_indexes[2]]

	best_results
}

find_best_offset_combination <- function(results, min_length = 10) {
	
	#Build an array with the sums of the collected R square values
	R_sq_sums = array(NA, c(length(results$early$R_sq),length(results$late$R_sq)));
	for (i in 1:length(results$early$R_sq)) {
		if (is.na(results$early$R_sq[i]) | is.nan(results$early$R_sq[i])) {
			next
		}
		for (j in 1:length(results$late$R_sq)) {
			if (is.na(results$late$R_sq[j]) | is.nan(results$late$R_sq[j])) {
				next
			}
			if ((j+i) > length(results$filt_init)) {
				next
			}
			
			R_sq_sums[i,j] = results$early$R_sq[i]+results$late$R_sq[j]
		}
	}
		
	#locate positions of the highest R squared value
	max_R_sq = max(R_sq_sums[! is.na(R_sq_sums)])
	if (max_R_sq == -Inf) {
		print(R_sq_sums)
	}
	
	highest_square_priority = -Inf
	
	for (i in 1:dim(R_sq_sums)[[1]]) {
		for (j in 1:dim(R_sq_sums)[[2]]) {
			if (is.na(R_sq_sums[i,j]) | is.nan(R_sq_sums[i,j])) {
				next
			}
			
			i_priority = i - min_length + 1
			if (i_priority <= 0) {
				i_priority = 0
			}
			j_priority = j - min_length + 1
			if (j_priority <= 0) {
				j_priority = 0
			}
			this_priority = i_priority^2 + j_priority^2
			
			if (  max_R_sq == R_sq_sums[i,j]
				& highest_square_priority < this_priority) {
				
				highest_square_priority = this_priority
				best_indexes = c(i,j)
			}
		}
	}
	best_indexes	
}

gather_linear_regions.boot <- function(results, 
	min_length = 10, col_lims = NaN, normed = 1, 
	log.trans = TRUE, boot.samp = NA) {

	sim_results <- list()
	
	#collect the entire set of adhesion signal values and lengths, excluding adhesions which don't live long enough
	all_ad_sig = c()
	all_length = c()
	for (i in 1:dim(results$exp_data)[[2]]) {
		temp = as.numeric(results$exp_data[i,])
		temp = temp[! is.nan(temp)]
		if (length(temp) >= min_length * 2) {
			all_length = c(all_length,length(temp))
			all_ad_sig = c(all_ad_sig,temp)	
		}
	}
	
	#produce fake adhesion signal and props variables by sampling from the values collected above
	sim_ad_sig <- array(NaN, dim = c(boot.samp, max(all_length) + 2))
	sim_props = list()
	for (i in 1:boot.samp) {
		this_length = sample(all_length,1)
		data = sample(all_ad_sig, this_length, replace = TRUE)
		sim_ad_sig[i,] <- c(NaN, data, array(NaN, dim = c(max(all_length) - length(data) + 1)))
		sim_props$death_status[i] = 1
	}
	
	sim_results <- gather_linear_regions(sim_ad_sig, sim_props, 
				       min_length = min_length, normed = normed, log.trans = log.trans, save.exp_data = FALSE)
}	

gather_correlations_from_dirs <- function (dirs, results, data_file='Area.csv',
	result.normed = TRUE, exp.normed = TRUE, result.log.trans = TRUE, 
	exp.log.trans = TRUE, results.file = NA, save.exp_data = TRUE) {
	
	corr_results = list()
	
	for (k in 1:length(dirs)) {
		if (is.na(dirs[[k]])) {
			next
		}
		
		print(dirs[[k]])
		
		exp_data <- as.matrix(read.table(file.path(dirs[[k]],data_file),header = FALSE, sep  = ','));
		
		corr_results[[k]] <- gather_correlations(results[[k]], exp_data, 
			result.normed = result.normed, exp.normed = exp.normed, 
			result.log.trans = result.log.trans, exp.log.trans = exp.log.trans, 
			results.file = results.file, save.exp_data = save.exp_data);
			
        if (! is.na(results.file)) {
            this_result = corr_results[[k]]
            save(this_result,file = file.path(dirs[[k]],results.file))
        }
	}
	
	corr_results
}

gather_correlations <- function(result, exp_data, result.normed = TRUE, 
	exp.normed = FALSE, result.log.trans = TRUE, exp.log.trans = FALSE, 
	results.file = NA, save.exp_data = TRUE) {

	corr_result = list()
	
	if (save.exp_data) {
		corr_result$exp_data = exp_data
	}
	
	count = 0
	for (i in 1:length(result$early$R_sq)) {
		data_1 = as.numeric(result$exp_data[i,])
		data_1 = data_1[! is.nan(data_1)]
		data_2 = as.numeric(exp_data[i,])
		data_2 = data_2[! is.nan(data_2)]		
		
		corr_result$early[i] = NA
		corr_result$late[i] = NA

		if (! is.na(result$early$R_sq[i])) {
			this_data_1 = data_1[1:result$early$offset[i]]
			this_data_2 = data_2[1:result$early$offset[i]]
			
			if (result.normed) {
				this_data_1 = this_data_1/this_data_1[1]
			}
			if (exp.normed) {
				this_data_2 = this_data_2/this_data_2[1]
			}
			
			if (result.log.trans) {
				this_data_1 = log(this_data_1)
			}
			if (exp.log.trans) {
				this_data_2 = log(this_data_2)
			}
			
			if (sd(this_data_1) != 0 & sd(this_data_2) != 0) {
				corr_result$early[i] = as.numeric(cor.test(this_data_1,this_data_2,use="all.obs")$estimate)
				corr_result$conf$early_lower[i] = as.numeric(cor.test(this_data_1,this_data_2,use="all.obs")$conf.int)[1]
				corr_result$conf$early_upper[i] = as.numeric(cor.test(this_data_1,this_data_2,use="all.obs")$conf.int)[2]
			}
		}
		if (! is.na(result$late$R_sq[i])) {
			this_data_1 = data_1[(length(data_1) - result$late$offset[i]):length(data_1)]
			this_data_2 = data_2[(length(data_2) - result$late$offset[i]):length(data_2)]

			if (result.normed) {
				this_data_1 = this_data_1[1]/this_data_1
			}
			if (exp.normed) {
				this_data_2 = this_data_2[1]/this_data_2
			}
			
			if (result.log.trans) {
				this_data_1 = log(this_data_1)
			}
			if (exp.log.trans) {
				this_data_2 = log(this_data_2)
			}
				
			if (sd(this_data_1) != 0 & sd(this_data_2) != 0) {
				corr_result$late[i] = as.numeric(cor.test(this_data_1,this_data_2,use="all.obs")$estimate)
				corr_result$conf$late_lower[i] = as.numeric(cor.test(this_data_1,this_data_2,use="all.obs")$conf.int)[1]
				corr_result$conf$late_upper[i] = as.numeric(cor.test(this_data_1,this_data_2,use="all.obs")$conf.int)[2]
			}
		}
	}
	corr_result
}

########################################
#Plotting Functions
########################################

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
    pdf(file.path(dir,file),width=8.5,height=8.5,pointsize=14);
    par(mfrow=c(2,2),bty='n', mar=c(5,4,1,1))
	
	#########################################################################
	#Plot 1 - Slope versus R squared (early)
	#########################################################################	
    plot(results$early_slope[results$early_R_sq != 0],
         results$early_R_sq[results$early_R_sq != 0],
         xlab='Formation Slope', ylab='R Squared', cex = 0.5,
         ylim = c(0,1)
        );

	#########################################################################
	#Plot 2 - R squared versus Slope (late)
	#########################################################################	
    plot(results$late_slope[results$late_R_sq != 0 & results$exp_props$death_status],
		 results$late_R_sq[results$late_R_sq != 0 & results$exp_props$death_status],
         xlab='Decay Slope (/min)', ylab='R Squared', cex = 0.5,
         ylim = c(0,1)
        )
    
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
	    pdf(file.path(dir,hist_file),width=8.5,height=8.5,pointsize=14);
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
		results = list(results)
	}
	for (i in 1:length(results)) {
		res = results[[i]];
		
		range = 1*(max(res$exp_props$starting_edge_dist) - min(res$exp_props$starting_edge_dist))
		
		death_status = res$exp_props$death_status[1:length(res$late_R_sq)]
		range_status = res$exp_props$starting_edge_dist[1:length(res$late_R_sq)]

		early_slopes <- c(early_slopes, res$early$slope[res$early$R_sq > r_cutoff & ! is.na(res$early$R_sq) & range_status < range])
		late_slopes  <- c(late_slopes, res$late$slope[res$late$R_sq > r_cutoff & ! is.na(res$late$R_sq) & death_status])
	}
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

plot_ad_seq <- function (results,index,type='early',...) {
	ad_seq = as.vector(results$exp_data[index,])
	ad_seq = t(ad_seq[!(is.nan(ad_seq))])
	
	if (type == 'early') {
	this_ad_seq = ad_seq[1:results$early$offset[index]];
	this_ad_seq = log(this_ad_seq/this_ad_seq[1]);
		
	x = c(0,results$early$offset[index]);
	y = c(results$early$slope[index]*x[1] + results$early$inter[index],
		   	  results$early$slope[index]*x[2] + results$early$inter[index])
		
	plot(1:results$early$offset[index],this_ad_seq,xlab='Time (minutes)',ylab='ln(Intensity/First Intensity)',
			 ylim=c(min(this_ad_seq,y),max(this_ad_seq,y)))
		
	lines(x,y,col='red',lwd=2)
	r_sq_val_str = sprintf('%.3f',results$early$R_sq[index])
	slope_val_str = sprintf('%.3f',results$early$slope[index])
	exp_str = paste('R^2=',r_sq_val_str,'\n Slope = ',slope_val_str,sep='')
	text(x[1]+3,0.5*max(this_ad_seq),paste('R^2 = ',sprintf('%.3f',results$early$R_sq[index]),'\n Slope = ',sprintf('%.3f',results$early$slope[index]),sep=''))
	}

	if (type == 'late') {
	this_ad_seq = ad_seq[(length(ad_seq) - results$late$offset[index]) : length(ad_seq)];
	this_ad_seq = log(this_ad_seq[1]/this_ad_seq);

	x = c(length(ad_seq) - results$late$offset[index],length(ad_seq));
	y = c(results$late$slope[index]*x[1] + results$late$inter[index],
		   	  results$late$slope[index]*x[2] + results$late$inter[index])
		
	plot((length(ad_seq) - results$late$offset[index]) : length(ad_seq),
		 this_ad_seq, xlab='Time (minutes)', ylab='ln(First Intensity/Intensity)',
		 ylim=c(min(this_ad_seq,y),max(this_ad_seq,y)))
		
	lines(x,y,col='red',lwd=2)
	text(x[1]+ 3,0.5*max(this_ad_seq),paste('R^2 = ',sprintf('%.3f',results$late$R_sq[index]),'\n Slope = ',sprintf('%.3f',results$late$slope[index]),sep=''))
	}
}

plot_overall_residuals <- function(results,dir,file='overall_residual_plot.pdf',window = 0.1) {
	
	resid = gather_exp_residuals(results)
	
	resid_win = gather_exp_win_residuals(resid,window=window)
	
	library(Hmisc)
	pdf(file = file.path(dir,file),width=12)
	par(mfcol=c(1,2))
	errbar(resid_win$x$early, resid_win$y$early, 
		   resid_win$y$early - resid_win$err$early, resid_win$y$early + resid_win$err$early, 
		   xlab='Scaled Linear Fit Position',ylab='Residuals')
	lines(c(0,1),c(0,0))
	
	errbar(resid_win$x$late, resid_win$y$late, 
		   resid_win$y$late - resid_win$err$late, resid_win$y$late + resid_win$err$late, 
		   xlab='Scaled Linear Fit Position',ylab='Residuals')
	lines(c(0,1),c(0,0))
	dev.off()
}

gather_exp_residuals <- function(results, min_R_sq = NA) {
	resid = list()
	for (i in 1:length(results)) {
		res = results[[i]]
		for (j in 1:length(res$early$R_sq)) {
			if (is.na(res$early$R_sq[j])) {	
				next
			}
			if (is.numeric(min_R_sq) & (res$early$R_sq[j] < min_R_sq)) {
				next
			}

			resid_list = res$early$residual[j][[1]]
			x = list(seq(0,1,1/(length(resid_list)-1)))
			resid$y$early = c(resid$y$early,resid_list)
			resid$x$early = c(resid$x$early,x)
		}
		for (j in 1:length(res$late$R_sq)) {
			if (is.na(res$late$R_sq[j])) {
				next
			}
			if (is.numeric(min_R_sq) & (res$late$R_sq[j] < min_R_sq)) {
				next
			}

			resid_list = res$late$residual[j][[1]]
			x = list(seq(0,1,1/(length(resid_list)-1)))
			resid$y$late = c(resid$y$late,resid_list)
			resid$x$late = c(resid$x$late,x)
		}
	}
	resid
}

gather_exp_win_residuals <- function(resid, window) {
	resid_win = list()
	for (i in seq(window,1,by=window)) {
		temp = c()
		for (j in 1:length(resid$x$early)) {
			this_resid_x = resid$x$early[[j]]
			for (k in 1:length(this_resid_x)) {
				if (this_resid_x[[k]] <= i & this_resid_x[[k]] >= i - window) {
					temp = c(temp,resid$y$early[[j]])
				}
			}
		}
		resid_win$x$early   = c(resid_win$x$early,i-(window/2))
		resid_win$y$early   = c(resid_win$y$early,mean(temp))
		resid_win$err$early = c(resid_win$err$early,sd(temp))
		
		temp = c()
		for (j in 1:length(resid$x$late)) {
			this_resid_x = resid$x$late[[j]]
			for (k in 1:length(this_resid_x)) {
				if (this_resid_x[[k]] <= i & this_resid_x[[k]] >= i - window) {
					temp = c(temp,resid$y$late[[j]])
				}
			}
		}
		resid_win$x$late   = c(resid_win$x$late,i-(window/2))
		resid_win$y$late   = c(resid_win$y$late,mean(temp))
		resid_win$err$late = c(resid_win$err$late,sd(temp))
	}
	resid_win
}

boxplot_with_points <- function(data, 
	colors=c('red','green','yellow','blue','pink','cyan','gray','orange','brown','purple'), 
	notch=F, names, range=1.5, ...) {
		
	box.data = boxplot(data,notch = notch,names = names,varwidth=T,range = range,...)
	for (i in 1:length(data)) {
		this_data = data[[i]]
		temp_data = this_data[this_data >= box.data$stat[1,i] & this_data <= box.data$stat[5,i]]
		points(jitter(array(0,dim=c(1,length(temp_data))),10)+i,temp_data,col=colors[[i]])
	}
}

hist_with_percents <- function(data, ...) {
	hist_data = hist(data, ...);

	for (i in 1:length(hist_data$counts)) {
		y_pos = 0;

		if ((hist_data$counts[i] + 0.05*max(hist_data$counts)) > 0.5*max(hist_data$counts)) {
			y_pos = hist_data$counts[i] - 0.05*max(hist_data$counts);
		} else {
			y_pos = hist_data$counts[i] + 0.05*max(hist_data$counts);
		}
		text(hist_data$mids[i], y_pos, sprintf('%.03f',hist_data$counts[i]/sum(hist_data$counts)), srt = 45)
	}
}


########################################
#Misc functions
########################################

filter_results <- function(results,needed_R_sq=0.9) {
	points = list()
	for (i in 1:length(results)) {
		res = results[[i]]

		early_filt = is.finite(res$early$R_sq) & res$early$R_sq > needed_R_sq & is.finite(res$early$slope) & res$early$slope > 0
		late_filt = is.finite(res$late$R_sq) & res$late$R_sq > needed_R_sq & is.finite(res$late$slope) & res$late$slope > 0
	
		points$early$slope = c(points$early$slope,res$early$slope[early_filt])
		points$late$slope = c(points$late$slope,res$late$slope[late_filt])
		points$early$R_sq = c(points$early$R_sq, res$early$R_sq[early_filt])
		points$late$R_sq = c(points$late$R_sq, res$late$R_sq[late_filt])

		points$early$stable_lifetime = c(points$early$stable_lifetime, res$stable_lifetime[early_filt])
		points$late$stable_lifetime = c(points$late$stable_lifetime, res$stable_lifetime[late_filt])
		points$early$stable_mean = c(points$early$stable_mean, res$stable_mean[early_filt])
		points$late$stable_mean = c(points$late$stable_mean, res$stable_mean[late_filt])
		points$early$stable_variance = c(points$early$stable_variance, res$stable_variance[early_filt])
		points$late$stable_variance = c(points$late$stable_variance, res$stable_variance[late_filt])

		points$early$lin_num = c(points$early$lin_num,which(early_filt))
		points$late$lin_num = c(points$late$lin_num,which(late_filt))
		
		points$early$exp_dir = c(points$early$exp_dir,rep(res$exp_dir,length(which(early_filt))))
		points$late$exp_dir = c(points$late$exp_dir,rep(res$exp_dir,length(which(late_filt))))
		
		points$ind_exp[[i]] = list(early_slope = res$early$slope[early_filt],
								   late_slope = res$late$slope[late_filt])
		
		if (any(names(res$exp_props) == 'starting_edge_dist')) {
			points$early$edge_dist = c(points$early$edge_dist,res$exp_props$starting_edge_dist[early_filt])
		}
		
		if (any(names(res$exp_props) == 'ending_edge_dist')) {
			points$late$edge_dist = c(points$late$edge_dist,res$exp_props$ending_edge_dist[late_filt])
		}
	}
	points$early = as.data.frame(points$early)
	points$late = as.data.frame(points$late)
	
	points
}


load_results <- function(dirs,file) {
	results = list()
	for (i in 1:length(dirs)) {
		this_file = file.path(dirs[i],file)
		if (file.exists(this_file)) {
			load(file.path(dirs[i],file))
			results[[i]] = this_result
		}
	}
	results
}

trim_args_list <- function(args) {
	for (i in 1:length(args)) {
		if (substr(args[i],0,1) == '-') {
			args[i] = substring(args[i],2)
		}
	}
	args
}

write_high_r_rows <- function(result, dir, file=c('early_R_sq.csv','late_R_sq.csv','neg_slope_R_sq.csv'), min_R_sq = 0.9) {
	if (! file.exists(dir)) {
		dir.create(dir,recursive=TRUE)
	}
	
	row_nums = which(is.finite(result$early$R_sq) & result$early$R_sq > min_R_sq)
	if (! is.null(row_nums)) {
		write.table(row_nums,file=file.path(dir,file[1]), row.names=FALSE, col.names=FALSE)
	}
	
	row_nums = which(is.finite(result$early$R_sq) & result$late$R_sq > min_R_sq)
	if (! is.null(row_nums)) {
		write.table(row_nums,file=file.path(dir,file[2]), row.names=FALSE, col.names=FALSE)
	}
	
	row_nums = which(is.finite(result$early$R_sq) & result$early$R_sq > min_R_sq & result$early$slope < 0)
	if (! is.null(row_nums)) {
		write.table(row_nums,file=file.path(dir,file[3]), row.names=FALSE, col.names=FALSE)
	}
}

gather_datafile_from_dirs <- function (dirs, data_file='Average_adhesion_signal.csv') {
	exp_data = list()
	for (k in 1:length(dirs)) {
		exp_data[[k]] <- read.table(file.path(dirs[[k]],data_file),header = TRUE, sep  = ',');
	}
	exp_data
}

get_legend_rect_points <- function(left_x,bottom_y,right_x,top_y,box_num) {
	left_x_seq = array(left_x,box_num)
	right_x_seq = array(right_x,box_num)
	bottom_y_seq = c()
	top_y_seq = c()
	
	for (i in 1:box_num) {
		bottom_y_seq = c(bottom_y_seq, (top_y - bottom_y)*((i-1)/11)+bottom_y)
		top_y_seq = c(top_y_seq,(top_y - bottom_y)*(i/11)+bottom_y)
	}
	rbind(left_x_seq,bottom_y_seq,right_x_seq,top_y_seq)
}	

find_birth_death_rate <- function(results) {
	birth_rate = list()
	death_rate = list()
	for (i in 1:length(results)) {
		birth_rate[[i]] = rep(0,dim(results[[i]]$exp_data)[[2]]-1);
		death_rate[[i]] = rep(0,dim(results[[i]]$exp_data)[[2]]-1);
		for (j in 1:dim(results[[i]]$exp_data)[[1]]) {
			data_line = is.na(as.vector(results[[i]]$exp_data[j,]));
			true_line = which(!data_line);
#		print(data_line)
			if (length(true_line) == 0) {
				next;
			}
			if (data_line[1]) {
				birth_rate[[i]][true_line[1]-1] = birth_rate[[i]][true_line[1]-1] + 1;
			}
			if (data_line[length(data_line)]) {
				death_rate[[i]][true_line[length(true_line)]] = death_rate[[i]][true_line[length(true_line)]] + 1;
			}		
		}
	}

	for (i in 1:length(birth_rate)) {
	#	print(cor(birth_rate[[i]],death_rate[[i]]))
	#	print(mean(birth_rate[[i]]-death_rate[[i]]))
	}
}

################################################################################
# Main Program
################################################################################

args <- commandArgs(TRUE)

if (length(args) != 0) {
	args <- trim_args_list(args)
	
	ave_results = gather_bilinear_models_from_dirs(args, 
		results.file=file.path('..','intensity_model.Rdata'))
		
	gather_bilinear_models_from_dirs(args, 
		data_file='Background_corrected_signal.csv', 
		results.file=file.path('..','corrected_intensity_model.Rdata'))

	gather_bilinear_models_from_dirs(args, 
		data_file='Shrunk_corrected_signal.csv', 
		results.file=file.path('..','shrunk_intensity_model.Rdata'))
		
	write_high_r_rows(ave_results[[1]],file.path(args[1],'..','for_vis'))
}