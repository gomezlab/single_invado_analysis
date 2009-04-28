################################################################################
#linear_regions.R: various functions used to find and plot the linear regions
#  and associated data from the focal adhesion identification/analysis programs
################################################################################

########################################
#Data fitting functions
########################################

gather_bilinear_models_from_dirs <- function (dirs, min_length = 10,
	data_file='Average_adhesion_signal.csv', col_lims = NA, 
	normed = TRUE, log.trans = TRUE, boot.samp = NA, results.file = NA,
	save.exp_data = TRUE, debug = FALSE) {
	
	results = list()
	
	for (i in 1:length(dirs)) {
		if (is.na(dirs[[i]])) {
			next;
		}
		if (! file.exists(file.path(dirs[[i]],data_file))) {
			next;
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
							boot.samp = boot.samp, save.exp_data = save.exp_data, debug=debug);
        regex_range = regexpr("time_series_[[:digit:]]",dirs[[i]])
        if (regex_range[1] == -1) {
        	results[[i]]$exp_dir = dirs[[i]];
        } else {
        	results[[i]]$exp_dir = substr(dirs[[i]], regex_range[1], regex_range[1] + attr(regex_range,'match.length'));
        }
        
        if (! is.na(results.file)) {
        	this_file = file.path(dirs[[i]],results.file)
        	if (! file.exists(dirname(this_file))) {
				dir.create(dirname(this_file),recursive=TRUE)
			}
            this_result = results[[i]];
            save(this_result,file = this_file);
        }
	}

	return(results);
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

		results$assembly$offset[i] = temp_results$assembly$offset
		results$assembly$R_sq[i]   = temp_results$assembly$R_sq
		results$assembly$p_val[i]   = temp_results$assembly$p_val
		results$assembly$inter[i]  = temp_results$assembly$inter
		results$assembly$slope[i]  = temp_results$assembly$slope
		results$assembly$fold_change[i]  = temp_results$assembly$fold_change
		results$assembly$residual[i]  = temp_results$assembly$residual
		
		results$disassembly$offset[i]  = temp_results$disassembly$offset
		results$disassembly$R_sq[i]    = temp_results$disassembly$R_sq
		results$disassembly$p_val[i]    = temp_results$disassembly$p_val
		results$disassembly$inter[i]   = temp_results$disassembly$inter
		results$disassembly$slope[i]   = temp_results$disassembly$slope
		results$disassembly$fold_change[i]  = temp_results$disassembly$fold_change
		results$disassembly$residual[i]   = temp_results$disassembly$residual
		
		#if either of the offset values are NA, then we weren't able 
		#to get a fit for one side of the data, don't calcudisassembly a 
		#stable lifetime because we don't know how long the adhesion 
		#was around before the movie starts or after it ends
		if (! is.na(results$disassembly$offset[i]) && ! is.na(results$assembly$offset[i])) {
			results$stable_data_set[[i]] = numeric_data_set[results$assembly$offset[i]:(length(numeric_data_set) - results$disassembly$offset[i])]
			
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
	if (length(results$assembly) > 0) {
		for (i in 1:length(results$assembly)) {
			for (j in (length(results$assembly[[i]]) + 1):desired_length) {
				results$assembly[[i]][j] = NA
			}
		}
	}
	if (length(results$disassembly) > 0) {
		for (i in 1:length(results$disassembly)) {
			for (j in (length(results$disassembly[[i]]) + 1):desired_length) {
				results$disassembly[[i]][j] = NA
			}
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
	resid = list(assembly = list(), disassembly = list())	
	results$filt_init = initial_data_set[! is.nan(initial_data_set)]
	if (length(results$filt_init) == 0) {
		return(c(NA, NA))
	}
	this_data_set = data.frame(y = results$filt_init, x = 1:length(results$filt_init))

	#Search the beginning of the sequence for a linear fit
	assembly_slope_calculated = FALSE;
	if (is.nan(initial_data_set[1]) & ! exp_props$split_birth_status) {
		assembly_slope_calculated = TRUE;
		for (j in min_length:dim(this_data_set)[[1]]) {
			assembly_subset = this_data_set[1:j,]
			if (normed) {
				assembly_subset$y = assembly_subset$y/assembly_subset$y[1]
			}
			if (log.trans) {
				assembly_subset$y = log(assembly_subset$y)
			}
				
			model <- lm(y ~ x, data = assembly_subset)
			summary <- summary(model);
			
			results$assembly$R_sq[j] = summary$adj.r.squared
			#dealing with a degerate case, where lm produce NaN for the R squared 
			#value when the data set is a flat line, see:
			#	>data <- data.frame(x = c(1,2,3), y = c(1,1,1))
			#	>summary(lm(y ~ x, data=data))
			if (is.nan(results$assembly$R_sq[j])) {
				results$assembly$R_sq[j] = 1
			}
			
			results$assembly$p_val[j] = summary$coefficients[2,4]
			results$assembly$length[j] = dim(assembly_subset)[[1]]
			results$assembly$offset[j] = j
			results$assembly$inter[j] = coef(model)[[1]]
			results$assembly$slope[j] = coef(model)[[2]]

			if (log.trans) {
				results$assembly$fold_change[j] = max(assembly_subset$y)
			} else {				
				results$assembly$fold_change[j] = max(assembly_subset$y)/min(assembly_subset$y)
			}
			resid$assembly[[j]] = as.numeric(resid(model))
		}
	} else {
		results$assembly$R_sq[1] = 0

		results$assembly$p_val[1] = NA		
		results$assembly$length[1] = NA
		results$assembly$offset[1] = NA
		results$assembly$inter[1] = NA
		results$assembly$slope[1] = NA
		results$assembly$fold_change[1] = NA

		resid$assembly[[1]] = NA
	}
	
	#Search the end of the sequence for a linear fit
	disassembly_slope_calculated = FALSE;
	if (is.nan(initial_data_set[length(initial_data_set)]) & exp_props$death_status) {
		disassembly_slope_calculated = TRUE;
		for (j in min_length:dim(this_data_set)[[1]]) {
			disassembly_subset = this_data_set[(dim(this_data_set)[[1]]-j):dim(this_data_set)[[1]],]
			if (normed) {
				disassembly_subset$y = disassembly_subset$y[1]/disassembly_subset$y
			}
			if (log.trans) {
				disassembly_subset$y = log(disassembly_subset$y)
			}
	
			model <- lm(y ~ x, data = disassembly_subset)
			summary <- summary(model);
			
			results$disassembly$R_sq[j] = summary$adj.r.squared
			#dealing with a degerate case, where lm produce NaN for the R squared 
			#value when the data set is a flat line, see:
			#	>data <- data.frame(x = c(1,2,3), y = c(1,1,1))
			#	>summary(lm(y ~ x, data=data))
			if (is.nan(results$disassembly$R_sq[j])) {
				results$disassembly$R_sq[j] = 1
			}
			
			results$disassembly$p_val[j] = summary$coefficients[2,4]
			results$disassembly$length[j] = dim(disassembly_subset)[[1]]
			results$disassembly$offset[j] = j
			results$disassembly$inter[j] = coef(model)[[1]]
			results$disassembly$slope[j] = coef(model)[[2]]

			if (log.trans) {
				results$disassembly$fold_change[j] = max(disassembly_subset$y)
			} else {						
				results$disassembly$fold_change[j] = max(disassembly_subset$y)/min(disassembly_subset$y)
			}
			resid$disassembly[[j]] = as.numeric(resid(model))
		}
	} else {
		results$disassembly$R_sq[1] = 0
		
		results$disassembly$p_val[1] = NA	
		results$disassembly$length[1] = NA
		results$disassembly$offset[1] = NA
		results$disassembly$inter[1] = NA
		results$disassembly$slope[1] = NA
		results$disassembly$fold_change[1] = NA

		resid$disassembly[[1]] = NA	
	}

	best_indexes = find_best_offset_combination(results, min_length = min_length)
	
	#With the R squared matrix calcudisassemblyd reset the r_sq componenets to NA, if needed since 
	#there were no fits calcudisassemblyd for them
	if (! assembly_slope_calculated) {
		results$assembly$R_sq[1] = NA
	}	
	if (! disassembly_slope_calculated) {
		results$disassembly$R_sq[1] = NA
	}
	
	best_results = list()
				
	best_results$assembly = as.data.frame(results$assembly)[best_indexes[1],]
	best_results$assembly$residual = resid$assembly[best_indexes[1]]
	best_results$disassembly = as.data.frame(results$disassembly)[best_indexes[2],]
	best_results$disassembly$residual = resid$disassembly[best_indexes[2]]

	best_results
}

find_best_offset_combination <- function(results, min_length = 10) {
	
	#Build an array with the sums of the collected R square values
	R_sq_sums = array(NA, c(length(results$assembly$R_sq),length(results$disassembly$R_sq)));
	for (i in 1:length(results$assembly$R_sq)) {
		if (is.na(results$assembly$R_sq[i]) | is.nan(results$assembly$R_sq[i])) {
			next
		}
		for (j in 1:length(results$disassembly$R_sq)) {
			if (is.na(results$disassembly$R_sq[j]) | is.nan(results$disassembly$R_sq[j])) {
				next
			}
			if ((j+i) > length(results$filt_init)) {
				next
			}
			
			R_sq_sums[i,j] = results$assembly$R_sq[i]+results$disassembly$R_sq[j]
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
			this_priority = i_priority + j_priority
			
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
	for (i in 1:length(result$assembly$R_sq)) {
		data_1 = as.numeric(result$exp_data[i,])
		data_1 = data_1[! is.nan(data_1)]
		data_2 = as.numeric(exp_data[i,])
		data_2 = data_2[! is.nan(data_2)]		
		
		corr_result$assembly[i] = NA
		corr_result$disassembly[i] = NA

		if (! is.na(result$assembly$R_sq[i])) {
			this_data_1 = data_1[1:result$assembly$offset[i]]
			this_data_2 = data_2[1:result$assembly$offset[i]]
			
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
				corr_result$assembly[i] = as.numeric(cor.test(this_data_1,this_data_2,use="all.obs")$estimate)
				corr_result$conf$assembly_lower[i] = as.numeric(cor.test(this_data_1,this_data_2,use="all.obs")$conf.int)[1]
				corr_result$conf$assembly_upper[i] = as.numeric(cor.test(this_data_1,this_data_2,use="all.obs")$conf.int)[2]
			}
		}
		if (! is.na(result$disassembly$R_sq[i])) {
			this_data_1 = data_1[(length(data_1) - result$disassembly$offset[i]):length(data_1)]
			this_data_2 = data_2[(length(data_2) - result$disassembly$offset[i]):length(data_2)]

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
				corr_result$disassembly[i] = as.numeric(cor.test(this_data_1,this_data_2,use="all.obs")$estimate)
				corr_result$conf$disassembly_lower[i] = as.numeric(cor.test(this_data_1,this_data_2,use="all.obs")$conf.int)[1]
				corr_result$conf$disassembly_upper[i] = as.numeric(cor.test(this_data_1,this_data_2,use="all.obs")$conf.int)[2]
			}
		}
	}
	corr_result
}

########################################
#Plotting Functions
########################################

plot_lin_reg_set <- function(results,dir,file='linear_regions.pdf', hist_file=NA) {
	assembly_slope = c()
	assembly_error = c()
	assembly_n = c()
	disassembly_slope  = c()
	disassembly_error  = c()
	disassembly_n = c()
	for (i in 1:10) {
		slope_ests = exp_set_slope_estimate(results,r_cutoff = ((i - 1)/10));
		assembly_slope[i] = slope_ests$assembly
		assembly_n[i] = slope_ests$assembly_n
		assembly_error[i] = slope_ests$assembly_sem
		disassembly_slope[i] = slope_ests$disassembly
		disassembly_n[i] = slope_ests$disassembly_n
		disassembly_error[i] = slope_ests$disassembly_sem
	}
	
	library(Hmisc)
    pdf(file.path(dir,file),width=8.5,height=8.5,pointsize=14);
    par(mfrow=c(2,2),bty='n', mar=c(5,4,1,1))
	
	#########################################################################
	#Plot 1 - Slope versus R squared (assembly)
	#########################################################################	
    plot(results$assembly_slope[results$assembly_R_sq != 0],
         results$assembly_R_sq[results$assembly_R_sq != 0],
         xlab='Formation Slope', ylab='R Squared', cex = 0.5,
         ylim = c(0,1)
        );

	#########################################################################
	#Plot 2 - R squared versus Slope (disassembly)
	#########################################################################	
    plot(results$disassembly_slope[results$disassembly_R_sq != 0 & results$exp_props$death_status],
		 results$disassembly_R_sq[results$disassembly_R_sq != 0 & results$exp_props$death_status],
         xlab='Decay Slope (/min)', ylab='R Squared', cex = 0.5,
         ylim = c(0,1)
        )
    
	#Plot 3
    errbar(seq(0,0.9,by=0.1), assembly_slope, assembly_slope - assembly_error, assembly_slope + assembly_error, xlab = 'R squared cutoff', ylab='Accumulation Rate (/min)')
    for (i in 1:10) {
    	x_points = seq(0,0.9,by=0.1)
    	text(x_points[i],assembly_slope[i]+assembly_error[i],assembly_n[i])
    }
	
	#Plot 4
    errbar(seq(0,0.9,by=0.1), disassembly_slope, disassembly_slope - disassembly_error, disassembly_slope + disassembly_error, xlab = 'R squared cutoff', ylab='Decay Rate (/min)')
    for (i in 1:10) {
    	x_points = seq(0,0.9,by=0.1)
    	text(x_points[i],disassembly_slope[i]+disassembly_error[i],disassembly_n[i])
    }
	
	if (! is.na(hist_file)) {
	    pdf(file.path(dir,hist_file),width=8.5,height=8.5,pointsize=14);
    	par(mfrow=c(2,2),bty='n', mar=c(5,4,1,1))
        
	    #Plot 1
    	hist(results$assembly_length, xlab = 'Linear Sequence Length (min)', main='')
    
	    #Plot 2
    	hist(results$disassembly_length, xlab = 'Linear Sequence Length (min)', main='')
    
	    #Plot 3
    	hist(results$assembly_slope[results$assembly_R_sq > 0.9], xlab = 'Slope (/min)', main='')

	    #Plot 4
    	hist(results$disassembly_slope[results$disassembly_R_sq > 0.9 & results$exp_props$death_status], xlab = 'Slope (/min)', main='')
    	dev.off()
    }
    
    dev.off()
}

exp_set_slope_estimate <- function(results,r_cutoff=0.9) {
	assembly_slopes <- c()
	disassembly_slopes <- c()
	
	if (! is.null(names(results))) {
		results = list(results)
	}
	for (i in 1:length(results)) {
		res = results[[i]];
		
		range = 1*(max(res$exp_props$starting_edge_dist) - min(res$exp_props$starting_edge_dist))
		
		death_status = res$exp_props$death_status[1:length(res$disassembly_R_sq)]
		range_status = res$exp_props$starting_edge_dist[1:length(res$disassembly_R_sq)]

		assembly_slopes <- c(assembly_slopes, res$assembly$slope[res$assembly$R_sq > r_cutoff & ! is.na(res$assembly$R_sq) & range_status < range])
		disassembly_slopes  <- c(disassembly_slopes, res$disassembly$slope[res$disassembly$R_sq > r_cutoff & ! is.na(res$disassembly$R_sq) & death_status])
	}
	slopes <- list(assembly = mean(assembly_slopes),
				   assembly_n = length(assembly_slopes),
				   assembly_cv = sd(assembly_slopes)/mean(assembly_slopes),
				   assembly_sem = sd(assembly_slopes)/sqrt(length(assembly_slopes)),
				   assembly_sd = sd(assembly_slopes),
				   disassembly = mean(disassembly_slopes),
				   disassembly_n = length(disassembly_slopes),
				   disassembly_cv = sd(disassembly_slopes)/mean(disassembly_slopes),
				   disassembly_sem = sd(disassembly_slopes)/sqrt(length(disassembly_slopes)),
				   disassembly_sd = sd(disassembly_slopes)
				  )
}

plot_ad_seq <- function (results,index,type='assembly',...) {
	ad_seq = as.vector(results$exp_data[index,])
	ad_seq = t(ad_seq[!(is.nan(ad_seq))])
	
	if (type == 'assembly') {
		this_ad_seq = ad_seq[1:results$assembly$offset[index]];
		this_ad_seq = log(this_ad_seq/this_ad_seq[1]);
		
		x = c(0,results$assembly$offset[index]);
		y = c(results$assembly$slope[index]*x[1] + results$assembly$inter[index],
			  results$assembly$slope[index]*x[2] + results$assembly$inter[index])
		
		plot(1:results$assembly$offset[index],this_ad_seq,xlab='Time (minutes)',ylab='ln(Intensity/First Intensity)',
				 ylim=c(min(this_ad_seq,y),max(this_ad_seq,y)))
		
		lines(x,y,col='red',lwd=2)
		r_sq_val_str = sprintf('%.3f',results$assembly$R_sq[index])
		slope_val_str = sprintf('%.3f',results$assembly$slope[index])
		exp_str = paste('R^2=',r_sq_val_str,'\n Slope = ',slope_val_str,sep='')
		text(x[1]+3,0.5*max(this_ad_seq), 
			 paste('R^2 = ',sprintf('%.3f',results$assembly$R_sq[index]),'\n Slope = ',sprintf('%.3f',results$assembly$slope[index]),sep=''))
	}

	if (type == 'disassembly') {
		this_ad_seq = ad_seq[(length(ad_seq) - results$disassembly$offset[index]) : length(ad_seq)];
		this_ad_seq = log(this_ad_seq[1]/this_ad_seq);

		x = c(length(ad_seq) - results$disassembly$offset[index],length(ad_seq));
		y = c(results$disassembly$slope[index]*x[1] + results$disassembly$inter[index],
		   	  results$disassembly$slope[index]*x[2] + results$disassembly$inter[index])
		
		plot((length(ad_seq) - results$disassembly$offset[index]) : length(ad_seq),
			 this_ad_seq, xlab='Time (minutes)', ylab='ln(First Intensity/Intensity)',
			 ylim=c(min(this_ad_seq,y),max(this_ad_seq,y)))
		
		lines(x,y,col='red',lwd=2)
		text(x[1]+ 3,0.5*max(this_ad_seq), 
			 paste('R^2 = ',sprintf('%.3f',results$disassembly$R_sq[index]),'\n Slope = ',sprintf('%.3f',results$disassembly$slope[index]),sep=''))
	}
	
	if (type == 'overall') {
		plot(0:(length(ad_seq)-1), ad_seq, xlab='Time (minutes)', ylab='Intensity',
			 ylim=c(0,1))
	}
}

plot_overall_residuals <- function(results,dir,file='overall_residual_plot.pdf',window = 0.1) {
	
	resid = gather_exp_residuals(results)
	
	resid_win = gather_exp_win_residuals(resid,window=window)
	
	library(Hmisc)
	pdf(file = file.path(dir,file),width=12)
	par(mfcol=c(1,2))
	errbar(resid_win$x$assembly, resid_win$y$assembly, 
		   resid_win$y$assembly - resid_win$err$assembly, resid_win$y$assembly + resid_win$err$assembly, 
		   xlab='Scaled Linear Fit Position',ylab='Residuals')
	lines(c(0,1),c(0,0))
	
	errbar(resid_win$x$disassembly, resid_win$y$disassembly, 
		   resid_win$y$disassembly - resid_win$err$disassembly, resid_win$y$disassembly + resid_win$err$disassembly, 
		   xlab='Scaled Linear Fit Position',ylab='Residuals')
	lines(c(0,1),c(0,0))
	dev.off()
}

gather_exp_residuals <- function(results, min_R_sq = NA) {
	resid = list()
	for (i in 1:length(results)) {
		res = results[[i]]
		for (j in 1:length(res$assembly$R_sq)) {
			if (is.na(res$assembly$R_sq[j])) {	
				next
			}
			if (is.numeric(min_R_sq) & (res$assembly$R_sq[j] < min_R_sq)) {
				next
			}

			resid_list = res$assembly$residual[j][[1]]
			x = list(seq(0,1,1/(length(resid_list)-1)))
			resid$y$assembly = c(resid$y$assembly,resid_list)
			resid$x$assembly = c(resid$x$assembly,x)
		}
		for (j in 1:length(res$disassembly$R_sq)) {
			if (is.na(res$disassembly$R_sq[j])) {
				next
			}
			if (is.numeric(min_R_sq) & (res$disassembly$R_sq[j] < min_R_sq)) {
				next
			}

			resid_list = res$disassembly$residual[j][[1]]
			x = list(seq(0,1,1/(length(resid_list)-1)))
			resid$y$disassembly = c(resid$y$disassembly,resid_list)
			resid$x$disassembly = c(resid$x$disassembly,x)
		}
	}
	resid
}

gather_exp_win_residuals <- function(resid, window) {
	resid_win = list()
	for (i in seq(window,1,by=window)) {
		temp = c()
		for (j in 1:length(resid$x$assembly)) {
			this_resid_x = resid$x$assembly[[j]]
			for (k in 1:length(this_resid_x)) {
				if (this_resid_x[[k]] <= i & this_resid_x[[k]] >= i - window) {
					temp = c(temp,resid$y$assembly[[j]])
				}
			}
		}
		resid_win$x$assembly   = c(resid_win$x$assembly,i-(window/2))
		resid_win$y$assembly   = c(resid_win$y$assembly,mean(temp))
		resid_win$err$assembly = c(resid_win$err$assembly,sd(temp))
		
		temp = c()
		for (j in 1:length(resid$x$disassembly)) {
			this_resid_x = resid$x$disassembly[[j]]
			for (k in 1:length(this_resid_x)) {
				if (this_resid_x[[k]] <= i & this_resid_x[[k]] >= i - window) {
					temp = c(temp,resid$y$disassembly[[j]])
				}
			}
		}
		resid_win$x$disassembly   = c(resid_win$x$disassembly,i-(window/2))
		resid_win$y$disassembly   = c(resid_win$y$disassembly,mean(temp))
		resid_win$err$disassembly = c(resid_win$err$disassembly,sd(temp))
	}
	resid_win
}

boxplot_with_points <- function(data, 
	colors=c('red','green','yellow','blue','pink','cyan','gray','orange','brown','purple'), 
	notch=F, names, range=1.5, inc.n.counts = TRUE, inc.points = TRUE, ...) {
		
	
	if (inc.n.counts) {
		for (i in 1:length(data)) {
			names[i] = paste(names[i], ' (n=', length(data[[i]]), ')', sep ='');
		}
	}
	
	box.data = boxplot(data,notch = notch,names = names,varwidth=T,range = range,...)
	if (inc.points) {
		for (i in 1:length(data)) {
			this_data = data[[i]]
			temp_data = this_data[this_data >= box.data$stat[1,i] & this_data <= box.data$stat[5,i]]
			points(jitter(array(0,dim=c(1,length(temp_data))),10)+i,temp_data,col=colors[[i]])
		}
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
filter_mixed_results <- function(results, corrected, min_R_sq=0.9, max_p_val = 0.05, pos_slope = TRUE) {
	points = list()
	for (i in 1:length(corrected)) {
		corr = corrected[[i]]
		res = results[[i]]

		assembly_filt = (is.finite(corr$assembly$R_sq) & corr$assembly$R_sq >= min_R_sq
					   & is.finite(corr$assembly$slope) & is.finite(corr$assembly$p_val) 
					   & corr$assembly$p_val <= max_p_val)
		disassembly_filt = (is.finite(corr$disassembly$R_sq) & corr$disassembly$R_sq >= min_R_sq 
						  & is.finite(corr$disassembly$slope) & is.finite(corr$disassembly$p_val) 
						  & corr$disassembly$p_val <= max_p_val)
		
		if (pos_slope) {
			assembly_filt = assembly_filt & corr$assembly$slope > 0
			disassembly_filt = disassembly_filt & corr$disassembly$slope > 0
		}				  
	
		points$assembly$slope = c(points$assembly$slope, corr$assembly$slope[assembly_filt])
		points$disassembly$slope = c(points$disassembly$slope, corr$disassembly$slope[disassembly_filt])
		points$assembly$R_sq = c(points$assembly$R_sq, corr$assembly$R_sq[assembly_filt])
		points$disassembly$R_sq = c(points$disassembly$R_sq, corr$disassembly$R_sq[disassembly_filt])
		points$assembly$p_val = c(points$assembly$p_val, corr$assembly$p_val[assembly_filt])
		points$disassembly$p_val = c(points$disassembly$p_val, corr$disassembly$p_val[disassembly_filt])
		points$assembly$offset = c(points$assembly$offset, corr$assembly$offset[assembly_filt])
		points$disassembly$offset = c(points$disassembly$offset, corr$disassembly$offset[disassembly_filt])

		points$assembly$stable_lifetime = c(points$assembly$stable_lifetime, res$stable_lifetime[assembly_filt])
		points$disassembly$stable_lifetime = c(points$disassembly$stable_lifetime, res$stable_lifetime[disassembly_filt])
		points$assembly$stable_mean = c(points$assembly$stable_mean, res$stable_mean[assembly_filt])
		points$disassembly$stable_mean = c(points$disassembly$stable_mean, res$stable_mean[disassembly_filt])
		points$assembly$stable_variance = c(points$assembly$stable_variance, res$stable_variance[assembly_filt])
		points$disassembly$stable_variance = c(points$disassembly$stable_variance, res$stable_variance[disassembly_filt])
		
		points$assembly$stable_cv = c(points$assembly$stable_cv, 
									  sqrt(res$stable_variance[assembly_filt])/res$stable_mean[assembly_filt])
		points$disassembly$stable_cv = c(points$disassembly$stable_cv, 
							  			 sqrt(res$stable_variance[disassembly_filt])/res$stable_mean[disassembly_filt])
		
		points$assembly$longevity = c(points$assembly$longevity, res$exp_props$longevity[assembly_filt])
		points$disassembly$longevity = c(points$disassembly$longevity, res$exp_props$longevity[disassembly_filt])

		points$assembly$average_speeds = c(points$assembly$average_speeds, res$exp_props$average_speeds[assembly_filt])
		points$disassembly$average_speeds = c(points$disassembly$average_speeds, res$exp_props$average_speeds[disassembly_filt])
		
		if (any(names(res$exp_props) == 'starting_edge_dist')) {
			points$assembly$edge_dist = c(points$assembly$edge_dist, res$exp_props$starting_edge_dist[assembly_filt])
		}
		
		if (any(names(res$exp_props) == 'ending_edge_dist')) {
			points$disassembly$edge_dist = c(points$disassembly$edge_dist, res$exp_props$ending_edge_dist[disassembly_filt])
		}
		
		points$assembly$lin_num = c(points$assembly$lin_num, which(assembly_filt))
		points$disassembly$lin_num = c(points$disassembly$lin_num, which(disassembly_filt))
		
		points$assembly$exp_dir = c(points$assembly$exp_dir, rep(res$exp_dir, length(which(assembly_filt))))
		points$disassembly$exp_dir = c(points$disassembly$exp_dir, rep(res$exp_dir, length(which(disassembly_filt))))

		points$assembly$exp_num = c(points$assembly$exp_num, rep(i,length(which(assembly_filt))))
		points$disassembly$exp_num = c(points$disassembly$exp_num, rep(i,length(which(disassembly_filt))))
		
		points$joint$birth_dist = c(points$joint$birth_dist, res$exp_props$starting_edge_dist[assembly_filt & disassembly_filt])
		points$joint$death_dist = c(points$joint$death_dist, res$exp_props$ending_edge_dist[assembly_filt & disassembly_filt])
	}
	points$assembly = as.data.frame(points$assembly)
	points$disassembly = as.data.frame(points$disassembly)
	points$joint = as.data.frame(points$joint)
	
	points
}

gather_offset_differences <- function(results_long, results_short, min_R_sq=0.9, max_p_val = 0.05, pos_slope = FALSE) {
	points = list()
	for (i in 1:length(results_long)) {
		long = results_long[[i]]
		short = results_short[[i]]

		assembly_filt = (is.finite(long$assembly$R_sq) & long$assembly$R_sq >= min_R_sq
					   & is.finite(long$assembly$slope) & is.finite(long$assembly$p_val) 
					   & long$assembly$p_val <= max_p_val)
		disassembly_filt = (is.finite(long$disassembly$R_sq) & long$disassembly$R_sq >= min_R_sq 
						  & is.finite(long$disassembly$slope) & is.finite(long$disassembly$p_val) 
						  & long$disassembly$p_val <= max_p_val)
		
		if (pos_slope) {
			assembly_filt = assembly_filt & long$assembly$slope > 0
			disassembly_filt = disassembly_filt & long$disassembly$slope > 0
		}				  
	
		points$assembly$offset = c(points$assembly$offset, long$a$off[assembly_filt] - short$a$off[assembly_filt])
		points$disassembly$offset = c(points$disassembly$offset, long$d$off[disassembly_filt] - short$d$off[disassembly_filt])
		points$assembly$R_sq = c(points$assembly$R_sq, long$a$R_sq[assembly_filt] - short$a$R_sq[assembly_filt])
		points$disassembly$R_sq = c(points$disassembly$R_sq, long$d$R_sq[disassembly_filt] - short$d$R_sq[disassembly_filt])
		
		points$assembly$stable_lifetime = c(points$assembly$stable_lifetime, 
											long$stable_lifetime[assembly_filt] - short$stable_lifetime[assembly_filt])
		points$disassembly$stable_lifetime = c(points$disassembly$stable_lifetime, 
											   long$stable_lifetime[disassembly_filt] - short$stable_lifetime[disassembly_filt])
	}
	points$assembly = as.data.frame(points$assembly)
	points$disassembly = as.data.frame(points$disassembly)
	
	points
}

filter_mixed_area <- function(area, corrected, min_R_sq=0.9, max_p_val = 0.05, pos_slope = FALSE) {
	points = list()
	for (i in 1:length(corrected)) {
		corr = corrected[[i]]
		res = area[[i]]

		assembly_filt = (is.finite(corr$assembly$R_sq) & corr$assembly$R_sq >= min_R_sq
					   & is.finite(corr$assembly$slope) & is.finite(corr$assembly$p_val) 
					   & corr$assembly$p_val <= max_p_val)
		disassembly_filt = (is.finite(corr$disassembly$R_sq) & corr$disassembly$R_sq > min_R_sq 
						  & is.finite(corr$disassembly$slope) & is.finite(corr$disassembly$p_val) 
						  & corr$disassembly$p_val <= max_p_val)
		
		if (pos_slope) {
			assembly_filt = assembly_filt & corr$assembly$slope > 0
			disassembly_filt = disassembly_filt & corr$disassembly$slope > 0
		}
	
		points$assembly$slope = c(points$assembly$slope, res$assembly$slope[assembly_filt])
		points$disassembly$slope = c(points$disassembly$slope, res$disassembly$slope[disassembly_filt])
		points$assembly$R_sq = c(points$assembly$R_sq, res$assembly$R_sq[assembly_filt])
		points$disassembly$R_sq = c(points$disassembly$R_sq, res$disassembly$R_sq[disassembly_filt])
		points$assembly$p_val = c(points$assembly$p_val, res$assembly$p_val[assembly_filt])
		points$disassembly$p_val = c(points$disassembly$p_val, res$disassembly$p_val[disassembly_filt])
		points$assembly$offset = c(points$assembly$offset, res$assembly$offset[assembly_filt])
		points$disassembly$offset = c(points$disassembly$offset, res$disassembly$offset[disassembly_filt])

		points$assembly$stable_lifetime = c(points$assembly$stable_lifetime, res$stable_lifetime[assembly_filt])
		points$disassembly$stable_lifetime = c(points$disassembly$stable_lifetime, res$stable_lifetime[disassembly_filt])
		points$assembly$stable_mean = c(points$assembly$stable_mean, res$stable_mean[assembly_filt])
		points$disassembly$stable_mean = c(points$disassembly$stable_mean, res$stable_mean[disassembly_filt])
		points$assembly$stable_variance = c(points$assembly$stable_variance, res$stable_variance[assembly_filt])
		points$disassembly$stable_variance = c(points$disassembly$stable_variance, res$stable_variance[disassembly_filt])
		
		points$assembly$stable_cv = c(points$assembly$stable_cv, 
									  sqrt(res$stable_variance[assembly_filt])/res$stable_mean[assembly_filt])
		points$disassembly$stable_cv = c(points$disassembly$stable_cv, 
							  			 sqrt(res$stable_variance[disassembly_filt])/res$stable_mean[disassembly_filt])
		
		points$assembly$longevity = c(points$assembly$longevity, res$exp_props$longevity[assembly_filt])
		points$disassembly$longevity = c(points$disassembly$longevity, res$exp_props$longevity[disassembly_filt])

		points$assembly$average_speeds = c(points$assembly$average_speeds, res$exp_props$average_speeds[assembly_filt])
		points$disassembly$average_speeds = c(points$disassembly$average_speeds, res$exp_props$average_speeds[disassembly_filt])
		
		if (any(names(res$exp_props) == 'starting_edge_dist')) {
			points$assembly$edge_dist = c(points$assembly$edge_dist, res$exp_props$starting_edge_dist[assembly_filt])
		}
		
		if (any(names(res$exp_props) == 'ending_edge_dist')) {
			points$disassembly$edge_dist = c(points$disassembly$edge_dist, res$exp_props$ending_edge_dist[disassembly_filt])
		}
		
		points$assembly$lin_num = c(points$assembly$lin_num, which(assembly_filt))
		points$disassembly$lin_num = c(points$disassembly$lin_num, which(disassembly_filt))
		
		points$assembly$exp_dir = c(points$assembly$exp_dir, rep(res$exp_dir, length(which(assembly_filt))))
		points$disassembly$exp_dir = c(points$disassembly$exp_dir, rep(res$exp_dir, length(which(disassembly_filt))))

		points$assembly$exp_num = c(points$assembly$exp_num, rep(i,length(which(assembly_filt))))
		points$disassembly$exp_num = c(points$disassembly$exp_num, rep(i,length(which(disassembly_filt))))
		
		points$joint$birth_dist = c(points$joint$birth_dist, res$exp_props$starting_edge_dist[assembly_filt & disassembly_filt])
		points$joint$death_dist = c(points$joint$death_dist, res$exp_props$ending_edge_dist[assembly_filt & disassembly_filt])
	}
	points$assembly = as.data.frame(points$assembly)
	points$disassembly = as.data.frame(points$disassembly)
	points$joint = as.data.frame(points$joint)
	
	points
}

gather_general_props <- function(results) {
	points = list()
	for (i in 1:length(results)) {
		res = results[[i]]
	
		points$longevity = c(points$longevity, res$exp_props$longevity)
		points$ending_edge = c(points$ending_edge, res$exp_props$ending_edge)
		points$starting_edge = c(points$starting_edge, res$exp_props$starting_edge)
		points$largest_area = c(points$largest_area, res$exp_props$largest_area)
		points$ad_sig = c(points$ad_sig, res$exp_props$ad_sig)
	}	
	points = as.data.frame(points)
	points
}

gather_single_image_props <- function(ind_results) {
	
	ind_data = list();
	
	for (i in 1:length(ind_results)) {
		res = ind_results[[i]]
		filt_by_area = res$Area >= min(res$Area)*3 & res$I_num == 1
		ind_data$Area = c(ind_data$Area, res$Area[filt_by_area]);
		ind_data$ad_sig = c(ind_data$ad_sig, res$Average_adhesion_signal[filt_by_area]);
		ind_data$axial_r = c(ind_data$axial_r, res$MajorAxisLength[filt_by_area]/res$MinorAxisLength[filt_by_area]);
	
		ind_data$cent_dist = c(ind_data$cent_dist, res$Centroid_dist_from_edge[filt_by_area]);
	}
	
	ind_data
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
	if (length(results) == 1) {
		results = results[[1]]
	}
	results
}

load_data_files <- function(dirs, files, headers, inc_exp_names = TRUE, debug = FALSE) {
	results = list()
	
	exp_names = list()
	all_files_present_dirs = c()	
	for (i in 1:length(dirs)) {
		seen_files = 0
		all_files_present = 1;
		for (j in 1:length(files)) {
			this_file = file.path(dirs[i],files[j])
		
			if (! file.exists(this_file)) {
				all_files_present = 0;
			}
		}
		if (all_files_present) {
			all_files_present_dirs = c(all_files_present_dirs, dirs[i]);
			
			regex_range = regexpr("time_series_[[:digit:]]",dirs[[i]])
	        if (regex_range[1] == -1) {
        		exp_names[[length(all_files_present_dirs)]] = dirs[[i]];
	        } else {
    	    	exp_names[[length(all_files_present_dirs)]] = substr(dirs[[i]], regex_range[1], regex_range[1] + attr(regex_range,'match.length'));
        	}
		}
	}
	
	if (debug) {
		print(all_files_present_dirs)
	}
	
	for (i in 1:length(files)) {
		results[[i]] = list()
	}
	
	header_array = array(headers, dim=c(1,length(all_files_present_dirs)))

	for (i in 1:length(all_files_present_dirs)) {
		for (j in 1:length(files)) { 
			results[[j]][[i]] = read.table(file.path(all_files_present_dirs[i],files[j]), header = header_array[j], sep=",")
		}
	}
	if (inc_exp_names) {
		results[[length(files) + 1]] = exp_names
	}
		
	if (length(results) == 1) {
		results = results[[1]]
	}
	results
}

find_col_conf_ints <- function(data, boot.samp = 100) {
	upper = c()
	lower = c()
	mean_vals = c()
	for (j in 1:dim(data)[[2]]) {
		this_col = data[,j];
		this_col = this_col[! is.na(this_col)]
		if (length(this_col) == 0) {
			upper = c(upper, NA)
			lower = c(lower, NA)
			mean_vals = c(mean_vals, NA)
		} else {
			boot_samp = boot(this_col, function(x,y) mean(x[y], na.rm=T), boot.samp)
			conf_int = as.numeric(quantile(boot_samp$t, c(0.025,0.975)))
			lower = c(lower, conf_int[1])
			upper = c(upper, conf_int[2])
			mean_vals = c(mean_vals, mean(this_col, na.rm=T))
		}
	}
	list(upper = upper, lower = lower, mean_vals = mean_vals)
}

trim_args_list <- function(args) {
	for (i in 1:length(args)) {
		if (substr(args[i],0,1) == '-') {
			args[i] = substring(args[i],2)
		}
	}
	args
}

write_assembly_disassembly_periods <- function(result, dir) {
	if (! file.exists(dir)) {
		dir.create(dir,recursive=TRUE)
	}
	
	row_nums = which(! is.na(result$assembly$offset))
	if (! is.null(row_nums)) {
		rows_and_offset = cbind(row_nums, result$assembly$offset[row_nums]);
		write.table(rows_and_offset,file=file.path(dir,'assembly_rows_offsets.csv'), sep=',', row.names=FALSE, col.names=FALSE)
	}
	
	row_nums = which(! is.na(result$disassembly$offset))
	if (! is.null(row_nums)) {
		rows_and_offset = cbind(row_nums, result$disassembly$offset[row_nums]);
		write.table(rows_and_offset,file=file.path(dir,'disassembly_rows_offsets.csv'), sep=',', row.names=FALSE, col.names=FALSE)
	}
		
}

write_high_r_rows <- function(result, dir, file=c('assembly_R_sq.csv','disassembly_R_sq.csv','neg_slope_R_sq.csv'), min_R_sq = 0.9) {
	if (! file.exists(dir)) {
		dir.create(dir,recursive=TRUE)
	}
	
	row_nums = which(is.finite(result$assembly$R_sq) & result$assembly$R_sq > min_R_sq)
	if (! is.null(row_nums)) {
		write.table(row_nums,file=file.path(dir,file[1]), row.names=FALSE, col.names=FALSE)
	}
	
	row_nums = which(is.finite(result$assembly$R_sq) & result$disassembly$R_sq > min_R_sq)
	if (! is.null(row_nums)) {
		write.table(row_nums,file=file.path(dir,file[2]), row.names=FALSE, col.names=FALSE)
	}
	
	row_nums = which(is.finite(result$assembly$R_sq) & result$assembly$R_sq > min_R_sq & result$assemblyassemblyassembly$slope < 0)
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

find_death_time <- function(alive_dead) {
	#alive_dead should be a sequence of T/F values, where T is alive and F is dead
	stopifnot(is.logical(alive_dead))
	
	#set the death time to NA, to indicate the default of no death observed
	death_time = NA
	
	#First deal with the posibility that the last entry is T, so the adhesion 
	#lives to the end of the experiment, thus, we don't see a death
	if (! alive_dead[length(alive_dead)]) {
		#If the first entry is T (alive), then we need to find the first dead
		#entry to find the death time
		if (alive_dead[1]) {
			dead_points = which(!alive_dead)
			if (any(dead_points)) {
				death_time = dead_points[1]
			}	
		} else {
			#We know the first entry is F (dead) and that the adhesion dies, now 
			#we search for the spot where the adhesion transistions from live 
			#to dead
			for (i in 2:length(alive_dead)) {
				if (alive_dead[i - 1] & !alive_dead[i]) {
					death_time = i
				}
			}
		}
	}
	death_time
}

find_birth_time <- function(alive_dead) {
	#alive_dead should be a sequence of T/F values, where T is alive and F is dead
	stopifnot(is.logical(alive_dead))
	
	#set the birth time to NA, to indicate the default of no birth observed
	birth_time = NA
	
	#all we need to search for is the first T (alive) entry, assuming the first 
	#entry is F (dead) otherwise, we don't know when birth occured
	if (! alive_dead[1]) {
		birth_time = which(alive_dead)[1]
	}
	birth_time
}

########################################
#Tests
########################################

#Birth/Death Rates
stopifnot(is.na(find_death_time(c(F,T,T,T))))
stopifnot(find_death_time(c(T,T,F,F)) == 3)
stopifnot(find_death_time(c(F,T,F,F)) == 3)
stopifnot(find_death_time(c(F,T,T,F)) == 4)

stopifnot(is.na(find_birth_time(c(T,T,T,F))))
stopifnot(find_birth_time(c(F,T,F,F)) == 2)


################################################################################
# Main Program
################################################################################

args <- commandArgs(TRUE)

if (length(args) != 0) {
	args <- trim_args_list(args)
    
    #each of the outputs of the following commands are saved to temp to avoid
    #writing the entire results contents to STDOUT, very useful when debugging
    #runs from the command line
	average_model = gather_bilinear_models_from_dirs(args,
		data_file='Average_adhesion_signal.csv',
		results.file=file.path('..','models','intensity_model.Rdata'))
	
	write_assembly_disassembly_periods(average_model[[1]],file.path(args[[1]],'..'))
		
	temp = gather_bilinear_models_from_dirs(args, 
		data_file='Background_corrected_signal.csv', 
		results.file=file.path('..','models','corrected_intensity_model.Rdata'))

	temp = gather_bilinear_models_from_dirs(args, 
		data_file='Shrunk_corrected_signal.csv', 
		results.file=file.path('..','models','shrunk_intensity_model.Rdata'))

	temp = gather_bilinear_models_from_dirs(args, 
		data_file='Area.csv', 
		results.file=file.path('..','models','area_model.Rdata'))
	
    temp = gather_bilinear_models_from_dirs(args, 
		data_file='Box_intensity.csv', 
		results.file=file.path('..','models','box_model.Rdata'))
}