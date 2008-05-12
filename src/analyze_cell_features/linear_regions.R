plot_ad_seq <- function (plot_params,this_ad_seq,num,dir='../plots/linear_regions/') {
	
	this_ad_seq = this_ad_seq[!(is.nan(this_ad_seq))];
	if (! file.exists(dir)) {
		dir.create(dir)
	}
	
	pdf(paste(dir,num,'.pdf',sep=''));
	plot(1:dim(this_ad_seq)[[2]],this_ad_seq,xlab='Time (minutes)',ylab='Average Paxillin Signal (AU)')
	
	early_x = c(0,plot_params$early_offset);
	early_y = c(early_x[[1]]*plot_params$early_slope + plot_params$early_inter,
				early_x[[2]]*plot_params$early_slope + plot_params$early_inter);
	lines(early_x,early_y,col='red',lwd=2)
	
	late_x = c(dim(this_ad_seq)[[2]] - plot_params$late_offset, dim(this_ad_seq)[[2]]);
	late_y = c(late_x[[1]]*plot_params$late_slope + plot_params$late_inter,
			   late_x[[2]]*plot_params$late_slope + plot_params$late_inter);
	lines(late_x, late_y, col='blue',lwd=2)
	
	#legend(late_x[[1]] + 1, late_y[[1]], paste(sprintf('%.1f',plot_params$late_R_sq),plot_params$late_slope,sep='\n'))
	#legend(early_x[[1]] + 1, early_y[[1]], paste(sprintf('%.1f',plot_params$early_R_sq),plot_params$early_slope,sep='\n'))
	dev.off();
}

gather_linear_regions <- function (dir, min_length=5) {

	ad_sig <- read.table(paste(dir,'Average_adhesion_signal.csv',sep=''),header = FALSE, sep  = ',');
	ad_props <- read.table(paste(dir,'../single_lin.csv',sep=''), header = TRUE, sep=',');

	m_results <- list(early_offset = array(min_length, dim = c(dim(ad_sig)[[1]])),
					  early_R_sq = array(0, dim = c(dim(ad_sig)[[1]])), 
#					  early_inter = array(NaN, dim = c(dim(ad_sig)[[1]])), 
#					  early_slope = array(NaN, dim = c(dim(ad_sig)[[1]])),
					  late_offset = array(min_length, dim = c(dim(ad_sig)[[1]])),
					  late_R_sq = array(0, dim = c(dim(ad_sig)[[1]]))
#					  late_inter = array(NaN, dim = c(dim(ad_sig)[[1]])), 
#					  late_slope = array(NaN, dim = c(dim(ad_sig)[[1]]))
					 );
					 
	for (i in 1:dim(ad_sig)[[1]]) {
		if ((i %% 1000) == 0) {
			#print(i)
		}
		if (ad_props$longevity[[i]] < min_length * 2) {
			next
		}
	
		temp = ad_sig[i,];
		temp = temp[!(is.nan(temp))];
		this_ad_sig = data.frame(y = t(temp[1,]), x = 1:dim(temp)[[2]]);
		names(this_ad_sig) <- c('y','x');
	
		cur_offset = m_results$early_offset[[i]]
		for (j in cur_offset:dim(this_ad_sig)[[1]]) {
			this_subset = this_ad_sig[1:cur_offset,]
		
			model <- lm(y ~ x, data = this_subset)
		
			summary <- summary(model);
			
			if (summary$r.squared > m_results$early_R_sq[[i]]) {
				m_results$early_R_sq[[i]] = summary$r.squared;
				m_results$early_length[[i]] = length(this_subset);
				m_results$early_offset[[i]] = cur_offset;
				m_results$early_inter[[i]] = coef(model)[[1]];
				m_results$early_slope[[i]] = coef(model)[[2]];
			}
		
			cur_offset = cur_offset + 1
		}
	
		cur_offset = m_results$late_offset[[i]]
		for (j in cur_offset:dim(this_ad_sig)[[1]]) {
			this_subset = this_ad_sig[(dim(this_ad_sig)[[1]]-cur_offset):dim(this_ad_sig)[[1]],]
		
			model <- lm(y ~ x, data = this_subset)
				
			summary <- summary(model);
			
			if (summary$r.squared > m_results$late_R_sq[[i]]) {				m_results$late_R_sq[[i]] = summary$r.squared;
				m_results$late_length[[i]] = length(this_subset);
				m_results$late_offset[[i]] = cur_offset;
				m_results$late_inter[[i]] = coef(model)[[1]];
				m_results$late_slope[[i]] = coef(model)[[2]];
			}
		
			cur_offset = cur_offset + 1
		}
	
	}

	for (i in 1:dim(ad_sig)[[1]]) {

		if (ad_props$longevity[[i]] < min_length * 2) {
			next
		}
		if (m_results$early_R_sq[[i]] < 0.9 || m_results$late_R_sq[[i]] < 0.9) {
			next
		}
#		if (m_results$early_slope[[i]] < 0 || m_results$late_slope[[i]] > 0) {
#			next
#		}
	
		#plot_ad_seq(m_results[i,],ad_sig[i,],i,dir=paste(dir,'../plots/linear_regions/',sep=''))
		count = count + 1;
	}

	m_results
}


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

all_data <- list();

#for (i in 1:length(dirs)) {
    for (i in 2:2) {
        dir = dirs[[i]];

        all_data$five_results[[i]] = gather_linear_regions(dir);
        all_data$ten_results[[i]] = gather_linear_regions(dir,10);
        five_results <- all_data$five_results[[i]];
        ten_results <- all_data$ten_results[[i]];

        five_early <- array(0,c(10,1));
        five_late <- array(0,c(10,1));
        ten_late <- array(0,c(10,1));
        ten_early <- array(0,c(10,1));
        five_early_error <- array(0,c(10,1));
        five_late_error <- array(0,c(10,1));
        ten_late_error <- array(0,c(10,1));
        ten_early_error <- array(0,c(10,1));

        for (i in 1:10) {
            five_set_early = five_results$early_R_sq > ((i-1)/10);		five_set_late = five_results$late_R_sq > ((i-1)/10);
            ten_set_early = ten_results$early_R_sq > ((i-1)/10);
            ten_set_late = ten_results$late_R_sq > ((i-1)/10);

            five_early[[i]] = mean(five_results$early_slope[five_set_early]);
            five_early_error[[i]] = sqrt(var(five_results$early_slope[five_set_early]/length(five_results$early_slope[five_set_early])));

            five_late[[i]] = mean(five_results$late_slope[five_set_late]);
            five_late_error[[i]] = sqrt(var(five_results$late_slope[five_set_late]/length(five_results$late_slope[five_set_late])));

            ten_early[[i]] = mean(ten_results$early_slope[ten_set_early]);
            ten_early_error[[i]] = sqrt(var(ten_results$early_slope[ten_set_early]/length(ten_results$early_slope[ten_set_early])));

            ten_late[[i]] = mean(ten_results$late_slope[ten_set_late]);
            ten_late_error[[i]] = sqrt(var(ten_results$late_slope[ten_set_late]/length(ten_results$late_slope[ten_set_late])));
        }

        pdf(paste(dir,'/../plots/linear_regions.pdf',sep=''),width=8.5,height=8.5,pointsize=14);
        par(mfrow=c(2,2),bty='n')

#Plot 1 - Slope versus R squared (early)
        all_x = c(five_results$early_slope[five_results$early_R_sq != 0],
                  ten_results$early_slope[ten_results$early_R_sq != 0]);
        all_y = c(five_results$early_R_sq[five_results$early_R_sq != 0],
                  ten_results$early_R_sq[ten_results$early_R_sq != 0]);

        plot(five_results$early_slope[five_results$early_R_sq != 0],
             five_results$early_R_sq[five_results$early_R_sq != 0],
             xlab='Early Slope', ylab='Early R Squared',
             xlim=c(min(all_x),max(all_x)),ylim=c(min(all_y),max(all_y)))
        points(ten_results$early_slope[ten_results$early_R_sq != 0],
               ten_results$early_R_sq[ten_results$early_R_sq != 0],
               col='red');

#Plot 2 - R squared versus Slope (late)
        all_x = c(five_results$late_slope[five_results$early_R_sq != 0],ten_results$late_slope[ten_results$early_R_sq != 0]);
        all_y = c(five_results$late_R_sq[five_results$early_R_sq != 0],ten_results$late_R_sq[ten_results$early_R_sq != 0]);

        plot(five_results$late_slope[five_results$early_R_sq != 0],	 five_results$late_R_sq[five_results$early_R_sq != 0],
                xlab='Late Slope', ylab='Late R Squared',
                xlim=c(min(all_x),max(all_x)),ylim=c(min(all_y),max(all_y)))
            points(ten_results$late_slope[ten_results$early_R_sq != 0],
                    ten_results$late_R_sq[ten_results$early_R_sq != 0],
                    col='red');

#Plot 3
        all_x = seq(0,0.9,by=0.1);
        all_y = c(five_early,ten_early)
            errbar(seq(0,0.9,by=0.1),five_early,
                    five_early - five_early_error,five_early + five_early_error, 
                    xlim=c(min(all_x),max(all_x)),ylim=c(min(all_y),max(all_y)),
                    xlab = 'R squared cutoff', ylab='Accumulation Rate (/min)')
            points(seq(0,0.9,by=0.1),ten_early,col='red')
            errbar(seq(0,0.9,by=0.1),ten_early, ten_early - ten_early_error, ten_early + ten_early_error, add=T, col='red')

#Plot 4
            all_x = seq(0,0.9,by=0.1);
        all_y = c(five_late,ten_late)
            errbar(seq(0,0.9,by=0.1),five_late,
                    five_late - five_late_error,five_late + five_late_error, 
                    xlim=c(min(all_x),max(all_x)),ylim=c(min(all_y),max(all_y)),
                    xlab = 'R squared cutoff', ylab='Degradation Rate (/min)')
            points(seq(0,0.9,by=0.1),ten_late,col='red')
            errbar(seq(0,0.9,by=0.1),ten_late, ten_late - ten_late_error, ten_late + ten_late_error, add=T, col='red')
            dev.off()

    }
