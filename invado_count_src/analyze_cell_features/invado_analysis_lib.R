################################################################################
# invado_analysis_lib.R: the R functions associated with processing the
# invaodopodia data
################################################################################

library(Hmisc);

get_non_nan_data <- function(files) {
    data = c();
    
    for (this_file in files) {
        temp = read.csv(this_file,header=F)
        temp = temp[!is.na(temp)];
        data = c(data,temp);

        if (any(temp < -500)) {
            print(this_file)
        }
    }
    
    return(data)
}

colMedians <- function(this_mat,na.rm=T) {
    medis = c()
    for (i in 1:dim(this_mat)[2]) {
        medis = c(medis,median(this_mat[,i],na.rm=na.rm))
    }
    return(medis)
}

colConfUpper <- function(this_mat) {
    upper = c()
    for (i in 1:dim(this_mat)[2]) {
        if (all(is.na(this_mat[,i]))) {
            upper = c(upper,NA);
        } else {
            temp = t.test(this_mat[,i],conf.level=0.95)
            if (is.nan(temp$conf.int[2])) {
                upper = c(upper,0)
            } else {
                upper = c(upper,temp$conf.int[2])
            }
        }
    }
    return(upper)
}

colConfLower <- function(this_mat) {
    lower = c()
    for (i in 1:dim(this_mat)[2]) {
        if (all(is.na(this_mat[,i]))) {
            lower = c(lower,NA);
        } else {
            temp = t.test(this_mat[,i],conf.level=0.95)
            if (is.nan(temp$conf.int[1])) {
                lower = c(lower,0)
            } else {
                lower = c(lower,temp$conf.int[1])
            }
        }
    }
    return(lower)
}

################################################################################
# Plotting
################################################################################

plot_single_exp_mat <- function(time,exp_mat,upper_conf=NA,lower_conf=NA,col=NA,...) {
    
    if (! is.na(col)) {
        col_light = rgb(t(col2rgb(col))/255,alpha=0.5);
    } else {
        col_light = NA;
    }
    
    for (i in 1:dim(exp_mat)[1]) {
        if (i == 1) {
            plot(time,exp_mat[i,],col=col_light,...)
        } else {
            lines(time,exp_mat[i,],col=col_light,...)
        }
    }
    
    lines(time,colMeans(exp_mat),col=col,lwd=3,...)
    
    if (! is.na(upper_conf[1]) && ! is.na(lower_conf[1])) {
        polygon(c(time,rev(time)), c(lower_conf,rev(upper_conf)),col=col_light,border=NA)
    }
}

plot_multiple_exp_mat <- function(time,exps,exp_confs,cols=NA,exp_names=NA,...) {
    
    #colors
    if (is.na(cols[1])) {
        cols = rainbow(length(names(exps)));
    }

    cols_light = list()
    for (exp_type in names(exps)) {
        cols_light[[exp_type]] = rgb(t(col2rgb(cols[[exp_type]]))/255,alpha=0.5)
    }

    #exps to include
    if (is.na(exp_names[1])) {
        exp_names = names(exps);
    }
    
    first = 1;
    for (exp_type in exp_names) {
        print(exp_type)
        if (first) {            
            plot(time,colMeans(exps[[exp_type]]),col=cols[[exp_type]],lwd=3,typ='l',...)
            first = 0;
        } else {
            lines(time,colMeans(exps[[exp_type]]),col=cols[[exp_type]],lwd=3,...)
        }
        polygon(c(time,rev(time)), c(exp_confs[[exp_type]]$upper,rev(exp_confs[[exp_type]]$lower)),
            col=cols_light[[exp_type]],border=cols_light[[exp_type]])
    }
}

################################################################################
# Main Program
################################################################################
