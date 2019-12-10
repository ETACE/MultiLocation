
library(parallel)
library(RSQLite)
library(mgcv)



gamplot <- function (x, residuals = FALSE, rug = TRUE, se = TRUE, pages = 0, 
    select = NULL, scale = -1, n = 100, n2 = 40, pers = FALSE, 
    theta = 30, phi = 30, jit = FALSE, xlab = NULL, ylab = NULL, 
    main = NULL, ylim = NULL, xlim = NULL, too.far = 0.1, all.terms = FALSE, 
    shade = FALSE, shade.col = "gray80", shift = 0, trans = I, 
    seWithMean = FALSE, ...) 
{
#par(ask=FALSE)
    sub.edf <- function(lab, edf) {
        pos <- regexpr(":", lab)[1]
        if (pos < 0) {
            pos <- nchar(lab) - 1
            lab <- paste(substr(lab, start = 1, stop = pos), 
                ",", round(edf, digits = 2), ")", sep = "")
        }
        else {
            lab1 <- substr(lab, start = 1, stop = pos - 2)
            lab2 <- substr(lab, start = pos - 1, stop = nchar(lab))
            lab <- paste(lab1, ",", round(edf, digits = 2), lab2, 
                sep = "")
        }
        lab
    }

    w.resid <- NULL
    if (length(residuals) > 1) {
        if (length(residuals) == length(x$residuals)) 
            w.resid <- residuals
        else warning("residuals argument to plot.gam is wrong length: ignored")
        partial.resids <- TRUE
    }
    else partial.resids <- residuals
    m <- length(x$smooth)
    order <- attr(x$pterms, "order")
    if (all.terms) 
        n.para <- sum(order == 1)
    else n.para <- 0
    if (m + n.para == 0) 
        stop("No terms to plot - nothing for plot.gam() to do.")
    if (se) {
        if (is.numeric(se)) 
            se2.mult <- se1.mult <- se
        else {
            se1.mult <- 2
            se2.mult <- 1
        }
        if (se1.mult < 0) 
            se1.mult <- 0
        if (se2.mult < 0) 
            se2.mult <- 0
    }
    else se1.mult <- se2.mult <- 1
    if (se && x$Vp[1, 1] <= 0) {
        se <- FALSE
        warning("No variance estimates available")
    }
    n.plots <- m + n.para
    if (pages > n.plots) 
        pages <- n.plots
    if (pages < 0) 
        pages <- 0
    if (pages != 0) {
        ppp <- n.plots%/%pages
        if (n.plots%%pages != 0) {
            ppp <- ppp + 1
            while (ppp * (pages - 1) >= n.plots) pages <- pages - 
                1
        }
        c <- trunc(sqrt(ppp))
        if (c < 1) 
            c <- 1
        r <- ppp%/%c
        if (r < 1) 
            r <- 1
        while (r * c < ppp) r <- r + 1
        while (r * c - ppp > c && r > 1) r <- r - 1
        while (r * c - ppp > r && c > 1) c <- c - 1
        #oldpar <- par(mfrow = c(r, c),ask=FALSE)
    }
    else {
        ppp <- 1
        #oldpar <- par()
    }
    if (partial.resids) {
        fv.terms <- predict(x, type = "terms")
        if (is.null(w.resid)) 
            w.resid <- x$residuals * sqrt(x$weights)
    }
    pd <- list()
    i <- 1
    if (m > 0) 
        for (i in 1:m) {
            if (x$smooth[[i]]$dim == 1) {
                raw <- x$model[x$smooth[[i]]$term]
                xx <- seq(min(raw), max(raw), length = n)
                if (x$smooth[[i]]$by != "NA") {
                  by <- rep(1, n)
                  dat <- data.frame(x = xx, by = by)
                  names(dat) <- c(x$smooth[[i]]$term, x$smooth[[i]]$by)
                }
                else {
                  dat <- data.frame(x = xx)
                  names(dat) <- x$smooth[[i]]$term
                }
                X <- PredictMat(x$smooth[[i]], dat)
                first <- x$smooth[[i]]$first.para
                last <- x$smooth[[i]]$last.para
                p <- x$coefficients[first:last]
                offset <- attr(X, "offset")
                if (is.null(offset)) 
                  fit <- X %*% p
                else fit <- X %*% p + offset
                if (se) {
                  if (seWithMean && inherits(attr(x$smooth[[i]], 
                    "qrc"), "qr")) {
                    X1 <- matrix(x$cmX, nrow(X), ncol(x$Vp), 
                      byrow = TRUE)
                    meanL1 <- x$smooth[[i]]$meanL1
                    if (!is.null(meanL1)) 
                      X1 <- X1/meanL1
                    X1[, first:last] <- X
                    se.fit <- sqrt(rowSums((X1 %*% x$Vp) * X1))
                  }
                  else se.fit <- sqrt(rowSums((X %*% x$Vp[first:last, 
                    first:last]) * X))
                }
                edf <- sum(x$edf[first:last])
                xterm <- x$smooth[[i]]$term
                if (is.null(xlab)) 
                  xlabel <- xterm
                else xlabel <- xlab
                if (is.null(ylab)) 
                  ylabel <- sub.edf(x$smooth[[i]]$label, edf)
                else ylabel <- ylab
                pd.item <- list(fit = fit, dim = 1, x = xx, ylab = ylabel, 
                  xlab = xlabel, raw = raw[[1]])
                if (partial.resids) {
                  pd.item$p.resid <- fv.terms[, length(order) + 
                    i] + w.resid
                }
                if (se) 
                  pd.item$se = se.fit * se1.mult
                pd[[i]] <- pd.item
                rm(pd.item)
            }
            else if (x$smooth[[i]]$dim == 2) {
                xterm <- x$smooth[[i]]$term[1]
                if (is.null(xlab)) 
                  xlabel <- xterm
                else xlabel <- xlab
                yterm <- x$smooth[[i]]$term[2]
                if (is.null(ylab)) 
                  ylabel <- yterm
                else ylabel <- ylab
                raw <- data.frame(x = as.numeric(x$model[xterm][[1]]), 
                  y = as.numeric(x$model[yterm][[1]]))
                n2 <- max(10, n2)
                xm <- seq(min(raw$x), max(raw$x), length = n2)
                ym <- seq(min(raw$y), max(raw$y), length = n2)
                xx <- rep(xm, n2)
                yy <- rep(ym, rep(n2, n2))
                if (too.far > 0) 
                  exclude <- exclude.too.far(xx, yy, raw$x, raw$y, 
                    dist = too.far)
                else exclude <- rep(FALSE, n2 * n2)
                if (x$smooth[[i]]$by != "NA") {
                  by <- rep(1, n2^2)
                  dat <- data.frame(x = xx, y = yy, by = by)
                  names(dat) <- c(xterm, yterm, x$smooth[[i]]$by)
                }
                else {
                  dat <- data.frame(x = xx, y = yy)
                  names(dat) <- c(xterm, yterm)
                }
                X <- PredictMat(x$smooth[[i]], dat)
                first <- x$smooth[[i]]$first.para
                last <- x$smooth[[i]]$last.para
                p <- x$coefficients[first:last]
                offset <- attr(X, "offset")
                if (is.null(offset)) 
                  fit <- X %*% p
                else fit <- X %*% p + offset
                fit[exclude] <- NA
                if (se) {
                  if (seWithMean && inherits(attr(x$smooth[[i]], 
                    "qrc"), "qr")) {
                    X1 <- matrix(x$cmX, nrow(X), ncol(x$Vp), 
                      byrow = TRUE)
                    meanL1 <- x$smooth[[i]]$meanL1
                    if (!is.null(meanL1)) 
                      X1 <- X1/meanL1
                    X1[, first:last] <- X
                    se.fit <- sqrt(rowSums((X1 %*% x$Vp) * X1))
                  }
                  else se.fit <- sqrt(rowSums((X %*% x$Vp[first:last, 
                    first:last]) * X))
                  se.fit[exclude] <- NA
                }
                edf <- sum(x$edf[first:last])
                if (is.null(main)) {
                  title <- sub.edf(x$smooth[[i]]$label, edf)
                }
                else title <- main
                pd.item <- list(fit = fit, dim = 2, xm = xm, 
                  ym = ym, ylab = ylabel, xlab = xlabel, title = title, 
                  raw = raw)
                if (is.null(ylim)) 
                  pd.item$ylim <- range(ym)
                else pd.item$ylim <- ylim
                if (is.null(xlim)) 
                  pd.item$xlim <- range(xm)
                else pd.item$xlim <- xlim
                if (se) 
                  pd.item$se = se.fit * se2.mult
                pd[[i]] <- pd.item
                rm(pd.item)
            }
            else {
                pd[[i]] <- list(dim = x$smooth[[i]]$dim)
            }
        } # Ende for-Schleife i in 1:m
return(pd)
}

present_value = function(data){

pv = 0

for(t in 1:length(data)){

pv = pv + discount_factor^t * data[t]

}

return(pv)

}




equivalent_profit_per_year = function(data){

pv = 0

for(t in 1:length(data)){

pv = pv + discount_factor^t * data[t]




}

pv = pv * (1-discount_factor)/(1-discount_factor^length(data)) 

return(pv)

}


stillRunning = TRUE

readkey <- function()
{
    cat("[press [enter] to continue]")
    number <- scan(n=1)
}
isEmpty= function(x){

return(length(x)==0)


}



detect_modes_innovators_max= function(){

    
    
  modes = c()
  
  
  
    
  for(p in 1:length(par))  {
  
  
  DATA <- c()

  for(r in 1:runs){


		
		
		eval(parse(text=paste("path ='",path_part1,par[p],path_part2,"'",sep="")))
	
		  print(path)
		 con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))



		firm = dbReadTable(con,"Firm")
  

		inno<-dbGetQuery(con,paste("SELECT currentInnoProb FROM Firm",sep=""))
		inno<-as.numeric(unlist(inno))

		inno = firm$currentInnoProb  
		inno[inno>0]=1

		da = firm$equilQuantity
		
		da[da>0]<-1
		 
		firm = data.frame(firm , isactive = da) 
		 
		datai <-  eval(parse(text=paste("data.frame( id = r, t = firm$X_ITERATION_NO, par=as.character(par[p]),firmID= firm$firmID, numLocations = firm$numLocations, inno= inno, mode = 1)")))
		DATA <-  rbind(DATA,datai)


		dbDisconnect(con)

		


  }


  
     mode1 = c()
     mode2= c()
     
   
      data_inno = DATA[DATA$inno == 1, ]
      data_imi = DATA[DATA$inno == 0, ]
     
     for(r in 1:runs){
     
     
    
   
     
     
	    if(max(as.numeric(data_inno[data_inno$id==r,]$numLocations))==5){
	    
		    mode1=c(mode1,r)
	    	    
	    }else{

		 mode2=c(mode2,r)
		
	    }
     
	   print(paste("Find max in run ",r," done",sep=""))
     }
     
     
     print(mode1)
     print(mode2)
     
     
     
      temp = list(list(mode1,mode2))
      modes = c(modes, temp )
     
     }
     
     
     
     
     
      pdf(paste(path_figs_base,"/modes.pdf",sep=""))
      
      num_mode1 = c()
      
      for(p in 1:length(par)){
      
	  num_mode1 = c(num_mode1, length(modes[[p]][[1]]))
      
      }
     
      plot(num_mode1,type="l")
      dev.off()
     
    return(modes)
   
   
   }
      



 detect_modes_imiators_min= function(){

    
    
  modes = c()
  
  
  
    
  for(p in 1:length(par))  {
  
  
	  RUNS = c()
	  
	  
	   DATA <- c()
  
  
  for(r in 1:runs){
  
  
  
	  
  
		
		eval(parse(text=paste("path ='",path_part1,par[p],path_part2,"'",sep="")))
		print(paste(path,"/run_",r,sep=""))
		notDone=FALSE
		
		if(stillRunning){
		
		  if(!file.exists(paste(path,"/run_",r,"/output.txt",sep=""))){
		  
		  notDone = TRUE
		  print("Not done yet")
		  
		  }else{
		  
		  
		  
		  da=read.csv(paste(path,"/run_",r,"/output.txt",sep=""))
		  
		 
		  da=da[,1]
		  da = paste(da,collapse=" ")
		  
		  print(grepl("Batch Done",da))
		  if(!grepl("Batch Done",da)){
		  
		  notDone = TRUE
		  print("Not done yet")
		  
		  }

		}
		
		
		}
		
		if(!notDone){
  
  
  
		  RUNS = c(RUNS, r)
  
  
		



		
		
		eval(parse(text=paste("path ='",path_part1,par[p],path_part2,"'",sep="")))
	
		  print(path)
		 con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))



		firm = dbReadTable(con,"Firm")
  

		inno<-dbGetQuery(con,paste("SELECT currentInnoProb FROM Firm",sep=""))
		inno<-as.numeric(unlist(inno))

		inno = firm$currentInnoProb  
		inno[inno>0]=1

		da = firm$equilQuantity
		
		da[da>0]<-1
		 
		firm = data.frame(firm , isactive = da) 
		 
		datai <-  eval(parse(text=paste("data.frame( id = r, t = firm$X_ITERATION_NO, par=as.character(par[p]),firmID= firm$firmID, numLocations = firm$numLocations, inno= inno, mode = 1)")))
		DATA <-  rbind(DATA,datai)


		dbDisconnect(con)

		


  }

  }
  
  
  
     mode1 = c()
     mode2= c()
     
   
      data_inno = DATA[DATA$inno == 1, ]
      data_imi = DATA[DATA$inno == 0, ]
      
     
     for(r in RUNS){
     
     
     
     
	    imi_exist = TRUE
	    
	    
	    data_im_r = data_imi[data_imi$id==r,]
	    
    
	    for(t in 1:itno){
	    
		  if(isEmpty(data_im_r[data_im_r$t==t,]$mode)){
		  imi_exist = FALSE
		   print("Mode1")
		  break
		  
		 
		  
		  }

	    }
     
     
	    if(imi_exist){
	    
		    mode2=c(mode2,r)
	    	    
	    }else{

		 mode1=c(mode1,r)
		
	    }
     
	   print(paste("Find max in run ",r," done",sep=""))
     }
     
     
     print(mode1)
     print(mode2)
     
     
     
      temp = list(list(mode1,mode2))
      modes = c(modes, temp )
     
     }
     
     
     
     
     
      pdf(paste(path_figs_base,"/modes.pdf",sep=""))
      
      num_mode1 = c()
      
      for(p in 1:length(par)){
      
	  num_mode1 = c(num_mode1, length(modes[[p]][[1]])/(length(modes[[p]][[1]])+length(modes[[p]][[2]])))
      
      }
     
      plot(num_mode1,type="l",xlab = "Entry Costs", ylab="rel. Frequency", xaxt="n", lwd=3)
      axis(1, at=1:length(num_mode1), labels=par)
      
      
      dev.off()
     
    return(modes)
   
   
   }
      
  entry_exit_rates = function()  {
  
  
	  av_inno_en = list()
	  av_imi_en = list()
	  av_strat_en = list()
	  
	  av_inno_ex = list()
	  av_imi_ex = list()
	  av_strat_ex = list()
	  
	  av_inno_sw = list()
	  av_imi_sw = list()
	  av_strat_sw = list()
  
  
	
	for(p in par){
	
	
		inno_en=c()
		imi_en=c()
		strat_en=c()
		
		inno_ex=c()
		imi_ex=c()
		strat_ex=c()
		
		inno_sw=c()
		imi_sw=c()
		strat_sw=c()


		eval(parse(text=paste("path ='",path_part1,p,path_part2,"'",sep="")))
		

		for(r in 1:runs){
		print(paste(path,"/run_",r,sep=""))
		
		
		notDone=FALSE
		
		if(stillRunning){
		
		  if(!file.exists(paste(path,"/run_",r,"/output.txt",sep=""))){
		  
		  notDone = TRUE
		  print("Not done yet")
		  
		  }else{
		  
		  
		  
		  da=read.csv(paste(path,"/run_",r,"/output.txt",sep=""))
		  
		 
		  da=da[,1]
		  da = paste(da,collapse=" ")
		  
		  print(grepl("Batch Done",da))
		  if(!grepl("Batch Done",da)){
		  
		  notDone = TRUE
		  print("Not done yet")
		  
		  }

		}
		
		
		}
		
		if(!notDone){
		
		
		
		
		con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
		
		
		
		#con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
		
		
		firms = dbReadTable(con,"Firm")
		exits = dbReadTable(con,"Exit")
		entry = dbReadTable(con,"Entry")
		
		
		
		
		id_max = max(firms$firmID)
		
		
		
		for(i in 1:id_max){
		
		
		
		
		    firm  = firms[firms$firmID==i,]
		    
		    exitsFirm  = exits[exits$firmID==i  & exits$exit == 1,]
		    entriesFirm = entry[entry$firmID==i & entry$entry == 1 & entry$marketEntry==0,]
		    
		  
		    
		    if(length(entriesFirm$entry)== 0 && length(exitsFirm$exit) == 0){
		    
		    
		    switch = 0
		    ent = 0
		    ex = 0
		    
		    
		    
		    }else if(length(entriesFirm$entry)> 0 && length(exitsFirm$exit) == 0){
		    
		    
		    switch = 0
		    ex = 0 
		    ent = sum (entriesFirm$entry) 
		    
		    
		    }else if(length(entriesFirm$entry)== 0 && length(exitsFirm$exit) > 0){
		    
		    
		    switch = 0
		    ex = sum (exitsFirm$exit)
		    ent = 0 
		    
		    
		    }else{
		    
		    
		    switch = 0
		    
		    for(e in 1: length(entriesFirm$entry)){
		    
		 
			  for(f in 1: length(exitsFirm$exit)){
			  
			   
		    
				if(exitsFirm$X_ITERATION_NO[f]==entriesFirm$X_ITERATION_NO[e]){
				
				    switch = switch + 1
				    
				    entriesFirm$entry[e]=0
				    exitsFirm$exit[f] = 0
				
				}

		      }
		    
		      }
		      
		      
		       ent = sum(entriesFirm$entry) 
			ex = sum (exitsFirm$exit)    
			
		      
		    
		    }
		    
		    ent = ent/ length(firm$firmID)
		    ex = ex   / length(firm$firmID)
		    switch = switch / length(firm$firmID)
		    
		    
		    #eval(parse(text=paste("data = ",FUN,"(firm$",variable,")",sep="")))
		   
		
		    
		    if(firm$firmStratgy[1] == 1){
		    
		    strat_ex = c(strat_ex,ex)
		    strat_en = c(strat_en,ent)
		    strat_sw = c(strat_sw,switch)
		    
		    }else{
		    
		    if(firm$innovator[1] == 1){
		    
		     
		    inno_ex = c(inno_ex,ex)
		    inno_en = c(inno_en,ent)
		    inno_sw = c(inno_sw,switch)
		    
		    }else{
		    
		     imi_ex = c(imi_ex,ex)
		    imi_en = c(imi_en,ent)
		    imi_sw = c(imi_sw,switch)
		    
		    }
		 
		
		
		
		
		}
		
		
		}
		
		}
		
		}
		
	
		
		
		av_inno_en= c(av_inno_en, list(inno_en))
		av_imi_en= c(av_imi_en, list(imi_en))
		av_strat_en= c(av_strat_en, list(strat_en))
		
		av_inno_ex= c(av_inno_ex, list(inno_ex))
		av_imi_ex= c(av_imi_ex, list(imi_ex))
		av_strat_ex= c(av_strat_ex, list(strat_ex))
		
		av_inno_sw= c(av_inno_sw, list(inno_sw))
		av_imi_sw= c(av_imi_sw, list(imi_sw))
		av_strat_sw= c(av_strat_sw, list(strat_sw))
		
		
		
	  }
	  
	  
	  
	  

	  if(!is.null(av_inno_en[[1]])){
	  pdf(paste(path_figs_base,"/rates_inno.pdf",sep=""))
	 
	  plot(par,lapply(av_inno_en, FUN=mean),type="l", ylim = range(lapply(av_inno_en, FUN=mean), lapply(av_inno_ex, FUN=mean),lapply(av_inno_sw, FUN=mean) ))
	  lines(par,lapply(av_inno_ex, FUN=mean),col=2)
	  lines(par,lapply(av_inno_sw, FUN=mean),col=3)
	
	  legend("topleft",legend=c("entry", "exit","switch"), lty=c(1,1,1),col=c(1:3))
	  dev.off()
	  
	  
	  }
	   
	  pdf(paste(path_figs_base,"/rates_imi.pdf",sep=""))
	 
	 
	   plot(par,lapply(av_imi_en, FUN=mean),type="l", ylim = range(lapply(av_imi_en, FUN=mean),lapply( av_imi_ex, FUN=mean),lapply(av_imi_sw, FUN=mean) ))
	   lines(par,lapply(av_imi_ex, FUN=mean),col=2)
	  lines(par,lapply(av_imi_sw, FUN=mean),col=3)
	 
	   legend("topleft",legend=c("entry", "exit","switch"), lty=c(1,1,1),col=c(1:3))
	  dev.off()
	   
	  pdf(paste(path_figs_base,"/rates_strat.pdf",sep=""))
	 
	   plot(par,lapply(av_strat_en, FUN=mean),type="l", ylim = range(lapply(av_strat_en, FUN=mean),lapply( av_strat_ex, FUN=mean),lapply(av_strat_sw , FUN=mean)))
	   lines(par,lapply(av_strat_ex, FUN=mean),col=2)
	  lines(par,lapply(av_strat_sw, FUN=mean),col=3)
	   legend("topleft",legend=c("entry", "exit","switch"), lty=c(1,1,1),col=c(1:3))
	  dev.off()
  
  
  
  
  }
     
     
     
     
     
     entry_exit_rates_modes = function()  {
  
  
	  av_inno_en_m1 = list()
	  av_imi_en_m1 = list()
	  av_strat_en_m1 = list()
	  
	  av_inno_ex_m1 = list()
	  av_imi_ex_m1 = list()
	  av_strat_ex_m1 = list()
	  
	  av_inno_sw_m1 = list()
	  av_imi_sw_m1 = list()
	  av_strat_sw_m1 = list()
	  
	  
	   av_inno_en_m2 = list()
	  av_imi_en_m2 = list()
	  av_strat_en_m2 = list()
	  
	  av_inno_ex_m2 = list()
	  av_imi_ex_m2 = list()
	  av_strat_ex_m2 = list()
	  
	  av_inno_sw_m2 = list()
	  av_imi_sw_m2 = list()
	  av_strat_sw_m2 = list()
  
  
  
	
	for(p in 1:length(par)){
	
	
		inno_en_m1=c()
		imi_en_m1=c()
		strat_en_m1=c()
		
		inno_ex_m1=c()
		imi_ex_m1=c()
		strat_ex_m1=c()
		
		inno_sw_m1=c()
		imi_sw_m1=c()
		strat_sw_m1=c()
		
		
		inno_en_m2=c()
		imi_en_m2=c()
		strat_en_m2=c()
		
		inno_ex_m2=c()
		imi_ex_m2=c()
		strat_ex_m2=c()
		
		inno_sw_m2=c()
		imi_sw_m2=c()
		strat_sw_m2=c()


		eval(parse(text=paste("path ='",path_part1,par[p],path_part2,"'",sep="")))
		
		print(p)
		

		for(r in modes[[p]][[1]]){
		print(paste(path,"/run_",r,sep=""))
		
		
		notDone=FALSE
		
		if(stillRunning){
		
		  if(!file.exists(paste(path,"/run_",r,"/output.txt",sep=""))){
		  
		  notDone = TRUE
		  print("Not done yet")
		  
		  }else{
		  
		  
		  
		  da=read.csv(paste(path,"/run_",r,"/output.txt",sep=""))
		  
		 
		  da=da[,1]
		  da = paste(da,collapse=" ")
		  
		  print(grepl("Batch Done",da))
		  if(!grepl("Batch Done",da)){
		  
		  notDone = TRUE
		  print("Not done yet")
		  
		  }

		}
		
		
		}
		
		if(!notDone){
		
		
		
		
		con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
		
		
		
		#con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
		
		
		firms = dbReadTable(con,"Firm")
		exits = dbReadTable(con,"Exit")
		entry = dbReadTable(con,"Entry")
		
		
		
		
		id_max = max(firms$firmID)
		
		
		
		for(i in 1:id_max){
		
		
		
		
		    firm  = firms[firms$firmID==i,]
		    
		    exitsFirm  = exits[exits$firmID==i  & exits$exit == 1,]
		    entriesFirm = entry[entry$firmID==i & entry$entry == 1 & entry$marketEntry==0,]
		    
		  
		    
		    if(length(entriesFirm$entry)== 0 && length(exitsFirm$exit) == 0){
		    
		    
		    switch = 0
		    ent = 0
		    ex = 0
		    
		    
		    
		    }else if(length(entriesFirm$entry)> 0 && length(exitsFirm$exit) == 0){
		    
		    
		    switch = 0
		    ex = 0 
		    ent = sum (entriesFirm$entry) 
		    
		    
		    }else if(length(entriesFirm$entry)== 0 && length(exitsFirm$exit) > 0){
		    
		    
		    switch = 0
		    ex = sum (exitsFirm$exit)
		    ent = 0 
		    
		    
		    }else{
		    
		    
		    switch = 0
		    
		    for(e in 1: length(entriesFirm$entry)){
		    
		 
			  for(f in 1: length(exitsFirm$exit)){
			  
			   
		    
				if(exitsFirm$X_ITERATION_NO[f]==entriesFirm$X_ITERATION_NO[e]){
				
				    switch = switch + 1
				    
				    entriesFirm$entry[e]=0
				    exitsFirm$exit[f] = 0
				
				}

		      }
		    
		      }
		      
		      
		       ent = sum(entriesFirm$entry) 
			ex = sum (exitsFirm$exit)    
			
		      
		    
		    }
		    
		    ent = ent/ length(firm$firmID)
		    ex = ex   / length(firm$firmID)
		    switch = switch / length(firm$firmID)
		    
		    
		    #eval(parse(text=paste("data = ",FUN,"(firm$",variable,")",sep="")))
		   
		
		    
		    if(firm$firmStratgy[1] == 1){
		    
		    strat_ex_m1 = c(strat_ex_m1,ex)
		    strat_en_m1 = c(strat_en_m1,ent)
		    strat_sw_m1 = c(strat_sw_m1,switch)
		    
		    }else{
		    
		    if(firm$innovator[1] == 1){
		    
		     
		    inno_ex_m1 = c(inno_ex_m1,ex)
		    inno_en_m1 = c(inno_en_m1,ent)
		    inno_sw_m1 = c(inno_sw_m1,switch)
		    
		    }else{
		    
		     imi_ex_m1 = c(imi_ex_m1,ex)
		    imi_en_m1 = c(imi_en_m1,ent)
		    imi_sw_m1 = c(imi_sw_m1,switch)
		    
		    }
		 
		
		
		
		
		}
		
		
		}
		
		}
		
		}
		
	
	
	
		for(r in modes[[p]][[2]]){
		print(paste(path,"/run_",r,sep=""))
		
		
		notDone=FALSE
		
		if(stillRunning){
		
		  if(!file.exists(paste(path,"/run_",r,"/output.txt",sep=""))){
		  
		  notDone = TRUE
		  print("Not done yet")
		  
		  }else{
		  
		  
		  
		  da=read.csv(paste(path,"/run_",r,"/output.txt",sep=""))
		  
		 
		  da=da[,1]
		  da = paste(da,collapse=" ")
		  
		  print(grepl("Batch Done",da))
		  if(!grepl("Batch Done",da)){
		  
		  notDone = TRUE
		  print("Not done yet")
		  
		  }

		}
		
		
		}
		
		if(!notDone){
		
		
		
		
		con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
		
		
		
		#con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
		
		
		firms = dbReadTable(con,"Firm")
		exits = dbReadTable(con,"Exit")
		entry = dbReadTable(con,"Entry")
		
		
		
		
		id_max = max(firms$firmID)
		
		
		
		for(i in 1:id_max){
		
		
		
		
		    firm  = firms[firms$firmID==i,]
		    
		    exitsFirm  = exits[exits$firmID==i  & exits$exit == 1,]
		    entriesFirm = entry[entry$firmID==i & entry$entry == 1 & entry$marketEntry==0,]
		    
		  
		    
		    if(length(entriesFirm$entry)== 0 && length(exitsFirm$exit) == 0){
		    
		    
		    switch = 0
		    ent = 0
		    ex = 0
		    
		    
		    
		    }else if(length(entriesFirm$entry)> 0 && length(exitsFirm$exit) == 0){
		    
		    
		    switch = 0
		    ex = 0 
		    ent = sum (entriesFirm$entry) 
		    
		    
		    }else if(length(entriesFirm$entry)== 0 && length(exitsFirm$exit) > 0){
		    
		    
		    switch = 0
		    ex = sum (exitsFirm$exit)
		    ent = 0 
		    
		    
		    }else{
		    
		    
		    switch = 0
		    
		    for(e in 1: length(entriesFirm$entry)){
		    
		 
			  for(f in 1: length(exitsFirm$exit)){
			  
			   
		    
				if(exitsFirm$X_ITERATION_NO[f]==entriesFirm$X_ITERATION_NO[e]){
				
				    switch = switch + 1
				    
				    entriesFirm$entry[e]=0
				    exitsFirm$exit[f] = 0
				
				}

		      }
		    
		      }
		      
		      
		       ent = sum(entriesFirm$entry) 
			ex = sum (exitsFirm$exit)    
			
		      
		    
		    }
		    
		    ent = ent/ length(firm$firmID)
		    ex = ex   / length(firm$firmID)
		    switch = switch / length(firm$firmID)
		    
		    
		    #eval(parse(text=paste("data = ",FUN,"(firm$",variable,")",sep="")))
		   
		
		    
		    if(firm$firmStratgy[1] == 1){
		    
		    strat_ex_m2 = c(strat_ex_m2,ex)
		    strat_en_m2 = c(strat_en_m2,ent)
		    strat_sw_m2 = c(strat_sw_m2,switch)
		    
		    }else{
		    
		    if(firm$innovator[1] == 1){
		    
		     
		    inno_ex_m2 = c(inno_ex_m2,ex)
		    inno_en_m2 = c(inno_en_m2,ent)
		    inno_sw_m2 = c(inno_sw_m2,switch)
		    
		    }else{
		    
		     imi_ex_m2 = c(imi_ex_m2,ex)
		    imi_en_m2 = c(imi_en_m2,ent)
		    imi_sw_m2 = c(imi_sw_m2,switch)
		    
		    }
		 
		
		
		
		
		}
		
		
		}
		
		}
		
		}
		
		
		av_inno_en_m1= c(av_inno_en_m1, list(inno_en_m1))
		av_imi_en_m1= c(av_imi_en_m1, list(imi_en_m1))
		av_strat_en_m1= c(av_strat_en_m1, list(strat_en_m1))
		
		av_inno_ex_m1= c(av_inno_ex_m1, list(inno_ex_m1))
		av_imi_ex_m1= c(av_imi_ex_m1, list(imi_ex_m1))
		av_strat_ex_m1= c(av_strat_ex_m1, list(strat_ex_m1))
		
		av_inno_sw_m1= c(av_inno_sw_m1, list(inno_sw_m1))
		av_imi_sw_m1= c(av_imi_sw_m1, list(imi_sw_m1))
		av_strat_sw_m1= c(av_strat_sw_m1, list(strat_sw_m1))
		
		
		
		av_inno_en_m2= c(av_inno_en_m2, list(inno_en_m2))
		av_imi_en_m2= c(av_imi_en_m2, list(imi_en_m2))
		av_strat_en_m2= c(av_strat_en_m2, list(strat_en_m2))
		
		av_inno_ex_m2= c(av_inno_ex_m2, list(inno_ex_m2))
		av_imi_ex_m2= c(av_imi_ex_m2, list(imi_ex_m2))
		av_strat_ex_m2= c(av_strat_ex_m2, list(strat_ex_m2))
		
		av_inno_sw_m2= c(av_inno_sw_m2, list(inno_sw_m2))
		av_imi_sw_m2= c(av_imi_sw_m2, list(imi_sw_m2))
		av_strat_sw_m2= c(av_strat_sw_m2, list(strat_sw_m2))
		
		
		
	  }
	 

	  
	  pdf(paste(path_figs_base,"/MODES_rates_inno.pdf",sep=""))
	  par(mfrow=c(2,1))
	  plot(par,lapply(av_inno_en_m1, FUN=mean),type="l", ylim = range(c(0.0,0.01)))
	  lines(par,lapply(av_inno_ex_m1, FUN=mean),col=2)
	  lines(par,lapply(av_inno_sw_m1, FUN=mean),col=3)
	  legend("topleft",legend=c("entry", "exit","switch"), lty=c(1,1,1),col=c(1:3))
	  
	  
	    plot(par,lapply(av_inno_en_m2, FUN=mean),type="l", ylim = range(c(0.0,0.01)))
	  lines(par,lapply(av_inno_ex_m2, FUN=mean),col=2)
	  lines(par,lapply(av_inno_sw_m2, FUN=mean),col=3)
	  legend("topleft",legend=c("entry", "exit","switch"), lty=c(1,1,1),col=c(1:3))
	  dev.off()
	   
	  pdf(paste(path_figs_base,"/MODES_rates_imi.pdf",sep=""))
	 
	 par(mfrow=c(2,1))
	   plot(par,lapply(av_imi_en_m1, FUN=mean),type="l", ylim = range(c(0.0,0.01)))
	   lines(par,lapply(av_imi_ex_m1, FUN=mean),col=2)
	  lines(par,lapply(av_imi_sw_m1, FUN=mean),col=3)
	 
	   legend("topleft",legend=c("entry", "exit","switch"), lty=c(1,1,1),col=c(1:3))
	   
	   
	      plot(par,lapply(av_imi_en_m2, FUN=mean),type="l", ylim = range(c(0.0,0.01)))
	   lines(par,lapply(av_imi_ex_m2, FUN=mean),col=2)
	  lines(par,lapply(av_imi_sw_m2, FUN=mean),col=3)
	 
	   legend("topleft",legend=c("entry", "exit","switch"), lty=c(1,1,1),col=c(1:3))
	  dev.off()
	   
	  pdf(paste(path_figs_base,"/MODES_rates_strat.pdf",sep=""))
	 par(mfrow=c(2,1))
	   plot(par,lapply(av_strat_en_m1, FUN=mean),type="l", ylim = range(c(0.0,0.01)))
	   lines(par,lapply(av_strat_ex_m1, FUN=mean),col=2)
	  lines(par,lapply(av_strat_sw_m1, FUN=mean),col=3)
	   legend("topleft",legend=c("entry", "exit","switch"), lty=c(1,1,1),col=c(1:3))
	   
	    plot(par,lapply(av_strat_en_m2, FUN=mean),type="l", ylim = range(c(0.0,0.01)))
	   lines(par,lapply(av_strat_ex_m2, FUN=mean),col=2)
	  lines(par,lapply(av_strat_sw_m2, FUN=mean),col=3)
	   legend("topleft",legend=c("entry", "exit","switch"), lty=c(1,1,1),col=c(1:3))
	  dev.off()
  
  
  
  
  }
  
  
  aggregate_character=function(variable="numImitators"){
  
  
	    total_data = list()
	    
	    
	     pdf(paste(path_figs_base,"/hist_",variable,".pdf",sep=""))
	     par(mfrow=c(2,2))
  
  
  for(p in 1:length(par)){
	
	
		data_full=c()


		eval(parse(text=paste("path ='",path_part1,par[p],path_part2,"'",sep="")))
		

		for(r in 1:runs){
		print(paste(path,"/run_",r,sep=""))
		
		
		notDone=FALSE
		
		if(stillRunning){
		
		  if(!file.exists(paste(path,"/run_",r,"/output.txt",sep=""))){
		  
		  notDone = TRUE
		  print("Not done yet")
		  
		  }else{
		  
		  
		  
		  da=read.csv(paste(path,"/run_",r,"/output.txt",sep=""))
		  
		 
		  da=da[,1]
		  da = paste(da,collapse=" ")
		  
		  print(grepl("Batch Done",da))
		  if(!grepl("Batch Done",da)){
		  
		  notDone = TRUE
		  print("Not done yet")
		  
		  }

		}
		
		
		}
		
		if(!notDone){
		
		
		
		
		con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
		
		
		data<-dbReadTable(con,"AggregatedData")
		
		
		   eval(parse(text=paste("da = data$",variable,sep="")))
		   
		   
		   data_full = c(data_full, as.numeric(da))
		   

		
		
		}
		
		
		}
  
  
 
	     hist(data_full, main=par[p])
	     
	     
	     total_data= c(total_data, list(data_full))
  
	}
	
	
	
	
	
	
  dev.off()
  
  pdf(paste(path_figs_base,"/agg_var_",variable,".pdf",sep=""))
  
    par(mfrow=c(2,1))
   plot(par,unlist(lapply(total_data, FUN=mean)),type="l",ylim = range(unlist(lapply(total_data, FUN=mean)) +unlist(lapply(total_data, FUN=sd)),unlist(lapply(total_data, FUN=mean)) -unlist(lapply(total_data, FUN=sd))))
	  lines( par,unlist(lapply(total_data, FUN=mean)) + unlist(lapply(total_data, FUN=sd)),lty=2)
	   lines( par,unlist(lapply(total_data, FUN=mean)) - unlist(lapply(total_data, FUN=sd)),lty=2)
  
  
  
  boxplot(total_data, names= par, main=variable)
  dev.off()

  
  
  
  
  }
  
  
  
  
  
  
  histogramm_firms=function(variable="numImitators", inno=TRUE){
  
  
	    total_data = list()
	    
	    
	     pdf(paste(path_figs_base,"/hist_firm_",variable,".pdf",sep=""))
	     par(mfrow=c(2,2))
  
  
  for(p in 1:length(par)){
	
	
		data_full=c()


		eval(parse(text=paste("path ='",path_part1,par[p],path_part2,"'",sep="")))
		

		for(r in 1:runs){
		print(paste(path,"/run_",r,sep=""))
		
		
		notDone=FALSE
		
		if(stillRunning){
		
		  if(!file.exists(paste(path,"/run_",r,"/output.txt",sep=""))){
		  
		  notDone = TRUE
		  print("Not done yet")
		  
		  }else{
		  
		  
		  
		  da=read.csv(paste(path,"/run_",r,"/output.txt",sep=""))
		  
		 
		  da=da[,1]
		  da = paste(da,collapse=" ")
		  
		  print(grepl("Batch Done",da))
		  if(!grepl("Batch Done",da)){
		  
		  notDone = TRUE
		  print("Not done yet")
		  
		  }

		}
		
		
		}
		
		if(!notDone){
		
		
		
		
		con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
		
		
		data<-dbReadTable(con,"Firm")
		
		if(inno){
		  data = data[as.numeric(data$innovator)== 1 , ]
		   
		 }else{
		    data = data[as.numeric(data$innovator)== 0, ]
		 }
		
		
		   eval(parse(text=paste("da = data$",variable,sep="")))
		   
		   
		   data_full = c(data_full, as.numeric(da))
		   

		
		
		}
		
		
		}
  
  
 
	     hist(data_full, main=par[p])
  
	}
	
	
	
	
	
	
  dev.off()
  
  }
  
  
  
  
  entry_characteristics = function(variable = "numInnovators", inno = TRUE, agg=FALLSE){
  
  
	    total_data = list()
  
  
  
  pdf(paste(path_figs_base,"/hist_entry_",variable,".pdf",sep=""))
	     par(mfrow=c(2,2))
  
	for(p in 1:length(par)){
	
	
		numfirms=c()


		eval(parse(text=paste("path ='",path_part1,par[p],path_part2,"'",sep="")))
		

		for(r in 1:runs){
		print(paste(path,"/run_",r,sep=""))
		
		
		notDone=FALSE
		
		if(stillRunning){
		
		  if(!file.exists(paste(path,"/run_",r,"/output.txt",sep=""))){
		  
		  notDone = TRUE
		  print("Not done yet")
		  
		  }else{
		  
		  
		  
		  da=read.csv(paste(path,"/run_",r,"/output.txt",sep=""))
		  
		 
		  da=da[,1]
		  da = paste(da,collapse=" ")
		  
		  print(grepl("Batch Done",da))
		  if(!grepl("Batch Done",da)){
		  
		  notDone = TRUE
		  print("Not done yet")
		  
		  }

		}
		
		
		}
		
		if(!notDone){
		
		
		
		
		con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
		
		
		data<-dbReadTable(con,"AggregatedData")
		
		entry = dbReadTable(con,"Entry")
		exit = dbReadTable(con,"Exit")
		
		firms = dbReadTable(con,"Firm")
		
		
		#eval(parse(text=paste("da = data$",variable,sep="")))
		 
		 if(inno){
		  entry = entry[as.numeric(entry$isInnovator)== 1 & entry$entry==1, ]
		   exit = exit[as.numeric(exit$isInnovator)== 1  & exit$exit==1, ]
		 }else{
		  entry = entry[as.numeric(entry$isInnovator) == 0 & entry$entry==1, ]
		   exit = exit[as.numeric(exit$isInnovator)== 0 & exit$exit==1, ]
		 }
		 
		 entry =   entry[ entry$entry == 1 & entry$marketEntry==0,]
		 
		
		 del=c()
		
		  if(length(entry[,1])>0){
		
		 for(i in 1:length(entry[,1])){
		 
		    for(j in 1:length(exit[,1])){
		    
			if(entry[i,]$X_ITERATION_NO == exit[j,]$X_ITERATION_NO && entry[i,]$firmID == exit[j,]$firmID){
			
			
			del= c(del,i)			
			break
			
			}
		    
			
		    
		    }
		 
		 
		 }
		 
		
		 
		 entry <- entry[-del,]
		   
		 }
		 
	
		
		  if(agg){  
		 for(t in as.numeric(entry$X_ITERATION_NO)){
		 
		  
		 
		   eval(parse(text=paste("numfirms = c(numfirms, as.numeric(  data[as.numeric(data$X_ITERATION_NO)==t,]$",variable,"))",sep="")))
		   
		   
		   }
		  }else{
		  
		  
		   for(t in 1:length(entry$X_ITERATION_NO)){
		   
		    eval(parse(text=paste("numfirms = c(numfirms, as.numeric(  entry[",t,",]$",variable,"))",sep="")))
		   
		  
		    
		   
		   }
		  
		  
		  
		  }  

		
		
		}
		
		
		
		
		}
  
	    
	    hist(numfirms, main=par[p])
	  
  
	}
	
	
	
	dev.off()
	
	
  
  
  }
  

  
     
average_character=function(variable="equilProfit", FUN= "present_value", filter_age = FALSE, filterTime=TRUE){


	    av_inno = list()
	    av_imi = list()
	    av_strat = list()


	    inno_gam = data.frame()
	    imi_gam = data.frame()
	    strat_gam= data.frame()
	    
	    
	    noInnovator = FALSE
	    noImitator = FALSE

	    

	for(p in 1:length(par)){
	
	
		inno=c()
		imi=c()
		strat=c()


		eval(parse(text=paste("path ='",path_part1,par[p],path_part2,"'",sep="")))
		

		for(r in 1:runs){
		print(paste(path,"/run_",r,sep=""))
		
		
		notDone=FALSE
		
		
		
			
		
		
		con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
	

		if(!dbExistsTable(con,"Firm")){
		  notDone= TRUE

                   print(paste("Not done :  ",path,sep=""))
		}else{

		 firms = dbReadTable(con,"Firm")

			
	
		if(max(firms$X_ITERATION_NO)<1000){

			
			print(max(firms$X_ITERATION_NO))
			notDone= TRUE
			print(paste("Not done :  ",path,sep=""))
		}

		}


		if(notDone){


		print(paste("NOt done: ",path,"/run_",r,"/iters.db",sep=""))

		}else{

		
		 if(variable=="profit_per_firm"){
		 
		  agg = dbReadTable(con,"AggregatedData")
		 
		 }
		
		
		
		#con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
		
		
	
	

	
		if(filterTime){
		
		
		firms = firms[firms$X_ITERATION_NO > 200,]
		
		}
		
		id_max = max(firms$firmID)
		
		
		for(i in 1:id_max){
		
		
		
		
		    firm  = firms[firms$firmID==i,]
		    
		    
		    #check if the firm is still alive at the end of the simulatiuon. if this is the case: skip
		    
		    
		   
		    
		    #if(max(firm$X_ITERATION_NO)!= itno){
		    {
		    
		  
		    if(variable=="locationCosts"){
		    
		       da1 = firm$locationCosts
		       da2 =  firm$numLocationsActive
		       
		       da = da1*da2
		    
		    
		    }else{
		    
		     eval(parse(text=paste("da = firm$",variable,sep="")))#
		     
		    }
		
		
		     
		     if(variable=="profit_per_firm"){
		     
		     
		     
			    it = firm$X_ITERATION_NO
			    
			    numfirms = agg[agg$X_ITERATION_NO==it,]$numActiveFirms
			    da = firm$equilProfit * numfirms / 5.0
		     
		   
			
		     
		     }
		     
		      
		      
		     if(!filter_age || (filter_age && length(da)>100)){
		      
		    
		    
		     eval(parse(text=paste("data = ",FUN,"(da)",sep="")))#
		     
		   
		   
		   
		    
		    if(firm$firmStratgy[1] == 1){
		    
		    strat = c(strat,data)
		    
		    strat_gam = rbind(strat_gam, data.frame(par=as.numeric(par[p]), y= data))
		    
		  
		    
		    }else{
		    
		    if(firm$innovator[1] == 1){
		    
		     
		    inno = c(inno,data)
		 
		  
		    
		    inno_gam = rbind(inno_gam, data.frame(par=as.numeric(par[p]), y= data))
		    
		    
		    
		    
		    
		    
		    }else{
		    
		     imi = c(imi,data)
		   
		     imi_gam = rbind(imi_gam, data.frame(par=as.numeric(par[p]), y= data))
		    
		    }
		 
		
		
		
		
		}
		
		
		
		}
		
		}
		}
		
		}
		}
		
		#print(paste("inno:",inno,sep=""))
		#readkey()
		
		
		  if(is.null(inno)){
		    
		    noInnovator = TRUE
		    
		    }
		    
		    
		    
		  if(is.null(imi)){
		    
		    noImitator = TRUE
		    
		    }
		    
		    
		    print(av_inno)
		
		
		av_inno= c(av_inno, list(inno))
		av_imi= c(av_imi, list(imi))
		av_strat= c(av_strat, list(strat))
		
		
		
		
		
	  }
	 
	 # print(av_inno)
	 
	 
	 if(!noInnovator){
	
	  
	  pdf(paste(path_figs_base,"/",variable,"_",FUN,"_inno.pdf",sep=""))
	  par(mfrow=c(2,1))
	 
	  
	  boxplot(av_inno, main = "Innovators", names =par)
	   print("BP12")
	  plot(par,unlist(lapply(av_inno, FUN=mean)),type="l",ylim = range(unlist(lapply(av_inno, FUN=mean)) +unlist(lapply(av_inno, FUN=sd)),unlist(lapply(av_inno, FUN=mean)) -unlist(lapply(av_inno, FUN=sd))))
	  lines( par,unlist(lapply(av_inno, FUN=mean)) + unlist(lapply(av_inno, FUN=sd)),lty=2)
	   lines( par,unlist(lapply(av_inno, FUN=mean)) - unlist(lapply(av_inno, FUN=sd)),lty=2)
	  dev.off()
	  
	  
	  }
	  
	  print("BP1")
	  
	  if(length(par)>8){
	  
	   pdf(paste(path_figs_base,"/GAM_",variable,"_",FUN,".pdf",sep=""))
	  
	   
	
	
	
	if(noInnovator){
	
	
	
	  mod1<-gam(y ~  s(par,k=6, bs="cr"), data = strat_gam)
	
	  mod3<-gam(y ~  s(par, k=6,bs="cr"), data =imi_gam)
	  print("BP13")
	
	
	m1 <- gamplot(mod1) # gam-Objekt Funktion gamplot übergeben
	
	m3 <- gamplot(mod3) 

	
	 par(mfrow=c(2,1))

	
	
	if(is.null(m3[[1]]$se)){
	
	m3[[1]]$se= 0
	
	
	}
	
	
	
	if(is.null(m1[[1]]$se)){
	
	m1[[1]]$se= 0
	
	
	}
	
	
	
	
	ra = range( m1[[1]]$fit+m1[[1]]$se+ mod1$coefficients[1], m1[[1]]$fit-m1[[1]]$se+ mod1$coefficients[1], m3[[1]]$fit+m3[[1]]$se+ mod3$coefficients[1], m3[[1]]$fit-m3[[1]]$se+ mod3$coefficients[1])
	
	
	print("BP14")

plot(m1[[1]]$x, m1[[1]]$fit+ mod1$coefficients[1],type="l",col=1, xlab=expression(kappa),ylab=variable,ylim=ra, main = "Strategic firm")# zeichnet den 1. Plot der gam-Modellierung

lines(m1[[1]]$x, m1[[1]]$fit-m1[[1]]$se+ mod1$coefficients[1],col=1,lty=2) # zeichnet das obere


lines(m1[[1]]$x, m1[[1]]$fit+m1[[1]]$se+ mod1$coefficients[1],col=1,lty=2) # zeichnet das obere




plot(m3[[1]]$x, m3[[1]]$fit+ mod3$coefficients[1],type="l",col=1, xlab=expression(kappa),ylab=variable,ylim=ra, main = "Imitators")# zeichnet den 1. Plot der gam-Modellierung

lines(m3[[1]]$x, m3[[1]]$fit-m3[[1]]$se+ mod3$coefficients[1],col=1,lty=2) # zeichnet das obere


lines(m3[[1]]$x, m3[[1]]$fit+m3[[1]]$se+ mod3$coefficients[1],col=1,lty=2) # zeichnet das obere

	
	}else if(noImitator){
	
	
	
	
	  mod1<-gam(y ~  s(par,k=6, bs="cr"), data = strat_gam)
	  mod2<-gam(y ~  s(par, k=6,bs="cr"), data = inno_gam)
	 
	  print("BP13")
	m1 <- gamplot(mod1) # gam-Objekt Funktion gamplot übergeben
	m2 <- gamplot(mod2) 
	
	
	if(is.null(m2[[1]]$se)){
	
	m2[[1]]$se= 0
	
	
	}
	
	
	
	if(is.null(m1[[1]]$se)){
	
	m1[[1]]$se= 0
	
	
	}
	
	
	
	 par(mfrow=c(2,1))
	
	
	ra = range( m1[[1]]$fit+m1[[1]]$se+ mod1$coefficients[1], m1[[1]]$fit-m1[[1]]$se+ mod1$coefficients[1],m2[[1]]$fit+m2[[1]]$se+ mod2$coefficients[1], m2[[1]]$fit-m1[[1]]$se+ mod2$coefficients[1])
	
	
	print("BP14")

plot(m1[[1]]$x, m1[[1]]$fit+ mod1$coefficients[1],type="l",col=1, xlab=expression(kappa),ylab=variable,ylim=ra, main = "Strategic firm")# zeichnet den 1. Plot der gam-Modellierung

lines(m1[[1]]$x, m1[[1]]$fit-m1[[1]]$se+ mod1$coefficients[1],col=1,lty=2) # zeichnet das obere


lines(m1[[1]]$x, m1[[1]]$fit+m1[[1]]$se+ mod1$coefficients[1],col=1,lty=2) # zeichnet das obere


plot(m2[[1]]$x, m2[[1]]$fit+ mod2$coefficients[1],type="l",col=1, xlab=expression(kappa),ylab=variable,ylim=ra, main = "Innovators")# zeichnet den 1. Plot der gam-Modellierung

lines(m2[[1]]$x, m2[[1]]$fit-m2[[1]]$se+ mod2$coefficients[1],col=1,lty=2) # zeichnet das obere


lines(m2[[1]]$x, m2[[1]]$fit+m2[[1]]$se+ mod2$coefficients[1],col=1,lty=2) # zeichnet das obere


	
	
	
	}else{
	
	
	
	  mod1<-gam(y ~  s(par,k=6, bs="cr"), data = strat_gam)
	  mod2<-gam(y ~  s(par, k=6,bs="cr"), data = inno_gam)
	  mod3<-gam(y ~  s(par, k=6,bs="cr"), data =imi_gam)
	  print("BP13")
	m1 <- gamplot(mod1) # gam-Objekt Funktion gamplot übergeben
	m2 <- gamplot(mod2) 
	m3 <- gamplot(mod3) 
	
	if(is.null(m3[[1]]$se)){
	
	m3[[1]]$se= 0
	
	
	}
	
	
	
	if(is.null(m1[[1]]$se)){
	
	m1[[1]]$se= 0
	
	
	}
	
	
	
	if(is.null(m2[[1]]$se)){
	
	m2[[1]]$se= 0
	
	
	}
	
	
	
	ra = range( m1[[1]]$fit+m1[[1]]$se+ mod1$coefficients[1], m1[[1]]$fit-m1[[1]]$se+ mod1$coefficients[1],m2[[1]]$fit+m2[[1]]$se+ mod2$coefficients[1], m2[[1]]$fit-m1[[1]]$se+ mod2$coefficients[1], m3[[1]]$fit+m3[[1]]$se+ mod3$coefficients[1], m3[[1]]$fit-m3[[1]]$se+ mod3$coefficients[1])
	
	
	print("BP14")
	
	
	 par(mfrow=c(3,1))

plot(m1[[1]]$x, m1[[1]]$fit+ mod1$coefficients[1],type="l",col=1, xlab=expression(kappa),ylab=variable,ylim=ra, main = "Strategic firm")# zeichnet den 1. Plot der gam-Modellierung

lines(m1[[1]]$x, m1[[1]]$fit-m1[[1]]$se+ mod1$coefficients[1],col=1,lty=2) # zeichnet das obere


lines(m1[[1]]$x, m1[[1]]$fit+m1[[1]]$se+ mod1$coefficients[1],col=1,lty=2) # zeichnet das obere


plot(m2[[1]]$x, m2[[1]]$fit+ mod2$coefficients[1],type="l",col=1, xlab=expression(kappa),ylab=variable,ylim=ra, main = "Innovators")# zeichnet den 1. Plot der gam-Modellierung

lines(m2[[1]]$x, m2[[1]]$fit-m2[[1]]$se+ mod2$coefficients[1],col=1,lty=2) # zeichnet das obere


lines(m2[[1]]$x, m2[[1]]$fit+m2[[1]]$se+ mod2$coefficients[1],col=1,lty=2) # zeichnet das obere


plot(m3[[1]]$x, m3[[1]]$fit+ mod3$coefficients[1],type="l",col=1, xlab=expression(kappa),ylab=variable,ylim=ra, main = "Imitators")# zeichnet den 1. Plot der gam-Modellierung

lines(m3[[1]]$x, m3[[1]]$fit-m3[[1]]$se+ mod3$coefficients[1],col=1,lty=2) # zeichnet das obere


lines(m3[[1]]$x, m3[[1]]$fit+m3[[1]]$se+ mod3$coefficients[1],col=1,lty=2) # zeichnet das obere
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	
	
	}


dev.off()
	 }
	 
	 
	 if(!noImitator){
	 
	    print("BP2")
	    
	    pdf(paste(path_figs_base,"/",variable,"_",FUN,"_imi.pdf",sep=""))
	    par(mfrow=c(2,1))
	    boxplot(av_imi, main = "Imitators", names =par)
	    plot(par,unlist(lapply(av_imi, FUN=mean)),type="l", ylim = range(unlist(lapply(av_imi, FUN=mean)) +unlist(lapply(av_imi, FUN=sd)),unlist(lapply(av_imi, FUN=mean)) -unlist(lapply(av_imi, FUN=sd))))
	    lines( par,unlist(lapply(av_imi, FUN=mean)) + unlist(lapply(av_imi, FUN=sd)),lty=2)
	    lines( par,unlist(lapply(av_imi, FUN=mean)) - unlist(lapply(av_imi, FUN=sd)),lty=2)
	    dev.off()
	  
	  
	  }
	
	  
	   
	  pdf(paste(path_figs_base,"/",variable,"_",FUN,"_strat.pdf",sep=""))
	  par(mfrow=c(2,1))
	  boxplot(av_strat, main = "Strategic Firm", names =par)
	   plot(par,unlist(lapply(av_strat, FUN=mean)),type="l", ylim = range(unlist(lapply(av_strat, FUN=mean)) +unlist(lapply(av_strat, FUN=sd)),unlist(lapply(av_strat, FUN=mean)) -unlist(lapply(av_strat, FUN=sd))))
	   lines( par,unlist(lapply(av_strat, FUN=mean)) + unlist(lapply(av_strat, FUN=sd)),lty=2)
	   lines( par,unlist(lapply(av_strat, FUN=mean)) - unlist(lapply(av_strat, FUN=sd)),lty=2)
	  dev.off()
	 

	

	

	for(i in 1: length(par)){




		pdf(paste(path_figs_base,"/",variable,"_",par[i],"histogramm.pdf",sep=""))

		hist(av_strat[[i]], main = par[i])
		dev.off()

	}


pdf(paste(path_figs_base,"/",variable,"_",par[i],"dist.pdf",sep=""))

plot(dist(av_strat[[1]]), col=1, type="l")


	for(i in 2: length(par)){



                lines(dist(av_strat[[i]]), col=i)


        }
 
	  dev.off()
	  
	  
	  file = paste(path_figs_base,"/wilcox_test_",variable,"_",FUN,".txt",sep="")
	 
	  cat("+++++++++++++++++++++++ \t Wilcoxon Test \t +++++++++++++++++++++++++++++\n\n", file=file)
	  
	  cat(paste("Considered variable:\t",variable,"\t applied function: \t",FUN,sep=""),file=file, append=TRUE)  
	  
	  cat("\n\n", file=file, append=TRUE)
	  cat("++++++++++++++++ \t Strategic firm:\t +++++++++++++++++++ \n\n",file=file, append=TRUE)  
	  for(i in 2:length(par)){
	    for(j in 1:(i-1)){
	     
	      cat(paste("Compare ",par[j]," with ",par[i],"\n",sep=""),file=file,append=TRUE)
	      sink(file,append=TRUE)
	      print( wilcox.test(av_strat[[j]],av_strat[[i]]))
	     sink()
	    }
	  }
	  
	  
	  if(!noInnovator){
	  
	   cat("++++++++++++++++ \t Innovator firm: \t +++++++++++++++++++\n\n",file=file,append=TRUE)  
	  for(i in 2:length(par)){
	    for(j in 1:(i-1)){
	     
	     cat(paste("Compare ",par[j]," with ",par[i],"\n",sep=""),file=file,append=TRUE)
	       sink(file,append=TRUE)
		print(wilcox.test(av_inno[[j]],av_inno[[i]]))
		sink()
	    }
	  }
	  
	  
	  }
	  
	  
	  if(noImitator){
	  
	 cat("++++++++++++++++ \t Imitator firm: \t +++++++++++++++++++\n\n",file=file,append=TRUE)
	  for(i in 2:length(par)){
	    for(j in 1:(i-1)){
	     
	      cat(paste("Compare ",par[j]," with ",par[i],"\n",sep=""),file=file,append=TRUE)
	        sink(file,append=TRUE)
	      print(wilcox.test(av_imi[[j]],av_imi[[i]]))
	      sink()
	    }
	  }
	  
	  
	  }
	  


}





check_if_done= function(){


 
  for(p2 in 1:length(par2)){
	
	
	  for(p in 1:length(par)){
	



		eval(parse(text=paste("path ='",path_part1,par[p],path_part2,par2[p2],path_part3,"'",sep="")))
		

		for(r in 1:runs){
		print(paste(path,"/run_",r,sep=""))
		
		
		notDone=FALSE
		
		
		
		
		
		
		con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
		
		agg= dbReadTable(con,"AggregatedData")
		
		
		print(max(agg$X_ITERATION_NO))
		
		
		
		}
		
		}}



}



CheckDatabase <- function(Connection) {



out <- tryCatch(
{
   dbConnect(dbDriver("SQLite"),Connection)

    out <-TRUE

},
error=function(cond) {
  out <- FALSE
}
  )    
  return(out)
}





get_data_parallel=function(parameterarray){


		print(parameterarray)

		p = as.numeric(parameterarray[[2]])
		p2 = as.numeric(parameterarray[[1]])
		
		 variable=parameterarray[[3]]
		 FUN=parameterarray[[4]] 
		 filter_age = as.logical(parameterarray[[5]]) 
		 bands=as.logical(parameterarray[[6]]) 
		
		

		DATA =c()



		eval(parse(text=paste("path ='",path_part1,par[p],path_part2,par2[p2],path_part3,"'",sep="")))
		
		if(par2[p2]=="3.0"){
		 eval(parse(text=paste("path ='",path_part1,par[p],"/locationCosts_0.01",path_part3,"'",sep="")))
		}


		for(r in 1:runs){
		print(paste(path,"/run_",r,sep=""))
		
		
		notDone=FALSE
		
			
		
		if(CheckDatabase(paste(path,"/run_",r,"/iters.db",sep=""))){	
		
		
		con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
		
		}else{

			 print(paste("Data base not found   ",path,"/run_",r,"/iters.db",sep=""))

			notDone = TRUE

		}
		
		#con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))





	
	if(!notDone){	
		if(dbExistsTable(con,"Firm") ){
		
		
		firms = dbReadTable(con,"Firm")
		
		  id_max = max(firms$firmID)
		 
		
		  if(!is.finite(id_max)|| max(firms$X_ITERATION_NO)<1000){
		 

		    print("NDD    NDD     NDD") 
		    notDone = TRUE
		  
		  
		  }
		
		}else{
		
		 notDone = TRUE
			print("Does not exists")
		}
		
	

	}
	
		if(!notDone){
		
		
		
		
		
		print(id_max)
		

		
		
		firms = firms[firms$X_ITERATION_NO > 200,]
		
		
		for(i in 1:id_max){
		
		
		
		
		    firm  = firms[firms$firmID==i,]
		    
		    
		    #check if the firm is still alive at the end of the simulatiuon. if this is the case: skip
		    
		    
		    
		    if(TRUE)
		    {
		    
		    
		     eval(parse(text=paste("da = firm$",variable,sep="")))#
		     
		     
		     
		    
		      
		      
		     if(!filter_age || (filter_age && length(da)>100)){
		      
		    
		    
		     eval(parse(text=paste("data = ",FUN,"(da)",sep="")))#
		     
		   
		    
		#	print(data)	   
		    
		    if(firm$firmStratgy[1] == 1){
		    
		    
		    
		      
		    
		      datai = data.frame(par1 = as.numeric(par[p]), par2 = as.numeric(par2[p2]), type=1 ,y = data )
		      
		      DATA = rbind(DATA, datai )
		    
		    
		    }else{
		    
		    if(firm$innovator[1] == 1){
		    
		     
		      datai = data.frame(par1 = as.numeric(par[p]), par2 = as.numeric(par2[p2]),type=2, y = data )
		      
		      DATA = rbind(DATA, datai )
		    
		    }else{
		    
		      datai = data.frame(par1 =as.numeric(par[p]), par2 = as.numeric(par2[p2]), type=3, y = data )
		      
		      DATA = rbind(DATA, datai )
		    
		    }
		 
		
		
		
		
		}
		
		
		
		}
		
		}
		}
		
		}
		
		
		
		
		
		}

    
	  return(DATA)

}

     
average_character_GAM_two_parameters =function(variable="totalProfit", FUN= "present_value", filter_age = FALSE, bands=FALSE, ylab = "", enlarge_range= FALSE, ylimspec=FALSE){


	    counter=0

	    DATA_inno = c()
	    DATA_imi = c()
	    DATA_strat = c()
	    
	    
	    
	
   pararray= list(c(1,1,variable, FUN, filter_age , bands))
   
  for(p2 in 1:length(par2)){
	
	
	  for(p in 1:length(par)){
	
	pararray = c(pararray,list(list(p2,p,variable, FUN, filter_age , bands)))
   
	}
     
     
     }

  

DATA = mclapply(pararray, get_data_parallel, mc.cores = 8)

data = DATA 

DATA = data[[1]]


#rint(DATA)

for(i in 2:length(data)){

DATA= rbind(DATA,data[[i]])


}

print(DATA)


    DATA_inno = DATA[DATA$type==2,]
    DATA_imi = DATA[DATA$type==3,]
    DATA_strat = DATA[DATA$type==1,]
	
	
	
	x =c() 
	eff = c()  
	sd_u=c()
	sd_l=c()
	tot = c()
	tot_sd_l = c()
	tot_sd_u = c()




	  
	for (i in 1:length(par2)){
	
	
#	print(DATA_strat[DATA_strat$par2==par2[i],])
	    model_1 <- gam( y ~   s(par1,bs="cr") , data = DATA_strat[DATA_strat$par2==as.numeric(par2[i]),])
	    m1 <- gamplot(model_1)
	    
	    
	    dat= seq(as.numeric(par[1]), as.numeric(par[length(par)]), 0.001)
	    predicted_data = predict(model_1,data.frame(par1=dat),se.fit=TRUE)

	    
	    
	    if(is.null(m1[[1]]$se))
	      m1[[1]]$se = 0
	
	    x= cbind(x,m1[[1]]$x )
	    eff = cbind(eff,m1[[1]]$fit )
	     sd_u = cbind(sd_u,m1[[1]]$fit + m1[[1]]$se )
	     sd_l = cbind(sd_l,m1[[1]]$fit - m1[[1]]$se )
	     
	     tot_sd_l = cbind(tot_sd_l,predicted_data$fit - predicted_data$se.fit )
	     tot_sd_u = cbind(tot_sd_u,predicted_data$fit + predicted_data$se.fit)
	     tot = cbind(tot,predicted_data$fit )
	     
	     
	     
	
	}
	
	
	pdf(paste(path_figs_base,"/",variable,"_",FUN,"_spline_eff_strat.pdf",sep=""))
	
	for (i in 1:length(par2)){
	
	  if(i==1){
	  
	     plot(x[,i], eff[,i],type="l",ylab="Effect",xlab=expression(kappa), ylim=      range( as.vector(sd_u),as.vector(sd_l)), main="Strategic firm")
	     lines(x[,i], sd_u[,i],lty=2, col=i)
	      lines(x[,i], sd_l[,i],lty=2, col=i)
	  
	  }else{
	  
	  
	  
	   lines(x[,i], eff[,i], col=i)
	    lines(x[,i], sd_u[,i],lty=2, col=i)
	      lines(x[,i], sd_l[,i],lty=2, col=i)
	  
	  
	  }
	
	
	
	
	}
	
  dev.off()
  
  
  
  
	pdf(paste(path_figs_base,"/",variable,"_",FUN,"_tot_eff_strat.pdf",sep=""))
	
	for (i in 1:length(par2)){
	
	  if(i==1){
	  
	  if(enlarge_range){
	  
	    mu =   0.2*(max( as.vector(tot)) - min( as.vector(tot)))
	  

	    if(ylimspec){

		ylim= c(0.755,0.8)

	    }else{
		
	  ylim = c(min( as.vector(tot)) - mu , mu + max( as.vector(tot)))
	  
	}
	  }else{
	  
	  ylim = range( as.vector(tot))
	  
	  }
	  
	     plot(dat, tot[,i],type="l",xlab=expression(kappa), ylim= ylim, ylab= ylab)
	      lines(dat, tot_sd_l[,i],lty=2, col=i)
	      lines(dat, tot_sd_u[,i],lty=2, col=i)
	    
	  
	  }else{
	  
	  
	  
	   lines(dat, tot[,i], col=i)
	lines(dat, tot_sd_l[,i],lty=2, col=i)
	      lines(dat, tot_sd_u[,i],lty=2, col=i)
	  
	  }
	
	
	
	
	}
	
  dev.off()

	
	
	
	
	
	pdf(paste(path_figs_base,"/",variable,"_",FUN,"_spline_eff_inno.pdf",sep=""))
	
	x =c() 
	eff = c()  
	sd_u=c()
	sd_l=c()
	tot = c()
	 tot_sd_l = c()
	tot_sd_u = c() 
	for (i in 1:length(par2)){
	
	    model_1 <- gam( y ~   s(par1,bs="cr") , data = DATA_inno[DATA_inno$par2==as.numeric(par2[i]),])
	    m1 <- gamplot(model_1)
	    
	     
	    dat= seq(as.numeric(par[1]), as.numeric(par[length(par)]), 0.001)
	    predicted_data = predict(model_1,data.frame(par1=dat),se.fit=TRUE)

	
	if(is.null(m1[[1]]$se))
	  m1[[1]]$se = 0
	
	    x= cbind(x,m1[[1]]$x )
	    eff = cbind(eff,m1[[1]]$fit )
	     sd_u = cbind(sd_u,m1[[1]]$fit + m1[[1]]$se )
	     sd_l = cbind(sd_l,m1[[1]]$fit - m1[[1]]$se )
	        tot_sd_l = cbind(tot_sd_l,predicted_data$fit - predicted_data$se.fit )
	     tot_sd_u = cbind(tot_sd_u,predicted_data$fit + predicted_data$se.fit)
	     tot = cbind(tot,predicted_data$fit )
	
	}
	
	
	for (i in 1:length(par2)){
	
	  if(i==1){
	  
	     plot(x[,i], eff[,i],type="l",ylab="Effect",xlab=expression(kappa), ylim=      range( as.vector(sd_u),as.vector(sd_l)), main="Innovators")
	    lines(x[,i], sd_u[,i],lty=2, col=i)
	      lines(x[,i], sd_l[,i],lty=2, col=i)
	  
	  }else{
	  
	  
	  
	   lines(x[,i], eff[,i], col=i)
	      lines(x[,i], sd_u[,i],lty=2, col=i)
	      lines(x[,i], sd_l[,i],lty=2, col=i)
	  
	  
	  }
	
	
	
	
	}
	
  dev.off()
  
  
  
  
	pdf(paste(path_figs_base,"/",variable,"_",FUN,"_tot_eff_inno.pdf",sep=""))
	
	for (i in 1:length(par2)){
	
	
	if(enlarge_range){
	  
	  mu =   0.2*(max( as.vector(tot)) - min( as.vector(tot)))
	  
	  ylim = c(min( as.vector(tot)) - mu , mu + max( as.vector(tot)))
	  
	  }else{
	  
	  ylim = range( as.vector(tot))
	  
	  }
	
	  if(i==1){
	  
	     plot(dat, tot[,i],type="l",xlab=expression(kappa), ylim=   ylim,ylab = ylab)
	       lines(dat, tot_sd_l[,i],lty=2, col=i)
	      lines(dat, tot_sd_u[,i],lty=2, col=i)
	  
	  }else{
	  
	  
	  
	   lines(dat, tot[,i], col=i)
	   lines(dat, tot_sd_l[,i],lty=2, col=i)
	      lines(dat, tot_sd_u[,i],lty=2, col=i)
	  
	  }
	
	
	
	
	}
	
  dev.off()
  
  
  
  
	
	
	pdf(paste(path_figs_base,"/",variable,"_",FUN,"_spline_eff_imi.pdf",sep=""))
	
	x =c() 
	eff = c()  
	sd_u=c()
	sd_l=c()
	tot = c()
	tot_sd_l = c()
	tot_sd_u = c()
	  
	for (i in 1:length(par2)){
	
	    model_1 <- gam( y ~   s(par1,bs="cr") , data = DATA_imi[DATA_imi$par2==as.numeric(par2[i]),])
	    m1 <- gamplot(model_1)
	    
	     
	    dat= seq(as.numeric(par[1]), as.numeric(par[length(par)]), 0.001)
	    predicted_data = predict(model_1,data.frame(par1=dat),se.fit=TRUE)

	    
	    if(is.null(m1[[1]]$se))
	  m1[[1]]$se = 0
	
	    x= cbind(x,m1[[1]]$x )
	    eff = cbind(eff,m1[[1]]$fit )
	    sd_u = cbind(sd_u,m1[[1]]$fit + m1[[1]]$se )
	     sd_l = cbind(sd_l,m1[[1]]$fit - m1[[1]]$se )
	       tot_sd_l = cbind(tot_sd_l,predicted_data$fit - predicted_data$se.fit )
	     tot_sd_u = cbind(tot_sd_u,predicted_data$fit + predicted_data$se.fit)
	     tot = cbind(tot,predicted_data$fit )
	
	}
	
	
	for (i in 1:length(par2)){
	
	  if(i==1){
	  
	     plot(x[,i], eff[,i],type="l",ylab="Effect",xlab=expression(kappa), ylim=      range( as.vector(sd_u),as.vector(sd_l)), main="Imitators")
	    lines(x[,i], sd_u[,i],lty=2, col=i)
	      lines(x[,i], sd_l[,i],lty=2, col=i)
	  
	  }else{
	  
	  
	  
	   lines(x[,i], eff[,i], col=i)
    lines(x[,i], sd_u[,i],lty=2, col=i)
	      lines(x[,i], sd_l[,i],lty=2, col=i)
	  
	  
	  }
	
	
	
	
	}
	
  dev.off()
  
  
  
	pdf(paste(path_figs_base,"/",variable,"_",FUN,"_tot_eff_imi.pdf",sep=""))
	
	for (i in 1:length(par2)){
	
	  if(i==1){
	  
	  if(enlarge_range){
	  
	  mu =   0.2*(max( as.vector(tot)) - min( as.vector(tot)))
	  
	  ylim = c(min( as.vector(tot)) - mu , mu + max( as.vector(tot)))
	  
	  }else{
	  
	  ylim = range( as.vector(tot))
	  
	  }
	  
	     plot(dat, tot[,i],type="l",xlab=expression(kappa), ylim=    ylim, ylab = ylab)
	       lines(dat, tot_sd_l[,i],lty=2, col=i)
	      lines(dat, tot_sd_u[,i],lty=2, col=i)
	  
	  }else{
	  
	  
	  
	   lines(dat, tot[,i], col=i)
	   lines(dat, tot_sd_l[,i],lty=2, col=i)
	      lines(dat, tot_sd_u[,i],lty=2, col=i)
	  
	  }
	
	
	
	
	}
	
  dev.off()
  
  
  print(paste("Number of missing db: ",counter, sep=""))

	      
	      
}





     
average_character_two_parameters =function(variable="equilProfit", FUN= "present_value", filter_age = FALSE, bands=FALSE){


	    DATA_inno = c()
	    DATA_imi = c()
	    DATA_strat = c()


	for(p in 1:length(par)){
	
	
	  for(p2 in 1:length(par2)){
	


		eval(parse(text=paste("path ='",path_part1,par[p],path_part2,par2[p2],path_part3,"'",sep="")))
		

		for(r in 1:runs){
		print(paste(path,"/run_",r,sep=""))
		
		
		notDone=FALSE
		
		
		if(!notDone){
		
		
		
		con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
		
		
		
		#con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
		
		
		firms = dbReadTable(con,"Firm")
		
		id_max = max(firms$firmID)
		
		
		firms = firms[firms$X_ITERATION_NO > 200,]
		
		
		for(i in 1:id_max){
		
		
		
		
		    firm  = firms[firms$firmID==i,]
		    
		    
		    #check if the firm is still alive at the end of the simulatiuon. if this is the case: skip
		    
		    
		    
		    if(TRUE)
		    {
		    
		    
		     eval(parse(text=paste("da = firm$",variable,sep="")))#
		     
		     
		     
		    
		      
		      
		     if(!filter_age || (filter_age && length(da)>100)){
		      
		    
		    
		     eval(parse(text=paste("data = ",FUN,"(da)",sep="")))#
		     
		   
		    
		   
		    
		    if(firm$firmStratgy[1] == 1){
		    
		    
		    
		      
		    
		      datai = data.frame(par1 = par[p], par2 = par2[p2], y = data )
		      
		      DATA_strat = rbind(DATA_strat, datai )
		    
		    
		    }else{
		    
		    if(firm$innovator[1] == 1){
		    
		     
		      datai = data.frame(par1 = par[p], par2 = par2[p2], y = data )
		      
		      DATA_inno = rbind(DATA_inno, datai )
		    
		    }else{
		    
		      datai = data.frame(par1 = par[p], par2 = par2[p2], y = data )
		      
		      DATA_imi = rbind(DATA_imi, datai )
		    
		    }
		 
		
		
		
		
		}
		
		
		
		}
		
		}
		}
		
		}
		}
		
		
		
		
		
		}
		
	  }
	 
	 
	
	
	 bp_data =c() 
	
	
	  
	 pdf(paste(path_figs_base,"/",variable,"_",FUN,"_TWO_strat.pdf",sep=""))
	
	 da1 = c()
	 da2= c()
	 da3 = c()
	 
	 for(p2 in 1:length(par2)){
	 
	   mean = c()
	   lo_ci= c()
	   hi_ci= c()
	 
	   for(p in 1:length(par)){
	   
	   
	      temp = data.frame( type = 1 ,p1=par[p], p2=par2[p2],  DATA_strat[DATA_strat$par1==par[p] & DATA_strat$par2 == par2[p2],]$y)
	      
	      bp_data =  rbind(bp_data, temp)
	   
	  
	      mean = c(mean, mean(DATA_strat[DATA_strat$par1==par[p] & DATA_strat$par2 == par2[p2],]$y ))
	      
	      
	      
	   
	       lo_ci = c(lo_ci, mean(DATA_strat[DATA_strat$par1==par[p] & DATA_strat$par2 == par2[p2],]$y )- sd(DATA_strat[DATA_strat$par1==par[p] & DATA_strat$par2 == par2[p2],]$y ))
	       hi_ci = c(hi_ci, mean(DATA_strat[DATA_strat$par1==par[p] & DATA_strat$par2 == par2[p2],]$y ) + sd(DATA_strat[DATA_strat$par1==par[p] & DATA_strat$par2 == par2[p2],]$y ))
	   
	   }
	   
	  
	    da1 = rbind(da1, mean)
	    da2 = rbind(da2, lo_ci)
	    da3 = rbind(da3, hi_ci)
	    
	    
	    }
	    
	     for(p2 in 1:length(par2)){
	  if(p2==1){
	    plot(par,da1[p2,],type="l",ylim = range(da1, da2, da3),lwd=2 , col= p2+1,xlab="Parameter", ylab= variable)
	    }else{
	     lines( par,da1[p2,],lty=1, col=p2+1,lwd=2)
	    
	    }
	     if(bands){
	      lines( par,da2[p2,],lty=2, col=p2+1)
	      lines( par,da3[p2,],lty=2, col=p2+1)
	    }
	  
	  
	  
	      }
	 
	      dev.off()
	      
	      
	      ############       Boxplot    ##############
	    
	
	dat= c()
	col= c()
	
	   for(p in 1:length(par)){    
	       
	 for(p2 in 1:length(par2)){
	 
		dat = rbind(dat,DATA_strat[DATA_strat$par1==par[p] & DATA_strat$par2 == par2[p2],]$y)
		col = c(col, p2+1)
		
		print(col)
	
	}
	 
	 }
	      
	   pdf(paste(path_figs_base,"/Boxplot_",variable,"_",FUN,"_TWO_strat.pdf",sep=""))
	boxplot(t(dat), col=col , main="Strategic firms", ylab= variable)
	dev.off()
	 
	 
	 if(!is.null(DATA_inno)){
	 
	 
	 pdf(paste(path_figs_base,"/",variable,"_",FUN,"_TWO_inno.pdf",sep=""))
	
	 da1 = c()
	 da2= c()
	 da3 = c()
	 
	 for(p2 in 1:length(par2)){
	 
	   mean = c()
	   lo_ci= c()
	   hi_ci= c()
	 
	   for(p in 1:length(par)){
	   
	   
	   
	  
	      mean = c(mean, mean(DATA_inno[DATA_inno$par1==par[p] & DATA_inno$par2 == par2[p2],]$y ))
	   
	       lo_ci = c(lo_ci, mean(DATA_inno[DATA_inno$par1==par[p] & DATA_inno$par2 == par2[p2],]$y )- sd(DATA_inno[DATA_inno$par1==par[p] & DATA_inno$par2 == par2[p2],]$y ))
	       hi_ci = c(hi_ci, mean(DATA_inno[DATA_inno$par1==par[p] & DATA_inno$par2 == par2[p2],]$y ) + sd(DATA_inno[DATA_inno$par1==par[p] & DATA_inno$par2 == par2[p2],]$y ))
	   
	   }
	   
	  
	    da1 = rbind(da1, mean)
	    da2 = rbind(da2, lo_ci)
	    da3 = rbind(da3, hi_ci)
	    
	    
	    }
	    
	     for(p2 in 1:length(par2)){
	  if(p2==1){
	    plot(par,da1[p2,],type="l",ylim = range(da1, da2, da3),lwd=2 ,col=p2+1, xlab="Parameter", ylab= variable)
	    }else{
	     lines( par,da1[p2,],lty=1, col=p2+1,lwd=2)
	    
	    }
	    
	    if(bands){
	      lines( par,da2[p2,],lty=2, col=p2+1)
	      lines( par,da3[p2,],lty=2, col=p2+1)
	    }
	  
	  
	  
	      }
	 
	      dev.off()
	      
	      
	      dat= c()
	    col= c()
	    
	      for(p in 1:length(par)){    
		  
	    for(p2 in 1:length(par2)){
	    
		    dat = rbind(dat,DATA_inno[DATA_inno$par1==par[p] & DATA_inno$par2 == par2[p2],]$y)
		    col = c(col, p2+1)
		    
		    print(col)
	
	}
	 
	 }
	      
	   pdf(paste(path_figs_base,"/Boxplot_",variable,"_",FUN,"_TWO_inno.pdf",sep="")) 
	boxplot(t(dat), col=col, main="Innovators", ylab= variable)
	      dev.off()
	      
	      }
	      
	      
	      if(!is.null(DATA_imi)){
	      
	      
	     pdf(paste(path_figs_base,"/",variable,"_",FUN,"_TWO_imi.pdf",sep=""))
	 
	 da1 = c()
	 da2= c()
	 da3 = c()
	 
	 for(p2 in 1:length(par2)){
	 
	   mean = c()
	   lo_ci= c()
	   hi_ci= c()
	 
	   for(p in 1:length(par)){
	   
	  
	      mean = c(mean, mean(DATA_imi[DATA_imi$par1==par[p] & DATA_imi$par2 == par2[p2],]$y ))
	   
	       lo_ci = c(lo_ci, mean(DATA_imi[DATA_imi$par1==par[p] & DATA_imi$par2 == par2[p2],]$y )- sd(DATA_imi[DATA_imi$par1==par[p] & DATA_imi$par2 == par2[p2],]$y ))
	       hi_ci = c(hi_ci, mean(DATA_imi[DATA_imi$par1==par[p] & DATA_imi$par2 == par2[p2],]$y ) + sd(DATA_imi[DATA_imi$par1==par[p] & DATA_imi$par2 == par2[p2],]$y ))
	   
	   }
	   
	  
	    da1 = rbind(da1, mean)
	    da2 = rbind(da2, lo_ci)
	    da3 = rbind(da3, hi_ci)
	    
	    
	    }
	    
	     for(p2 in 1:length(par2)){
	  if(p2==1){
	    plot(par,da1[p2,],type="l",ylim = range(da1, da2, da3),lwd=2 ,col=p2+1, xlab="Parameter", ylab= variable)
	    }else{
	     lines( par,da1[p2,],lty=1, col=p2+1,lwd=2)
	    
	    }
	    
	     if(bands){
	    lines( par,da2[p2,],lty=2, col=p2+1)
	    lines( par,da3[p2,],lty=2, col=p2+1)
	 
	  }
	  
	  
	      }
	 
	      dev.off()
	      
	      
	        
	      dat= c()
	    col= c()
	    
	      for(p in 1:length(par)){    
		  
	    for(p2 in 1:length(par2)){
	    
		    dat = rbind(dat,DATA_imi[DATA_imi$par1==par[p] & DATA_imi$par2 == par2[p2],]$y)
		    col = c(col, p2+1)
		    
		    print(col)
	
	}
	 
	 }
	      
	   pdf(paste(path_figs_base,"/Boxplot_",variable,"_",FUN,"_TWO_imi.pdf",sep="")) 
	boxplot(t(dat), col=col, main="Imitators", ylab= variable)
	      dev.off()
	
	 
	 
	 }
	 
	 
	 
	 
	  
	  

}


 
 
 scatter_plot_with_parameter=function(variable="equilProfit", FUN= "present_value", parameter= "marketEntryHazardRate"){


	    av_inno = list()
	    av_imi = list()
	    av_strat = list()
	    
	    
	    av_inno_par = list()
	    av_imi_par = list()
	    av_strat_par = list()


	    

	for(p in 1:length(par)){
	
	
		inno=c()
		imi=c()
		strat=c()
		
		inno_par=c()
		imi_par=c()
		strat_par=c()


		eval(parse(text=paste("path ='",path_part1,par[p],path_part2,"'",sep="")))
		

		for(r in 1:runs){
		print(paste(path,"/run_",r,sep=""))
		
		
		notDone=FALSE
		
		if(stillRunning){
		
		  if(!file.exists(paste(path,"/run_",r,"/output.txt",sep=""))){
		  
		  notDone = TRUE
		  print("Not done yet")
		  
		  }else{
		  
		  
		  
		  da=read.csv(paste(path,"/run_",r,"/output.txt",sep=""))
		  
		 
		  da=da[,1]
		  da = paste(da,collapse=" ")
		  
		  print(grepl("Batch Done",da))
		  if(!grepl("Batch Done",da)){
		  
		  notDone = TRUE
		  print("Not done yet")
		  
		  }

		}
		
		
		}
		
		if(!notDone){
		
		
		
		
		con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
		
		
		
		#con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
		
		
		firms = dbReadTable(con,"Firm")
		
		
		aggVars = dbReadTable(con, "AggregatedData")
		
		id_max = max(firms$firmID)
		
		
		  eval(parse(text=paste("para = aggVars$",parameter,"[1]",sep="")))
		
		
		for(i in 1:id_max){
		
		
		
		
		    firm  = firms[firms$firmID==i,]
		    
		    
		    #check if the firm is still alive at the end of the simulatiuon. if this is the case: skip
		    
		    
		   
		    
		    if(max(firm$X_ITERATION_NO)!= itno){
		    
		   
		     eval(parse(text=paste("data = ",FUN,"(firm$",variable,")",sep="")))
		   
		   
		    
		    if(firm$firmStratgy[1] == 1){
		    
		    strat = c(strat,data)
		    strat_par = c(strat_par,para)
		    
		    
		    }else{
		    
		    if(firm$innovator[1] == 1){
		    
		     
		    inno = c(inno,data)
		    inno_par = c(inno_par,para)
		    }else{
		    
		     imi = c(imi,data)
		    imi_par = c(imi_par,para)
		    }
		 
		
		
		
		
		}
		
		
		
		}
		
		
		}
		
		}
		
		
		
		
		
		}
		
		
		#print(inno)
		
		x11()
		plot(inno,inno_par )
		x11()
		plot(imi,imi_par )
		x11()
		plot(strat,strat_par )
		
		#readkey()
		
	  }
	 
	 
	 # print(av_inno)
	
	 
	  
	  

}

get_data_parallel_exit_entry=function(parameterarray){

		p = parameterarray[2]
		p2 = parameterarray[1]

		DATA =c()
		
		eval(parse(text=paste("path ='",path_part1,par[p],path_part2,par2[p2],path_part3,"'",sep="")))

                if(par2[p2]=="3.0"){


		 eval(parse(text=paste("path ='",path_part1,par[p],"/locationCosts_0.01",path_part3,"'",sep="")))



		}

		

		for(r in 1:runs){
		print(paste(path,"/run_",r,sep=""))
		
		
		notDone=FALSE
		
		
		
		
		
		
		con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
		
		
		
		#con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
		
		
		if(dbExistsTable(con,"Firm")){
		
		
		firms = dbReadTable(con,"Firm")
		
		  id_max = max(firms$firmID)
		  
	
		
		  if(!is.finite(id_max)){
		  
		    notDone = TRUE
		  
		  
		  }
		
		}else{
		
		 notDone = TRUE
		
		}
		
		
		if(!notDone){
		
		
		
		    con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
		    
		    
		    
		    #con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
		    
		    
		    firms = dbReadTable(con,"Firm")
		    exits = dbReadTable(con,"Exit")
		    entry = dbReadTable(con,"Entry")
		    
		    
		    
		    
		    
		id_max = max(firms$firmID)
		
		
		
		for(i in 1:id_max){
		
		
		
		
		    firm  = firms[firms$firmID==i,]
		    
		    exitsFirm  = exits[exits$firmID==i  & exits$exit == 1,]
		    entriesFirm = entry[entry$firmID==i & entry$entry == 1 & entry$marketEntry==0,]
		    
		  
		    
		    if(length(entriesFirm$entry)== 0 && length(exitsFirm$exit) == 0){
		    
		    
		    switch = 0
		    ent = 0
		    ex = 0
		    
		    
		    
		    }else if(length(entriesFirm$entry)> 0 && length(exitsFirm$exit) == 0){
		    
		    
		    switch = 0
		    ex = 0 
		    ent = sum (entriesFirm$entry) 
		    
		    
		    }else if(length(entriesFirm$entry)== 0 && length(exitsFirm$exit) > 0){
		    
		    
		    switch = 0
		    ex = sum (exitsFirm$exit)
		    ent = 0 
		    
		    
		    }else if(length(entriesFirm$entry)> 0 && length(exitsFirm$exit) > 0){
		    
		    
		    switch = 0
		    
		    for(e in 1: length(entriesFirm$entry)){
		    
		 
			  for(f in 1: length(exitsFirm$exit)){
			  
			   
			   
			    if(!is.na(exitsFirm$X_ITERATION_NO[f])&&!is.na(entriesFirm$X_ITERATION_NO[e]) ){
		    
				if(exitsFirm$X_ITERATION_NO[f]==entriesFirm$X_ITERATION_NO[e]){
				
				    switch = switch + 1
				    
				    entriesFirm$entry[e]=0
				    exitsFirm$exit[f] = 0
				
				}

		      }}
		    
		      }
		      
		      
		       ent = sum(entriesFirm$entry) 
			ex = sum (exitsFirm$exit)    
			
		      
		    
		    }
		    
		    ent = ent/ length(firm$firmID)
		    ex = ex   / length(firm$firmID)
		    switch = switch / length(firm$firmID)
		    
		    
		    #eval(parse(text=paste("data = ",FUN,"(firm$",variable,")",sep="")))
		   
		
		    
		    if(firm$firmStratgy[1] == 1){
		    
		       datai = data.frame(par1 = as.numeric(par[p]), par2 = as.numeric(par2[p2]), type=1, sw = switch, en = ent, ex = ex )
		      
		       DATA = rbind(DATA, datai )
		    
		    
		    }else{
		    
		    if(firm$innovator[1] == 1){
		    
		     
		    datai = data.frame(par1 = as.numeric(par[p]), par2 = as.numeric(par2[p2]), type=2, sw = switch, en = ent, ex = ex )
		      
		       DATA = rbind(DATA, datai )
		    
		    }else{
		    
		     datai = data.frame(par1 = as.numeric(par[p]), par2 = as.numeric(par2[p2]),type=3, sw = switch, en = ent, ex = ex )
		      
		       DATA = rbind(DATA, datai )
		    
		    }
		 
		
		
		
		}
		
		
		
		
		
		}
		
		
		}
		
		
		
		}

  return(DATA)


}
      
  entry_exit_rates_GAMs = function(ymin = 0.0, ymax =0.015)  {
  
  
	  DATA_inno = c()
	  DATA_imi = c()
	  DATA_strat = c()
	  
	 
   pararray= list(c(1,1))
   
  for(p2 in 1:length(par2)){
	
	
	  for(p in 1:length(par)){
	
	pararray = c(pararray,list(c(p2,p)))
   
	}
     
     
     }

     
  # print(pararray)  

DATA = mclapply(pararray, get_data_parallel_exit_entry, mc.cores =1)

data = DATA 

DATA = data[[1]]

for(i in 2:length(data)){

DATA= rbind(DATA,data[[i]])


}

    DATA_inno = DATA[DATA$type==2,]
    DATA_imi = DATA[DATA$type==3,]
    DATA_strat = DATA[DATA$type==1,]
	
	
	
	pdf(paste(path_figs_base,"/entry_rate_strat.pdf",sep=""))
	
	x =c() 
	eff = c()  
	sd_u=c()
	sd_l=c()
	  
	for (i in 1:length(par2)){
	
	    model_1 <- gam( en ~   s(par1,bs="cr") , data = DATA_strat[DATA_strat$par2==as.numeric(par2[i]),])
	    m1 <- gamplot(model_1)
	
	    x= cbind(x,m1[[1]]$x )
	    eff = cbind(eff,m1[[1]]$fit + model_1$coefficients[1] )
	    
	
	}
	
	print(x[,1])
	print(eff[,1])
	
	
	for (i in 1:length(par2)){
	
	  if(i==1){
	  
	     plot(x[,i], eff[,i],type="l",ylab="Entry Rate",xlab=expression(kappa), ylim=      range( as.vector(eff)), main="Strategic firm")
	    
	  
	  }else{
	  
	  
	  
	   lines(x[,i], eff[,i], col=i)

	  }
	
	
	
	
	}
	
  dev.off()
  
  
  
  	pdf(paste(path_figs_base,"/exit_rate_strat.pdf",sep=""))
	
	x =c() 
	eff = c()  
	sd_u=c()
	sd_l=c()
	  
	for (i in 1:length(par2)){
	
	    model_1 <- gam( ex ~   s(par1,bs="cr") , data = DATA_strat[DATA_strat$par2==as.numeric(par2[i]),])
	    m1 <- gamplot(model_1)
	
	    x= cbind(x,m1[[1]]$x )
	    eff = cbind(eff,m1[[1]]$fit + model_1$coefficients[1] )
	    
	
	}
	
	print(x[,1])
	print(eff[,1])
	
	
	for (i in 1:length(par2)){
	
	  if(i==1){
	  
	     plot(x[,i], eff[,i],type="l",ylab="Exit Rate",xlab=expression(kappa), ylim=      range( as.vector(eff)), main="Strategic firm")
	    
	  
	  }else{
	  
	  
	  
	   lines(x[,i], eff[,i], col=i)

	  }
	
	
	
	
	}
	
  dev.off()
  
  
  
  pdf(paste(path_figs_base,"/switching_rate_strat.pdf",sep=""))
	
	x =c() 
	eff = c()  
	sd_u=c()
	sd_l=c()
	  
	for (i in 1:length(par2)){
	
	    model_1 <- gam( sw ~   s(par1,bs="cr") , data = DATA_strat[DATA_strat$par2==as.numeric(par2[i]),])
	    m1 <- gamplot(model_1)
	
	    x= cbind(x,m1[[1]]$x )
	    eff = cbind(eff,m1[[1]]$fit + model_1$coefficients[1] )
	    
	
	}
	
	print(x[,1])
	print(eff[,1])
	
	
	for (i in 1:length(par2)){
	
	  if(i==1){
	  
	     plot(x[,i], eff[,i],type="l",ylab="Switching Rate", ylim= c(ymin, ymax), xlab=expression(kappa))
	    
	  
	  }else{
	  
	  
	  
	   lines(x[,i], eff[,i], col=i)

	  }
	
	
	
	
	}
	
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
	
	pdf(paste(path_figs_base,"/entry_rate_inno.pdf",sep=""))
	
	x =c() 
	eff = c()  
	sd_u=c()
	sd_l=c()
	  
	for (i in 1:length(par2)){
	
	    model_1 <- gam( en ~   s(par1,bs="cr") , data = DATA_inno[DATA_inno$par2==as.numeric(par2[i]),])
	    m1 <- gamplot(model_1)
	
	    x= cbind(x,m1[[1]]$x )
	    eff = cbind(eff,m1[[1]]$fit + model_1$coefficients[1] )
	    
	
	}
	
	print(x[,1])
	print(eff[,1])
	
	
	for (i in 1:length(par2)){
	
	  if(i==1){
	  
	     plot(x[,i], eff[,i],type="l",ylab="Entry Rate",xlab=expression(kappa), ylim=      range( as.vector(eff)), main="Innovators")
	    
	  
	  }else{
	  
	  
	  
	   lines(x[,i], eff[,i], col=i)

	  }
	
	
	
	
	}
	
  dev.off()
  
  
  
  	pdf(paste(path_figs_base,"/exit_rate_inno.pdf",sep=""))
	
	x =c() 
	eff = c()  
	sd_u=c()
	sd_l=c()
	  
	for (i in 1:length(par2)){
	
	    model_1 <- gam( ex ~   s(par1,bs="cr") , data = DATA_inno[DATA_inno$par2==as.numeric(par2[i]),])
	    m1 <- gamplot(model_1)
	
	    x= cbind(x,m1[[1]]$x )
	    eff = cbind(eff,m1[[1]]$fit + model_1$coefficients[1] )
	    
	
	}
	
	print(x[,1])
	print(eff[,1])
	
	
	for (i in 1:length(par2)){
	
	  if(i==1){
	  
	     plot(x[,i], eff[,i],type="l",ylab="Exit Rate",xlab=expression(kappa), ylim=      range( as.vector(eff)), main="Innovators")
	    
	  
	  }else{
	  
	  
	  
	   lines(x[,i], eff[,i], col=i)

	  }
	
	
	
	
	}
	
  dev.off()
  
  
  
  pdf(paste(path_figs_base,"/switching_rate_inno.pdf",sep=""))
	
	x =c() 
	eff = c()  
	sd_u=c()
	sd_l=c()
	  
	for (i in 1:length(par2)){
	
	    model_1 <- gam( sw ~   s(par1,bs="cr") , data = DATA_inno[DATA_inno$par2==as.numeric(par2[i]),])
	    m1 <- gamplot(model_1)
	
	    x= cbind(x,m1[[1]]$x )
	    eff = cbind(eff,m1[[1]]$fit + model_1$coefficients[1] )
	    
	
	}
	
	print(x[,1])
	print(eff[,1])
	
	
	for (i in 1:length(par2)){
	
	  if(i==1){
	  
	     plot(x[,i], eff[,i],type="l",ylab="Switching Rate",xlab=expression(kappa), ylim=      range( as.vector(eff)), main="Innovators firm")
	    
	  
	  }else{
	  
	  
	  
	   lines(x[,i], eff[,i], col=i)

	  }
	
	
	
	
	}
	
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
	
	pdf(paste(path_figs_base,"/entry_rate_imi.pdf",sep=""))
	
	x =c() 
	eff = c()  
	sd_u=c()
	sd_l=c()
	  
	for (i in 1:length(par2)){
	
	    model_1 <- gam( en ~   s(par1,bs="cr") , data = DATA_imi[DATA_imi$par2==as.numeric(par2[i]),])
	    m1 <- gamplot(model_1)
	
	    x= cbind(x,m1[[1]]$x )
	    eff = cbind(eff,m1[[1]]$fit + model_1$coefficients[1] )
	    
	
	}
	
	print(x[,1])
	print(eff[,1])
	
	
	for (i in 1:length(par2)){
	
	  if(i==1){
	  
	     plot(x[,i], eff[,i],type="l",ylab="Entry Rate",xlab=expression(kappa), ylim=      range( as.vector(eff)), main="Imitators")
	    
	  
	  }else{
	  
	  
	  
	   lines(x[,i], eff[,i], col=i)

	  }
	
	
	
	
	}
	
  dev.off()
  
  
  
  	pdf(paste(path_figs_base,"/exit_rate_imi.pdf",sep=""))
	
	x =c() 
	eff = c()  
	sd_u=c()
	sd_l=c()
	  
	for (i in 1:length(par2)){
	
	    model_1 <- gam( ex ~   s(par1,bs="cr") , data = DATA_imi[DATA_imi$par2==as.numeric(par2[i]),])
	    m1 <- gamplot(model_1)
	
	    x= cbind(x,m1[[1]]$x )
	    eff = cbind(eff,m1[[1]]$fit + model_1$coefficients[1] )
	    
	
	}
	
	print(x[,1])
	print(eff[,1])
	
	
	for (i in 1:length(par2)){
	
	  if(i==1){
	  
	     plot(x[,i], eff[,i],type="l",ylab="Exit Rate",xlab=expression(kappa), ylim=      range( as.vector(eff)), main="Imitators firm")
	    
	  
	  }else{
	  
	  
	  
	   lines(x[,i], eff[,i], col=i)

	  }
	
	
	
	
	}
	
  dev.off()
  
  
  
  pdf(paste(path_figs_base,"/switching_rate_imi.pdf",sep=""))
	
	x =c() 
	eff = c()  
	sd_u=c()
	sd_l=c()
	  
	for (i in 1:length(par2)){
	
	    model_1 <- gam( sw ~   s(par1,bs="cr") , data = DATA_imi[DATA_imi$par2==as.numeric(par2[i]),])
	    m1 <- gamplot(model_1)
	
	    x= cbind(x,m1[[1]]$x )
	    eff = cbind(eff,m1[[1]]$fit + model_1$coefficients[1] )
	    
	
	}
	
	print(x[,1])
	print(eff[,1])
	
	
	for (i in 1:length(par2)){
	
	  if(i==1){
	  
	     plot(x[,i], eff[,i],type="l",ylab="Switching Rate",xlab=expression(kappa), ylim=      range( as.vector(eff)), main="Imitators")
	    
	  
	  }else{
	  
	  
	  
	   lines(x[,i], eff[,i], col=i)

	  }
	
	
	
	
	}
	
  dev.off()



  }
     

 
 
 average_character_modes=function(variable="equilProfit", FUN= "present_value"){


	    av_inno_m1 = list()
	    av_imi_m1 = list()
	    av_strat_m1 = list()


	     av_inno_m2 = list()
	    av_imi_m2 = list()
	    av_strat_m2 = list()
	    

	for(p in 1:length(par)){
	
	
		inno_m1=c()
		imi_m1=c()
		strat_m1=c()
		
		
		inno_m2=c()
		imi_m2=c()
		strat_m2=c()



		eval(parse(text=paste("path ='",path_part1,par[p],path_part2,"'",sep="")))
		

		for(r in modes[[p]][[1]]){
		print(paste(path,"/run_",r,sep=""))
		
		
		notDone=FALSE
		
		if(stillRunning){
		
		  if(!file.exists(paste(path,"/run_",r,"/output.txt",sep=""))){
		  
		  notDone = TRUE
		  print("Not done yet")
		  
		  }else{
		  
		  
		  
		  da=read.csv(paste(path,"/run_",r,"/output.txt",sep=""))
		  
		 
		  da=da[,1]
		  da = paste(da,collapse=" ")
		  
		  print(grepl("Batch Done",da))
		  if(!grepl("Batch Done",da)){
		  
		  notDone = TRUE
		  print("Not done yet")
		  
		  }

		}
		
		
		}
		
		if(!notDone){
		
		
		
		
		con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
		
		
		
		#con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
		
		
		firms = dbReadTable(con,"Firm")
		
		id_max = max(firms$firmID)
		
		
		for(i in 1:id_max){
		
		
		
		
		    firm  = firms[firms$firmID==i,]
		    
		    
		    #check if the firm is still alive at the end of the simulatiuon. if this is the case: skip
		    
		    
		   
		    
		    if(max(firm$X_ITERATION_NO)!= itno){
		    
		   
		     eval(parse(text=paste("data = ",FUN,"(firm$",variable,")",sep="")))
		   
		   
		    
		    if(firm$firmStratgy[1] == 1){
		    
		    strat_m1 = c(strat_m1,data)
		    
		    
		    }else{
		    
		    if(firm$innovator[1] == 1){
		    
		     
		    inno_m1 = c(inno_m1,data)
		    
		    }else{
		    
		     imi_m1 = c(imi_m1,data)
		    
		    }
		 
		
		
		
		
		}
		
		
		
		}
		
		
		}
		
		}
		}
		
		
		
		
		
		for(r in modes[[p]][[2]]){
		print(paste(path,"/run_",r,sep=""))
		
		
		notDone=FALSE
		
		if(stillRunning){
		
		  if(!file.exists(paste(path,"/run_",r,"/output.txt",sep=""))){
		  
		  notDone = TRUE
		  print("Not done yet")
		  
		  }else{
		  
		  
		  
		  da=read.csv(paste(path,"/run_",r,"/output.txt",sep=""))
		  
		 
		  da=da[,1]
		  da = paste(da,collapse=" ")
		  
		  print(grepl("Batch Done",da))
		  if(!grepl("Batch Done",da)){
		  
		  notDone = TRUE
		  print("Not done yet")
		  
		  }

		}
		
		
		}
		
		if(!notDone){
		
		
		
		
		con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
		
		
		
		#con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
		
		
		firms = dbReadTable(con,"Firm")
		
		id_max = max(firms$firmID)
		
		
		for(i in 1:id_max){
		
		
		
		
		    firm  = firms[firms$firmID==i,]
		    
		    
		    #check if the firm is still alive at the end of the simulatiuon. if this is the case: skip
		    
		    
		   
		    
		    if(max(firm$X_ITERATION_NO)!= itno){
		    
		   
		     eval(parse(text=paste("data = ",FUN,"(firm$",variable,")",sep="")))
		   
		   
		    
		    if(firm$firmStratgy[1] == 1){
		    
		    strat_m2 = c(strat_m2,data)
		    
		    
		    }else{
		    
		    if(firm$innovator[1] == 1){
		    
		     
		    inno_m2 = c(inno_m2,data)
		    
		    }else{
		    
		     imi_m2 = c(imi_m2,data)
		    
		    }
		 
		
		
		
		
		}
		
		
		
		}
		
		
		}
		
		}
		}
		
		
		#print(inno)
		
		
		av_inno_m1= c(av_inno_m1, list(inno_m1))
		av_imi_m1= c(av_imi_m1, list(imi_m1))
		av_strat_m1= c(av_strat_m1, list(strat_m1))
		
		av_inno_m2= c(av_inno_m2, list(inno_m2))
		av_imi_m2= c(av_imi_m2, list(imi_m2))
		av_strat_m2= c(av_strat_m2, list(strat_m2))
		
		
		
	  }
	 
	 
	 # print(av_inno)
	 

	  
	  pdf(paste(path_figs_base,"/MODES_",variable,"_",FUN,"_inno.pdf",sep=""))
	  par(mfrow=c(2,2))
	  boxplot(av_inno_m1, main = "Innovators", names =par)
	  plot(par,unlist(lapply(av_inno_m1, FUN=mean)),type="l", ylim = range(unlist(lapply(av_inno_m1, FUN=mean)) + unlist(lapply(av_inno_m1, FUN=sd)),unlist(lapply(av_inno_m1, FUN=mean)) - unlist(lapply(av_inno_m1, FUN=sd))))
	  lines(unlist(lapply(av_inno_m1, FUN=mean)) + unlist(lapply(av_inno_m1, FUN=sd)),lty=2)
	  lines(unlist(lapply(av_inno_m1, FUN=mean)) - unlist(lapply(av_inno_m1, FUN=sd)),lty=2)
	  boxplot(av_inno_m2, main = "Innovators", names =par)
	  plot(par,unlist(lapply(av_inno_m2, FUN=mean)),type="l", ylim= range(unlist(lapply(av_inno_m2, FUN=mean)) + unlist(lapply(av_inno_m2, FUN=sd)),unlist(lapply(av_inno_m2, FUN=mean)) - unlist(lapply(av_inno_m2, FUN=sd))))
	  lines(unlist(lapply(av_inno_m2, FUN=mean)) + unlist(lapply(av_inno_m2, FUN=sd)),lty=2)
	  lines(unlist(lapply(av_inno_m2, FUN=mean)) - unlist(lapply(av_inno_m2, FUN=sd)),lty=2)
	  dev.off()
	   
	  pdf(paste(path_figs_base,"/MODES_",variable,"_",FUN,"_imi.pdf",sep=""))
	  par(mfrow=c(2,2))
	  boxplot(av_imi_m1, main = "Imitators", names =par)
	   plot(par,unlist(lapply(av_imi_m1, FUN=mean)),type="l", ylim=range(unlist(lapply(av_imi_m1, FUN=mean)) + unlist(lapply(av_imi_m1, FUN=sd)),unlist(lapply(av_imi_m1, FUN=mean)) - unlist(lapply(av_imi_m1, FUN=sd))))
	   lines(unlist(lapply(av_imi_m1, FUN=mean)) + unlist(lapply(av_imi_m1, FUN=sd)),lty=2)
	   lines(unlist(lapply(av_imi_m1, FUN=mean)) - unlist(lapply(av_imi_m1, FUN=sd)),lty=2)
	   boxplot(av_imi_m2, main = "Imitators", names =par)
	   plot(par,unlist(lapply(av_imi_m2, FUN=mean)),type="l", ylim = range(unlist(lapply(av_imi_m2, FUN=mean)) + unlist(lapply(av_imi_m2, FUN=sd)),unlist(lapply(av_imi_m2, FUN=mean)) - unlist(lapply(av_imi_m2, FUN=sd))))
	   lines(unlist(lapply(av_imi_m2, FUN=mean)) + unlist(lapply(av_imi_m2, FUN=sd)),lty=2)
	  lines(unlist(lapply(av_imi_m2, FUN=mean)) - unlist(lapply(av_imi_m2, FUN=sd)),lty=2)
	  dev.off()
	   
	  pdf(paste(path_figs_base,"/MODES_",variable,"_",FUN,"_strat.pdf",sep=""))
	  par(mfrow=c(2,2))
	  boxplot(av_strat_m1, main = "Strategic Firm", names =par)
	   plot(par,unlist(lapply(av_strat_m1, FUN=mean)),type="l",ylim=range(unlist(lapply(av_strat_m1, FUN=mean)) + unlist(lapply(av_strat_m1, FUN=sd)),unlist(lapply(av_strat_m1, FUN=mean)) - unlist(lapply(av_strat_m1, FUN=sd))))
	   lines(unlist(lapply(av_strat_m1, FUN=mean)) + unlist(lapply(av_strat_m1, FUN=sd)),lty=2)
	  lines(unlist(lapply(av_strat_m1, FUN=mean)) - unlist(lapply(av_strat_m1, FUN=sd)),lty=2)
	    boxplot(av_strat_m2, main = "Strategic Firm", names =par)
	   plot(par,unlist(lapply(av_strat_m2, FUN=mean)),type="l", ylim = range(unlist(lapply(av_strat_m2, FUN=mean)) + unlist(lapply(av_strat_m2, FUN=sd)),unlist(lapply(av_strat_m2, FUN=mean)) - unlist(lapply(av_strat_m2, FUN=sd))))
	   lines(unlist(lapply(av_strat_m2, FUN=mean)) + unlist(lapply(av_strat_m2, FUN=sd)),lty=2)
	  lines(unlist(lapply(av_strat_m2, FUN=mean)) - unlist(lapply(av_strat_m2, FUN=sd)),lty=2)
	  dev.off()
	  
	  

}



scenarios = c("1","2","3","4")


for(sc in scenarios){




path_figs_base = "./figs_strategy"
dir.create(path_figs_base )

path_figs_base = "./figs_strategy"
dir.create(path_figs_base )



path_base= "./its_strategy/"

path_figs_base = paste(path_figs_base,"/Scenario_",sc,"",sep="")
dir.create(path_figs_base )


path_part1 = paste(path_base, "industryScenario_",sc,"/strategyParameter_",sep="")  


path_part2= paste("/enteringCosts_3.0",sep="")


path_part3 = paste("/",sep="")

#ar=c("-0.15","-0.125","-0.1","-0.075","-0.05","-0.025","0.0","0.025","0.05","0.075","0.1","0.125","0.15","0.175","0.2","0.225","0.25","0.275","0.3","0.325","0.35","0.375","0.4","0.425","0.45")

par=c("-0.15","-0.125","-0.1","-0.075","-0.05","-0.025","0.0","0.025","0.05","0.075","0.1","0.125","0.15","0.175","0.2","0.225","0.25","0.275","0.3")


par2=c("3.0")

runs=200
itno=1000

filter_age = TRUE

discount_factor = 0.9967


entry_exit_rates_GAMs()


average_character("totalProfit",FUN="equivalent_profit_per_year", filter_age=filter_age )
average_character("numLocationsActive","mean", filter_age=filter_age)
average_character_GAM_two_parameters("directLinksInnovators",FUN="mean", filter_age=filter_age, ylab ="Locational Links" )
average_character_GAM_two_parameters("directLinksImitators",FUN="mean", filter_age=filter_age , ylab = "Locational Links")
average_character_GAM_two_parameters("marketProfit",FUN="equivalent_profit_per_year", filter_age=filter_age , ylab = "Market Profit", enlarge_range=TRUE)
