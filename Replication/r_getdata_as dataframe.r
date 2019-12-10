readkey <- function()
{
    cat("[press [enter] to continue]")
    number <- scan(n=1)
}


library(RSQLite)

path_figs_base = "./figs_industry_rho/figs_industry/"
dir.create(path_figs_base )


#path_figs_base = paste(path_figs_base,"/test2",sep="")
#dir.create(path_figs_base )


path_base= "./its_industry/"

runs=200

itno=1000

par=c("0.01")


# Retrieve aggegate data

data_frame_aggegate = function(){

DATA <- c()

  for(p in 1:length(par)){

  for(r in 1:runs){
  
  		path= paste(path_base,"fractionEffectivityImitators_0.05/rateLocationDecision_0.02/locationCosts_",par[p],"/",sep="")

		 con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
      
		data<-dbReadTable(con,"AggregatedData")
		
		if(r==1){
		
		  print("List of Variables:")
		  print(" ")
		  print(dbListFields(con,"AggregatedData"))
		 
		}


		datai <-  data.frame( r = r, t = 1:length(data), par = par[p],data)
		
		DATA = rbind(DATA, datai)
	
		dbDisconnect(con)

  }
}
return (DATA)

}



#retrieve firm data

data_frame_firm = function(){

DATA <- c()

  for(p in 1:length(par)){

  for(r in 1:runs){
  
  
	
  
  		path= paste(path_base,"fractionEffectivityImitators_0.05/rateLocationDecision_0.02/locationCosts_",par[p],"/",sep="")

		 con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
      
		data<-dbReadTable(con,"Firm")
		
		
		if(r==1){
		
		  print("List of Variables:")
		  print(" ")
		  print(dbListFields(con,"Firm"))
		 
		}

		datai <-  data.frame( r = r, t = 1:length(data), par = par[p],data)
		
		DATA = rbind(DATA, datai)
	
		dbDisconnect(con)

  }
}
return (DATA)

}



#retrieve firm data

data_frame_location = function(){

DATA <- c()

  for(p in 1:length(par)){

  for(r in 1:runs){
  
  		path= paste(path_base,"fractionEffectivityImitators_0.05/rateLocationDecision_0.02/locationCosts_",par[p],"/",sep="")

		 con<-dbConnect(dbDriver("SQLite"),paste(path,"/run_",r,"/iters.db",sep=""))
      
		data<-dbReadTable(con,"Location")
		
		
		
		if(r==1){
		
		  print("List of Variables:")
		  print(" ")
		  print(dbListFields(con,"Location"))
		 
		}

		datai <-  data.frame( r = r, t = 1:length(data), par = par[p],data)
		
		DATA = rbind(DATA, datai)
	
		dbDisconnect(con)

  }
}
return (DATA)

}


