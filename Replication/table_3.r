## This script generates a linear model that explains the average number of locations with the firm type and the length of the life cycle
##
##  Note: the number of batch runs is 200. If an different number has been used, the variable  no_runs has to be changed accordingly. 
##        Furthermore: We use a relative path starting from the top level folder using the parameters as given in run_exp_industry.sh
##	  If different parameters/ values have been used, the path has to adjusted. 

library(RSQLite)
  

#################

no_runs = 200
  
  DATA = c()

for(r in 1:no_runs){

  print(r)
 
  con<-dbConnect(dbDriver("SQLite"),paste("./its_industry/fractionEffectivityImitators_0.05/rateLocationDecision_0.02/locationCosts_0.01//run_",r,"/iters.db",sep=""))
  
  
  firms = dbReadTable(con,"Firm") 
  ids = unique(firms$firmID)
  
  for(i in ids ){
  
    firm = firms[firms$firmID==i,]
    
    
    imi =  1
    
    if(firm$innovator[1]==1){
      imi = 0
    }
    
    av_locations = mean(firm$numLocationsActive)
    age = max(firm$age)
  
  
    DATA = rbind(DATA , data.frame(r= r, id=i,imitator = imi, av_locs = av_locations, age = age))
  
  }
 
  dbDisconnect(con)
  
  }

  print(summary(lm(av_locs ~ imitator + age, data= DATA)))