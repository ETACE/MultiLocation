
## This script generates a logit model that explains the entry decision of a new-to-the-market firm depedning omn the total number of firms in a location, the number of innovators and imitators, and rero firms in location
##
##  Note: the number of batch runs is 200. If an different number has been used, the variable  no_runs has to be changed accordingly. 
##        Furthermore: We use a relative path starting from the top level folder using the parameters as given in run_exp_industry.sh
##	  If different parameters/ values have been used, the path has to adjusted. 

library(RSQLite)
library(mlogit)
  
  no_runs = 200
  
  
   DATA_INNO = c()
   DATA_IMI = c()

for(r in 1:no_runs){


  print(r)
  
  
 
  con<-dbConnect(dbDriver("SQLite"),paste("./its_industry_rho/fractionEffectivityImitators_0.05/rateLocationDecision_0.02/locationCosts_0.01//run_",r,"/iters.db",sep=""))
  #con<-dbConnect(dbDriver("SQLite"),paste("/home/pharting/workspace/DawidColombo/Paderborn/its_industry/locationCosts_0.015/rateLocationDecision_0.02/enteringCosts_3.0/run_",r,"/iters.db",sep=""))
  entry=dbReadTable(con,"Entry")
  
  
  entry = entry[entry$marketEntry==1,]
  
 
  entry_INNO = entry[entry$isInnovator==1,]
  entry_IMI = entry[entry$isInnovator==0,]
  
  empty_inno = entry_INNO$numCompetitors
  empty_inno[empty_inno>0]<-1
  

  
  
  DATA_INNO = rbind(DATA_INNO, data.frame( entry = entry_INNO$entry , no_firms = entry_INNO$numCompetitors, no_inno = entry_INNO$numInnovators  , no_imi= entry_INNO$numImitators  , non_empty= empty_inno))
  
  
  empty_imi = entry_IMI$numCompetitors
  empty_imi[empty_imi>0]<-1
  
  DATA_IMI = rbind(DATA_IMI, data.frame(entry = entry_IMI$entry , no_firms = entry_IMI$numCompetitors, no_inno = entry_IMI$numInnovators  , no_imi= entry_IMI$numImitators  , non_empty= empty_imi))
  
  
  dbDisconnect(con)
  
  }
  
  
  locID = c()
  
  id= c()
  for(i in 1:sum(DATA_INNO$entry)){
  
    id = c(id, rep(i,5))
    
    locID = c(locID,  sample(5))
 
  }
  
  DATA_INNO = data.frame(id=id, locID = locID, DATA_INNO)
  
   locID = c()
  id= c()
  for(i in 1:sum(DATA_IMI$entry)){
  
    id = c(id, rep(i,5))
  locID = c(locID,  sample(5))
  }
  
  DATA_IMI = data.frame(id=id, locID = locID ,DATA_IMI)
  
  
  
  datainno = mlogit.data(DATA_INNO, choice = "entry", shape="long",alt.var="locID", varying= c(4:7))
  dataimi = mlogit.data(DATA_IMI, choice = "entry", shape="long",alt.var="locID", varying= c(4:7))
  
  #datainno = mlogit.data(DATA_INNO, choice = "entry", shape="long", alt.var ="locID")
  
  
  print(summary(mlogit(entry ~ no_firms+ non_empty |-1, datainno)))
  print(summary(mlogit(entry ~ no_inno + no_imi + non_empty |-1, datainno)))
  
  print(summary(mlogit(entry ~  no_firms+ non_empty |-1, dataimi)))
  print(summary(mlogit(entry ~ no_inno + no_imi + non_empty |-1, dataimi)))
  