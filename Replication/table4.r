
  library(RSQLite)
  
 

#################


no_runs = 200
  
  DATA = c()

for(r in 1:no_runs){


  print(r)
 

  con<-dbConnect(dbDriver("SQLite"),paste("./its_industry_rho/fractionEffectivityImitators_0.05/rateLocationDecision_0.02/locationCosts_0.01//run_",r,"/iters.db",sep=""))
  
  #con<-dbConnect(dbDriver("SQLite"),paste("/home/pharting/Desktop/test/locationCosts_0.015/run_",r,"/iters.db",sep=""))
  exits=dbReadTable(con,"Exit")
  
  agg=dbReadTable(con,"aggregatedData")

  
  firms = dbReadTable(con,"Firm") 
  
  loc = dbReadTable(con,"Location")
  
  its = unique(exits$X_ITERATION_NO)
  
  
  
  for(i in its){
  
  if(i>1){
   
  firmsT = firms[firms$X_ITERATION_NO == i-1,]
  
  links = sum(firmsT$directLinksInnovators) + sum(firmsT$directLinksImitators) 
  
  
  emptyLocations = sum(loc[loc$X_ITERATION_NO == i-1,]$numberFirmsInLocation==0)
  
  print(emptyLocations)
  
  data = exits[exits$X_ITERATION_NO == i,]
  
  av_quality= agg[agg$X_ITERATION_NO == i,]$qualityIndex
  numInnoTotal= agg[agg$X_ITERATION_NO == i,]$numInnovators
  numImiTotal= agg[agg$X_ITERATION_NO == i,]$numImitators
  
  
  
  
  ids = unique(data$firmID)
  
  
  for(j in ids){
  
  
      
  
	linkind = firmsT[firmsT$firmID==j,]$directLinksInnovators + firmsT[firmsT$firmID==j,]$directLinksImitators
	
	numLocations = firmsT[firmsT$firmID==j,]$numLocationsActive
	
 
	firm = data[data$firmID==j,]
	
	
	
	
	if(sum(firm$exit)==0){
	
		    index = which(firm$npv==max(firm$npv))[1]
		    
		    relquality = firm[index,]$qualityProduct / av_quality
		    
		    dummy = 0
		    
		    if(firm[index,]$numInnovators>0 &&  firm[index,]$numImitators ==0 ){
		    
		    dummy = 1
		    
		    }
		    
		    print(numLocations)
		    
		    if(firm[index,]$isInnovator){
		  
			temp = data.frame(exit = firm[index,]$exit, innovator =firm[index,]$isInnovator ,  numinno = firm[index,]$numInnovators-1, numimi = firm[index,]$numImitators, relQuality = relquality, numInnoTotal = numInnoTotal-1, numImiTotal= numImiTotal, numFirmsLocation =  firm[index,]$numInnovators-1 + firm[index,]$numImitators, totalNumFirms = numImiTotal + numInnoTotal -1, dummy = dummy, links = links, linkind = linkind - firm[index,]$numInnovators+1 - firm[index,]$numImitators, numLocations = numLocations, emptyLocations =emptyLocations)  
		   
		   }else{
		   
		   
			temp = data.frame(exit = firm[index,]$exit, innovator =firm[index,]$isInnovator ,  numinno = firm[index,]$numInnovators, numimi = firm[index,]$numImitators-1, relQuality = relquality, numInnoTotal = numInnoTotal, numImiTotal= numImiTotal-1, numFirmsLocation =  firm[index,]$numInnovators + firm[index,]$numImitators -1, totalNumFirms = numImiTotal + numInnoTotal -1, dummy = dummy, links = links, linkind = linkind - firm[index,]$numInnovators - firm[index,]$numImitators +1, numLocations =numLocations, emptyLocations =emptyLocations)  
		   
		   
		   }
		   
	
	}else{
	
	 fi = firm[firm$exit==1,]
	 
		   relquality = fi$qualityProduct / av_quality
		   
		   
		   dummy = 0
		    
		    if(fi$numInnovators>0 &&  fi$numImitators ==0 ){
		    
		    dummy = 1
		    
		    }
		    
		    
		  if(fi$isInnovator==1){
		   temp  = data.frame(exit = fi$exit, innovator= fi$isInnovator, numinno = fi$numInnovators-1  , numimi= fi$numImitators, relQuality = relquality, numInnoTotal = numInnoTotal -1, numImiTotal= numImiTotal  ,  numFirmsLocation =  fi$numInnovators-1 + fi$numImitators, totalNumFirms = numImiTotal + numInnoTotal -1,dummy=dummy, links = links , linkind = linkind - fi$numInnovators +1 - fi$numImitators,numLocations =numLocations, emptyLocations =emptyLocations)  
		   }else{
		    temp  = data.frame(exit = fi$exit, innovator= fi$isInnovator, numinno = fi$numInnovators  , numimi= fi$numImitators-1, relQuality = relquality, numInnoTotal = numInnoTotal, numImiTotal= numImiTotal-1,numFirmsLocation =  fi$numInnovators-1 + fi$numImitators, totalNumFirms = numImiTotal + numInnoTotal -1,dummy=dummy , links = links, linkind = linkind - fi$numInnovators+1 - fi$numImitators,numLocations =numLocations, emptyLocations =emptyLocations)  
		   
		   }
		   	 
	
	}
	
	
	
	
	DATA = rbind(DATA, temp)
  
  
  }
  
  
  
  
  }
  
  }
  
  
  }
  
  DATAINNO = DATA[DATA$innovator==1,]
  DATAIMI = DATA[DATA$innovator==0,]
  
  

  print("Imitators:")
  print(summary(glm(exit ~ numinno + numimi + numinno*numimi , family = binomial(link="logit"),data=DATAIMI)))
  print("")
  print(summary(glm(exit ~ numinno + numimi + numinno*numimi + relQuality , family = binomial(link="logit"),data=DATAIMI)))

  print(summary(glm(exit ~ numinno + numimi + numinno*numimi + relQuality + totalNumFirms  + numLocations + emptyLocations, family = binomial(link="logit"),data=DATAIMI)))

  print("Innovators:")
  print(summary(glm(exit ~ numinno + numimi + numinno*numimi , family = binomial(link="logit"),data=DATAINNO)))
  print("")
  print(summary(glm(exit ~ numinno + numimi + numinno*numimi + relQuality , family = binomial(link="logit"),data=DATAINNO)))
 
  print(summary(glm(exit ~ numinno + numimi + numinno*numimi + relQuality + totalNumFirms  + numLocations + emptyLocations, family = binomial(link="logit"),data=DATAINNO)))


  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  