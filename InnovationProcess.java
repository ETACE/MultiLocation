

/*This class holds the innovation process*/
public class InnovationProcess implements Cloneable{
	
	public InnovationProcess clone()  {
        try {
			return (InnovationProcess) super.clone();
		} catch (CloneNotSupportedException e) {
			
			e.printStackTrace();
			return null;
		}
		
        
	}
	
	int timeToMarket;
	
	double qualityConcect;
	
	
	/*Constructor*/
	InnovationProcess(int locID, double prevQualityConcept, boolean innovation){
		
		
		if(innovation){
			
			//locationid = locid; // In which location
			timeToMarket = Model.timeToMarket;  // Set the time until the innovation is introduced
			qualityConcect = prevQualityConcept*(1+Model.random()*Model.maxQualityConceptProgress);  // Set the new quality
			
		}else{ // if it is the imitation of the best quality , then...
			
			if(locID==0){ // In Montecarlo analysis without random draw
				
				timeToMarket = Model.timeToMarket;  // Set the time until the innovation is introduced
				qualityConcect = prevQualityConcept ;  // Set the new quality
			}else if(locID==-1){
				
				
				timeToMarket = Model.timeToMarket+Model.timeToImitation;  // Set the time until the innovation is introduced
				qualityConcect = prevQualityConcept ;  // Set the new quality
				
			}else{
							
			
				timeToMarket = Model.timeToMarket;  // Set the time until the innovation is introduced
				qualityConcect = prevQualityConcept   + Model.randomNormal(0, Model.sdNoiseImitationConcept);  // Set the new quality
			
		}
		
	}
		
	}
	
	
	
	
	
	/*Constructor*/
	InnovationProcess(int locID, double prevQualityConcept, boolean innovation, double random){
		
		
		if(innovation){  //Called for innovations in the Montecarlo analysis
			
			//locationID = locID; // In which location
			timeToMarket = Model.timeToMarket;  // Set the time until the innovation is introduced
			qualityConcect = prevQualityConcept*(1+random);  // Set the new quality
			
		}else{ // if it is the imitation of the best quality , then... Never used -- can be deleted 
			
			if(locID==0){
				//locationID = locID; // In which location
				timeToMarket = Model.timeToMarket;  // Set the time until the innovation is introduced
				qualityConcect = prevQualityConcept ;  // Set the new quality
			}else if(locID==-1){
				
				//locationID = locID; // In which location
				timeToMarket = Model.timeToMarket+Model.timeToImitation;  // Set the time until the innovation is introduced
				qualityConcect = prevQualityConcept ;  // Set the new quality
				
			}else{
							
				//locationID = locID; // In which location
				timeToMarket = Model.timeToMarket;  // Set the time until the innovation is introduced
				qualityConcect = prevQualityConcept   + random;  // Set the new quality
			
		}
		
	}
		
	}
	

}
