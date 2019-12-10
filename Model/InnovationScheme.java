import java.util.ArrayList;



class InnovationScheme{
		
		int firmID;
		
		
		
		ArrayList<InnovationItem> innovationList;
		
		InnovationScheme( int id, double innovationProbability){
			
			firmID = id;

			//Set up the times of innovations;
	
			innovationList = new ArrayList<InnovationItem>();
			
			for(int t=0; t < Model.timeHorizon; t++){
		
				
					innovationList.add(new InnovationItem(Model.random(),Model.random()*Model.maxQualityConceptProgress));
	
				
			}
			
		}
		
		
		
	}