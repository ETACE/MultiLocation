import java.util.ArrayList; 

/**
 * This class defines the innovation scheme for a particular location. Thus, it defines at which point in
 * time, which firm imitates which one within location locID
 * */

public class ImitationScheme {
	

	int locationID;
	
	ArrayList<Firm> firmList;

	ArrayList<Imitation> imitationsInLocationList;
	
	/*Constructor*/
	ImitationScheme( int locID, ArrayList<Firm> fiList){
		
		locationID = locID;
		firmList = fiList;
		imitationsInLocationList = new ArrayList<Imitation>();
		
		//Set up the times of imitations;
		// i -> imitator
		for(int i=0; i < firmList.size(); i++){
			// j -> imitated
			for(int j=0; j < firmList.size(); j++){

				if(i!=j){
					//Add imitation item
					imitationsInLocationList.add(new Imitation(firmList.get(i).firmID,firmList.get(j).firmID,firmList.get(i).currentImiProb));

				}
				
			}
	
			imitationsInLocationList.add(new Imitation(firmList.get(i).firmID,0,firmList.get(i).currentImiProb)); // For the public knowledge stock (not used anymore)
		}
		
	}
	
	
	/*Internal class for imitation items; defines at which time exactly imitation takes place */
	class Imitation{
		
		
		int imitatorID;
		int imitatedID;
		ArrayList<Integer> imitationList;
		
		
		/*Constructor*/
		Imitation(int imitator, int imitated, double imitationProbability){
			
			imitatorID = imitator;
			imitatedID = imitated;
			
			imitationList = new ArrayList<Integer>();
			
			for(int t=0; t< Model.timeHorizon;t++){
	
				//Random draws for imitation at time t
				if(Model.random() < imitationProbability){
					
					imitationList.add(t);
					
				}

			}
	
		}
		
		
	}

}
