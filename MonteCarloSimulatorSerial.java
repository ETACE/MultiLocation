
import java.io.File;
import java.util.ArrayList;

/**This Class is used to compute the NPV for the different location strategies. This is done under the assumption that 
 * the location profiles remain constant over the MC except for the location choice under study. Furthermore, the innovation and 
 * imitation activities of all firms follow a schedule that is the same for the r-th MC run regardless of the strategy under study. 
 * This is done to reduce the influence of random noise on the decision making of firms  */

public class MonteCarloSimulatorSerial {

	Firm aFirm;
	protected ArrayList<FirmReference> competitorList;
	protected ArrayList<Firm> firmList;
	protected ArrayList<Location> locationList;
	
	protected ArrayList<InnoContainer> innovationScheme;
	protected ArrayList<ImiContainer> imitationScheme;
	boolean pendingEntry = false;
	int entryID;  
	SQLMontecarlo dataBase;

	/* Constructor */
	MonteCarloSimulatorSerial(Firm bFirm, ArrayList<FirmReference> compList) {

		aFirm = bFirm.clone();
		aFirm.competitorList = deepCopyCompetitorList(aFirm.competitorList);
		competitorList = deepCopyCompetitorList(compList);
		innovationScheme = new ArrayList<InnoContainer>();
		imitationScheme = new ArrayList<ImiContainer>();
		firmList = new ArrayList<Firm>();
		locationList = new ArrayList<Location>();
		pendingEntry = false;

	}
	
	
	/**This method  sets up the innovation and imitation schemes for time horizon X montecarlo runs
	 * Note that we use a fixed scheme for each batch of simulations. This means that in each scenario, the r-th run has the same innovation and
	 * imitation profile. This is to reduce the influence of the stochastic noise*/
	void setup(){
		
		
		initMontecarlo();
		
		for(int r =0; r < Model.numMontecarloRuns;r++){
			
			// For innovations
			innovationScheme.add(new InnoContainer());
			
			for(int i=0; i < firmList.size();i++){
		
					innovationScheme.get(innovationScheme.size()-1).listOfEvents.add(new InnovationScheme(firmList.get(i).firmID,firmList.get(i).currentInnoProb));  

			}
			
			//For imitations
			imitationScheme.add(new ImiContainer());
			
			for(int i =0; i < Model.locationList.size();i++){
				
				imitationScheme.get(imitationScheme.size()-1).listOfEvents.add(new ImitationScheme(Model.locationList.get(i).locationID, firmList));

			}
			
			
			
		}
		


		/*Debug**************************************/
		if(Model.printDebugLocation){  
			
			/*First: check if data base for that firm already exists*/
			
			File f = new File("montecarlo_"+aFirm.firmID+".db");
			if(!f.exists()) { 
				
				
				dataBase = new SQLMontecarlo("montecarlo_"+aFirm.firmID+".db",aFirm.firmID);
				
				dataBase.createEntryExitTable();
				dataBase.createEntryTable();
				dataBase.createExitTable();
				dataBase.createBaselineTable();
			    
				
				
			}else{
	
				dataBase = new SQLMontecarlo("montecarlo_"+aFirm.firmID+".db");
				
			}
			
		}

	}
	
	/**This method  is used to initialize the monte carlo simulations. To do so we first delete the old lists and then set up the firms 
	 * in a way that changing the memory of the instances during the MC does not affect the memory of the actual agents. (deep cpying)*/
	void initMontecarlo() {

		firmList.clear();
		locationList.clear();
		
		
		/*Create deep copies of firms and list components. This is to avoid that the actual firm data is overwritten if the methods of
		 * of the agents are called.*/
		
		//Start with firm under consideration
		Firm copyFirm = (Firm) aFirm.clone();

		copyFirm.activeLocationList = deepCopyLocationList(aFirm.activeLocationList);
		copyFirm.passiveLocationList = deepCopyLocationList(aFirm.passiveLocationList);

		copyFirm.competitorList = deepCopyCompetitorList(aFirm.competitorList);

		copyFirm.pendingInnovationsList  = deepCopyPendingInnovationsList(aFirm.pendingInnovationsList);
		
		
		firmList.add(copyFirm);

		// Add competitors
		for (int i = 0; i < competitorList.size(); i++) {

			Firm bFirm = (Firm) Model.firmList.get(
					Model.returnFirmIndex(competitorList.get(i).firmID)).clone();

			bFirm.activeLocationList = deepCopyLocationList(Model.firmList
					.get(Model.returnFirmIndex(competitorList.get(i).firmID)).activeLocationList);
			bFirm.passiveLocationList = deepCopyLocationList(Model.firmList
					.get(Model.returnFirmIndex(competitorList.get(i).firmID)).passiveLocationList);

			bFirm.competitorList = deepCopyCompetitorList((Model.firmList
					.get(Model.returnFirmIndex(competitorList.get(i).firmID)).competitorList));
			
			bFirm.pendingInnovationsList  = deepCopyPendingInnovationsList(Model.firmList
					.get(Model.returnFirmIndex(competitorList.get(i).firmID)).pendingInnovationsList);

			firmList.add(bFirm);

		}

		//Edit competitor lists
		for (int i = 0; i < firmList.size(); i++) {

			Firm aFirm = firmList.get(i).clone();

			aFirm.competitorList.clear();

			for (int j = 0; j < firmList.size(); j++) {

				if (i != j) {

					Firm bFirm = firmList.get(j).clone();
					aFirm.competitorList.add(new FirmReference(bFirm.firmID,
							bFirm.quality, bFirm.qualityConcept,bFirm.innovator, bFirm.maxPendingQualityConcept));

				}

			}

		}
		
		
		// Add location lists
		for(int i=0; i <Model.locationList.size();i++){
			
			
			locationList.add(Model.locationList.get(i).clone());
			
	
		}

	}
	
	
	/**This method  is used to run the MC simulation under the current location profile of the firm*/
	double monteCarloBase() {
		

		pendingEntry = false; 

		double npv= 0;

		//To store the NPVs
		ArrayList<Double> tempNPVList = new ArrayList<Double>(); 

		/* For n Montecarlo runs */
		for (int r = 0; r < Model.numMontecarloRuns; r++) {

			//For each run the MC simulator has to be initialized
			initMontecarlo();
			
			// Call method to run the simulator for one run
			double tempNV = MontecarloRun(0,Model.timeHorizon, r+1, 0,0,0);
			
			// Add NPV to array
			tempNPVList.add(tempNV);
			
			if(Double.isNaN(tempNV)){
				System.out.println("monte Carlo base: gives NAN");
				System.exit(1);	
			}
		}
			
		/*Return either mean or average*/
		if(Model.monteCarloMedian){
			npv = Model.returnMedianDoubleList(tempNPVList);	
		}else{	
			npv = Model.returnMeanDoubleList(tempNPVList);	
		}

		return npv;
			
	
	}
	
	
	
	/**This method  is used to run the MC simulation under the assumtption that the firm additionally enters
	 * location with ID locatioID*/
	double monteCarloEntry(int locationID){
		
		/*Check entry*/
		pendingEntry = true; 
		entryID = locationID;
		
		double npv = 0;
		//To store the NPVs
		ArrayList<Double> tempNPVList = new ArrayList<Double>(); 
		
		/* For n Montecarlo runs */
		for (int r = 0; r < Model.numMontecarloRuns; r++) {
			
			//For each run the MC simulator has to be initialized
			initMontecarlo();
	
		try {
			/*Since the location profile of the firm changes, we have to adjust the active and passive location list (temporarily)*/
			for (int i = 0; i < firmList.get(0).passiveLocationList
					.size(); i++) {

				if (firmList.get(0).passiveLocationList.get(i).locationID == locationID) {

					firmList.get(0).activeLocationList.add(firmList
							.get(0).passiveLocationList.get(i));

					firmList.get(0).passiveLocationList.remove(i);
					break;
		
				}
			}
				
			/*Update the competitor arrays*/
			
			for(int i=1; i < firmList.size();i++){
				
				for(int j=0; j < firmList.get(i).activeLocationList.size();j++){
					
					if(firmList.get(i).activeLocationList.get(j).locationID==locationID){
						
						firmList.get(i).activeLocationList.get(j).competitorList.add(new FirmReference(firmList.get(0).firmID, null,firmList.get(0).qualityConcept,firmList.get(0).innovator,firmList.get(0).maxPendingQualityConcept ));
		
					}
	
				}

			}

			// Call method to run the simulator for one run
		double tempNV = MontecarloRun(0,Model.timeHorizon,r+1, 1, entryID, 0);
		
		//Add NPV to array
		tempNPVList.add(tempNV);
		
		if(Double.isNaN(tempNV)){
			System.out.println("monte Carlo entry: gives NAN");
			System.exit(1);	
			
		}

			
		} catch (IndexOutOfBoundsException e) {

			System.out.println("Index out of bound exception: ");

		}

	}
		
	/*Return either mean or median*/
	
	if(Model.monteCarloMedian){
		
		npv = Model.returnMedianDoubleList(tempNPVList);
		
	}else{
		
		npv = Model.returnMeanDoubleList(tempNPVList);
		
	}

	return npv;
}
	
	

/**This method  is used to run the MC simulation under the assumption that the firm additionally exits
 * location with ID locatioID*/
double monteCarloExit(int locationID){

	pendingEntry = false; 
	
	double npv=0;

	//Array to store the NPVS
	ArrayList<Double> tempNPVList = new ArrayList<Double>(); 
	
	/* For n Montecarlo runs */

	for (int r = 0; r < Model.numMontecarloRuns; r++) {
	
		//For each run the MC simulator has to be initialized
	initMontecarlo();

	/*Since the location profile of the firm changes, we have to adjust the active and passive location list (temporarily)*/
	try {

		/*Deactivate location*/
		for (int i = 0; i < firmList.get(0).activeLocationList.size(); i++) {

			if (firmList.get(0).activeLocationList.get(i).locationID == locationID) {

				firmList.get(0).passiveLocationList.add(firmList
						.get(0).activeLocationList.get(i));

				firmList.get(0).activeLocationList.remove(i);
				break;
				
				}

			}

		/*Adjust competitor list*/
		
		for(int i=1; i < firmList.size();i++){
			
			for(int j=0; j < firmList.get(i).activeLocationList.size();j++){
				
				if(firmList.get(i).activeLocationList.get(j).locationID==locationID){
					
					
					for(int k =0; k < firmList.get(i).activeLocationList.get(j).competitorList.size();k++ ){
						
						if(firmList.get(i).activeLocationList.get(j).competitorList.get(k).firmID==firmList.get(0).firmID){
							
							firmList.get(i).activeLocationList.get(j).competitorList.remove(k);
							break;
							
						}
							
					}
		
				}
			}
			
			
			}
			
			/*Run MC run*/
			double tempNV = MontecarloRun(0,Model.timeHorizon, r+1,2,0,locationID);
			
			//Add NPV to list
			tempNPVList.add(tempNV);
			
			if(Double.isNaN(tempNV)){
				System.out.println("monte Carlo exit: gives NAN");
				System.exit(1);	
				
			}

		
		} catch (IndexOutOfBoundsException e) {

			System.out.println("Index out of bound exception: ");

		}
	}
	
	// Return either mean or median of NPVs
	if(Model.monteCarloMedian){
		
		npv = Model.returnMedianDoubleList(tempNPVList);
		
	}else{
		
		npv = Model.returnMeanDoubleList(tempNPVList);
		
	}

	return npv;
	
}


	
/**This method  is used to run the MC simulation under the assumption that the firm switches from idLocationExit to idLocationEntry
 * location with ID locatioID*/
double monteCarloEntryExit(int idLocationExit, int idLocationEntry){
		
	
	pendingEntry = true; 
	entryID = idLocationEntry;

	double npv = 0;
	
	//Array to store the NPVS
	ArrayList<Double> tempNPVList = new ArrayList<Double>(); 
	
	/* For n Montecarlo runs */
	for (int r = 0; r < Model.numMontecarloRuns; r++) {
	
		// Init the MC simulator
		initMontecarlo();
		
		
		/*Since the location profile of the firm changes, we have to adjust the active and passive location list (temporarily)*/
		try {

			/*Firm under consideration*/
			
			/*Exit*/
			for (int i = 0; i < firmList.get(0).activeLocationList
					.size(); i++) {

				if (firmList.get(0).activeLocationList.get(i).locationID == idLocationExit) {

					firmList.get(0).passiveLocationList.add(firmList
							.get(0).activeLocationList.get(i));

					firmList.get(0).activeLocationList.remove(i);
					
					break;
				}
			}
					
			/*Entry*/		
			for (int i = 0; i < firmList.get(0).passiveLocationList
					.size(); i++) {

				if (firmList.get(0).passiveLocationList.get(i).locationID == idLocationEntry) {

					firmList.get(0).activeLocationList.add(firmList
							.get(0).passiveLocationList.get(i));

					firmList.get(0).passiveLocationList.remove(i);
					break;
			
				}
			}
				
			//Adjust the competitors
				
			for(int i=1; i < firmList.size();i++){
				
				for(int j=0; j < firmList.get(i).activeLocationList.size();j++){
					
					if(firmList.get(i).activeLocationList.get(j).locationID==idLocationEntry){
						
						firmList.get(i).activeLocationList.get(j).competitorList.add(new FirmReference(firmList.get(0).firmID, null,firmList.get(0).qualityConcept ,firmList.get(0).innovator,firmList.get(0).maxPendingQualityConcept));
							
					}
					
					if(firmList.get(i).activeLocationList.get(j).locationID==idLocationExit){
						
						
						for(int k =0; k < firmList.get(i).activeLocationList.get(j).competitorList.size();k++ ){
							
							if(firmList.get(i).activeLocationList.get(j).competitorList.get(k).firmID==firmList.get(0).firmID){
								
								firmList.get(i).activeLocationList.get(j).competitorList.remove(k);
								break;	
							}	
						}	
					}	
				}			
			}
				
		//Call MC simulator
		double tempNV =MontecarloRun(0,Model.timeHorizon,r+1,3,idLocationEntry, idLocationExit);
			
		//Add NPV to array
		tempNPVList.add(tempNV);
		
		if(Double.isNaN(tempNV)){
			System.out.println("monte Carlo entry exit: gives NAN");
			System.exit(1);	
			
		}

			
		} catch (IndexOutOfBoundsException e) {

			System.out.println("Index out of bound exception: ");

		}
	}
	
	//Return either mean or median
	if(Model.monteCarloMedian){
		
		npv = Model.returnMedianDoubleList(tempNPVList);
		
	}else{
		
		npv = Model.returnMeanDoubleList(tempNPVList);
		
	}
	
	return npv;		
}


/**This method  is executes a generic Montecarlo run*/
	double MontecarloRun(int start, int end, int run, int isIn, int entryID, int exitID) {

		double npv = 0;

		// MC-run over time horizon
		for (int t = start; t < end; t++) {

			Model.negativeOutput = false;

			
			//Go through all steps for all firms in each iteration
				for (int i = 0; i < firmList.size(); i++) {

					firmList.get(i).firmSetup();
					firmList.get(i).determineOutput();
		
				}
		
				if(Model.negativeOutput){
					
					boolean solve = solveNegativeOutput();
					
					if(!solve){
						/*If the negative output is not solved, it means the negative output is at the side of the considered firm*/
						//System.out.println("time t "+t);
						//System.out.println("Quantity "+firmList.get(0).equilQuantity+" Profit "+firmList.get(0).equilProfit+" Quality "+firmList.get(0).quality+" Concept Quality "+firmList.get(0).qualityConcept+" "+firmList.get(0).qualityContributionKnowledge);
						break;
						
					}
				}
					
					
			
			/*Adjust competitor arrays*/
			for(int i=0; i< firmList.size(); i++){
				
				Firm aFirm = (Firm) firmList.get(i);
				 for(int j=0; j < aFirm.competitorList.size();j++){
					 for(int k=0; k < firmList.size();k++){
						 if(aFirm.competitorList.get(j).firmID==firmList.get(k).firmID){
							 
							 aFirm.competitorList.get(j).output = firmList.get(k).equilQuantity;
							 break;
						 }
	 
					 }
					 
					 
				 }
				
			}	
					
			// Innovation and imitation steps of all firms according to the schedule	
			for (int i = 0; i < firmList.size(); i++) {	
				
				firmList.get(i).determinePriceProfits();

				innovationProcess(i, run, t);
				imitationProcess(i, run, t);

			}

			//Update quality
			for (int i = 0; i < firmList.size(); i++) {

				firmList.get(i).updateQuality();

			}
			

			updateFirmArrays();
		
			// Add npv
			npv += Math.pow(Model.discontFactor, t)
					* firmList.get(0).equilProfit;
			
			
			/*For debugging*/
			if(Model.printDebugLocation){  
				
				if(isIn == 0){
					
					for(int i=0; i < firmList.size();i++){
					dataBase.insertbaseline(run,Model.iteration, t, npv, firmList.get(i));
					}
					
				}else if(isIn == 1){
					
					for(int i=0; i < firmList.size();i++){
						dataBase.insertEntry(run,Model.iteration, t, entryID, npv, firmList.get(i));
						}
					
				}else if(isIn == 2){
					
					for(int i=0; i < firmList.size();i++){
						dataBase.insertEntry(run,Model.iteration, t, exitID, npv, firmList.get(i));
						}
					
				}else{
					
					for(int i=0; i < firmList.size();i++){
						dataBase.insertEntryExit(run,Model.iteration, t,entryID, exitID, npv, firmList.get(i));
						}
					
				}
			}
			
		
		}
		return npv;

	}

	/**This method  determines the innovation steps during the MC. The idea is that there is an exogeneous scheme of probabilities
	 * that is the same for a run r. Based on this, the firm checks, given the innovation probability, whether or not a innovation has been successful in t of r. In
	 * this case the quality is updated.*/	
void innovationProcess(int indexFirm, int numOfRun, int time){
		
			
		for(int i=0; i < innovationScheme.get(numOfRun-1).listOfEvents.size();i++ ){
			
			if(innovationScheme.get(numOfRun-1).listOfEvents.get(i).firmID==firmList.get(indexFirm).firmID){
				
				//Check if innovation has been realized given the innovation probability and the determined random draw. If this is the case: update quality
					if( innovationScheme.get(numOfRun-1).listOfEvents.get(i).innovationList.get(time).prob < firmList.get(indexFirm).innovationProbability())
						{
						firmList.get(indexFirm).pendingInnovationsList.add(new InnovationProcess(0,Math.max(firmList.get(indexFirm).tempQualityConcept,firmList.get(indexFirm).maxPendingQualityConcept),true, innovationScheme.get(numOfRun-1).listOfEvents.get(i).innovationList.get(time).randomdraw));
					
						break;
					}
	
				break;
			}
					
		}
	
	}
	

/**This method  determines the imitation steps during the MC. Since the location profile is fix, we define for each firm and location a imitation plan (consisting of points in time) at which 
 * time the firm imitates the quality from the other firms in the location. So we have to go through the active location list, 
 * search for each competitor and check if for the competitor at time t an imitation is scheduled. These checks require a long list of nested for-loops and if-conditions  */	
void imitationProcess(int indexFirm, int numOfRun, int time){
		

	for(int i=0; i < firmList.get(indexFirm).activeLocationList.size();i++){
		
		for(int j=0; j < imitationScheme.get(numOfRun-1).listOfEvents.size(); j++){

			//Find location
			if(firmList.get(indexFirm).activeLocationList.get(i).locationID==imitationScheme.get(numOfRun-1).listOfEvents.get(j).locationID){
				
				
				for(int k=0; k < imitationScheme.get(numOfRun-1).listOfEvents.get(j).imitationsInLocationList.size();k++){
					
					// Find firm to be imitated
					if(imitationScheme.get(numOfRun-1).listOfEvents.get(j).imitationsInLocationList.get(k).imitatorID==firmList.get(indexFirm).firmID){
						
						for(int t=0; t < imitationScheme.get(numOfRun-1).listOfEvents.get(j).imitationsInLocationList.get(k).imitationList.size();t++){
							
							//Check iof imitation is scheduled
							if(imitationScheme.get(numOfRun-1).listOfEvents.get(j).imitationsInLocationList.get(k).imitationList.get(t)==time){
								
								
									for(int l=0; l < firmList.get(indexFirm).activeLocationList.get(i).competitorList.size();l++){
										
										
										if(firmList.get(indexFirm).activeLocationList.get(i).competitorList.get(l).firmID==imitationScheme.get(numOfRun-1).listOfEvents.get(j).imitationsInLocationList.get(k).imitatedID){
					
												/*Imitation of concept quality*/
											
												if(Model.modelModeTimeToMarket)
													firmList.get(indexFirm).pendingInnovationsList.add(new InnovationProcess(0,firmList.get(indexFirm).activeLocationList.get(i).competitorList.get(l).maxPendingQualityConcept,false));
												else
													firmList.get(indexFirm).pendingInnovationsList.add(new InnovationProcess(0,firmList.get(indexFirm).activeLocationList.get(i).competitorList.get(l).qualityConcept,false));
								
											break;	
											}
										}
								
									break;
								}
							}
						}
					}
				break;	
				
				}	
			}
		}	
	}
		
		

	/**
	 * This method makes a deep copy of the FirmLocation object 
	 */
	ArrayList<FirmLocation> deepCopyLocationList(
			ArrayList<FirmLocation> firmLocationList) {

		ArrayList<FirmLocation> copyFirmLocationList = new ArrayList<FirmLocation>();

		for (int i = 0; i < firmLocationList.size(); i++) {

			FirmLocation aLocation = firmLocationList.get(i).clone();

			aLocation.competitorList = deepCopyCompetitorList(aLocation.competitorList);

			copyFirmLocationList.add(aLocation);

		}

		return copyFirmLocationList;

	}

	ArrayList<FirmReference> deepCopyCompetitorList(
			ArrayList<FirmReference> competitorList) {

		ArrayList<FirmReference> copyCompetitorList = new ArrayList<FirmReference>();

		for (int i = 0; i < competitorList.size(); i++) {

			copyCompetitorList.add(competitorList.get(i).clone());

		}

		return copyCompetitorList;

	}
	
	
	
	ArrayList<InnovationProcess> deepCopyPendingInnovationsList(
			ArrayList<InnovationProcess> pendingInnovationsList) {

		ArrayList<InnovationProcess> copyPendingInnovationsList = new ArrayList<InnovationProcess>();

		for (int i = 0; i < pendingInnovationsList.size(); i++) {

			copyPendingInnovationsList.add(pendingInnovationsList.get(i).clone());

		}

		return copyPendingInnovationsList;

	}
	
	
	

/**
 * This method is used to clean up 
 */

void destructMontecarloSimulator() {
	
	
	if(Model.printDebugLocation)
		dataBase.commit();

	aFirm = null;
	competitorList = null;
	firmList = null;	
	innovationScheme = null;
	imitationScheme = null;


}

	
/**
 * This method is used to handle those cases where the output of a firm is negative (which cannot be excluded given the 
 * market model)
 */
boolean solveNegativeOutput(){

	while(Model.negativeOutput){
		
		Model.negativeOutput = false;
		double tempqual = 99999.0;
		int tempID = 0;
			
		for(int i=0; i< firmList.size(); i++){
	
			if(firmList.get(i).active && firmList.get(i).equilQuantity<tempqual){
					
					tempqual = firmList.get(i).equilQuantity;
					tempID = firmList.get(i).firmID;

			}
		
		}
		
		for(int i=0; i< firmList.size(); i++){
			
			if(firmList.get(i).firmID==tempID){
				
				/*If index != 0 then we do not consider the firm that is doing the location decision; otherwise return false*/
				if(i!=0){
					/*One of the other firms has a negative profit*/
					firmList.remove(i);
					break;
				}else{
					
					//Firm itself has negative output
					return false;
					
					
				}
		
			}
			
		}

		for(int i=0; i< firmList.size(); i++){
			
			Firm aFirm = (Firm) firmList.get(i);
			if(aFirm.active){
	
						for(int k=0; k < aFirm.competitorList.size();k++ ){
							
							if(tempID==aFirm.competitorList.get(k).firmID){
								
								aFirm.competitorList.remove(k);	
								break;
							}
	
						}
			
				aFirm.determineOutput();
	
			}
		}
		
	
	}

return true;

}

	
/**
 * This method is used to update the firm arrays
 */
public void updateFirmArrays() {

	

	for (int i = 0; i < firmList.size(); i++) {

		Firm aFirm = (Firm) firmList.get(i);

		/*Update competitor list (by rewriting it)*/
		aFirm.competitorList.clear();

		for (int j = 0; j < firmList.size(); j++) {

			if (firmList.get(j).firmID != aFirm.firmID) {

				aFirm.competitorList.add(new FirmReference(
						firmList.get(j).firmID, firmList.get(j).quality,
						 firmList.get(j).qualityConcept, firmList.get(j).innovator,firmList.get(j).maxPendingQualityConcept));

			}

		}

		/* Update competitor list in FirmLocation object */
		for (int l = 0; l < aFirm.activeLocationList.size(); l++) {

			for (int j = 0; j < aFirm.activeLocationList.get(l).competitorList
					.size(); j++) {

				for (int k = 0; k < firmList.size(); k++) {

					if (aFirm.activeLocationList.get(l).competitorList
							.get(j).firmID == firmList.get(k).firmID) {

						aFirm.activeLocationList.get(l).competitorList
								.get(j).qualityConcept = firmList.get(k).qualityConcept;
						
						aFirm.activeLocationList.get(l).competitorList
						.get(j).maxPendingQualityConcept = firmList.get(k).maxPendingQualityConcept;
						break;

					}

				}

			}
		}

	}

}
	
	
	
}
