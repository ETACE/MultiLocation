import java.awt.Color;
import java.io.File;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Collections;
import uchicago.src.sim.engine.SimInit;
import uchicago.src.sim.engine.SimpleModel;
import uchicago.src.sim.gui.RectNetworkItem;
import uchicago.src.sim.network.DefaultEdge;
import uchicago.src.sim.network.Node;
import uchicago.src.sim.util.Random;


/**
 * This java code is related to the paper "R&D location strategies" authored by Luca Colombo, Herbert Dawid, and Philipp Harting 
 * 
 * 
 *Model.java is the central java class of the java model in which we define the structure of the model and the sequence of events.
 *For more information see the Readme file  
 */

public class Model extends SimpleModel{

	/*Declaration of model parameters*/
	
	/*++++Switches +++++*/
	
	/*Debugging and simulations*/
	
	public static boolean printDebug = false; /*Write data to text file*/
	public static boolean printDebugLocation= false; /*Write data to text file*/
	
	
	/*Model modes*/
	public static boolean marketEntryExit = true;  //+++++ Switch off for strategy analysis ++++
	public static boolean isMaxNumberLocation = true;  //+++++ Switch off for strategy analysis ++++

	public static boolean strategyExperiment = true;  //+++++ Switch on for strategy analysis ++++
	public static boolean modelModeTimeToMarket = false; //
	public static boolean strategyInnovator = true ;
	public static boolean monteCarloMedian = false;


	public static int industryScenario = 4; // Scenario 1: INS; Scenario 2: INW; Scenario 3: IMS; Scenario 4: IMW
	public static int maxNumberLoations = 10;
	public static boolean starRegionScenario = false; // if true, then one region has higher academic activity than other regions
	public static double acadamecicActivity = 1.0;
	public static boolean fixedNumLocations = true; // if this is true, there is  a flexible number of locations where at any point in time there is one empty location
	 
	/*++++ Parameters +++++*/
	
	/*Number of agents*/
	public static int numFirms = 6;   //Initial number of firms
	public static int numLocations = 5;  // Number of locations
	
	
	public static double marketEntryNu = 3.0;  
	public static int targetNumberOfFirms = 6;

	
	/*Model parameters*/
	public static int minNumFirms = 1; // Min num of firms
	public static int transitionPhase = 200;
	public static int timeHorizon = 100;  // length of NPV horizon in Location decision
	public static int numMontecarloRuns = 100; // Num of MC runs in location decision
	public static int numInitialLocations = 1;
	public static int parTimeToMarket = 40;
	public static int timeToMarket = 40;
	//public static int knowledgeDiffusionPeriod = 20 ; // before 010317: 40;
	public static int timeToImitation = 40;
	public double enteringCosts = 3.0; // Costs entering a location:
	public double exitingCosts = 0.0; // Costs exiting a location:
	public static double locationCosts = 0.015; // costs of operating in a location
	public static double reservationPrice = 2.0;  // a in the model description
	public static double productDifferentiation = 0.7; // Degree of product differentiation; gamma in the model description
	
	public static double maxQualityConceptProgress= 0.07; //0.01; // nu in the model description 
	public static double sdNoiseImitationKnowledge = 0.0; // sigma_a in the model description
	public static double sdNoiseImitationConcept= 0.0; // sigma_lambda in the model description 
	public static double discontFactor= 0.9967; // Discount factor used to compute the NPV
	public static double rateLocationDecision = 0.02;
	public static double sigma = 0.9;
	public static double householdBudget = 5; // Consumption budget of households
	public static double marketEntryHazardRate = 0.042;  // Hazard rate of market entry

	public double fractionEffectivityInnovators = 0.05;  // Determines the efficiency of innovators wrt imitation
	public double fractionEffectivityImitators = 0.05; // Determines the efficiency of imitators wrt innovation
	public static double marketEntryCosts = 0.5;  // Entry costs of firms in the market. 
	public static double marketExitHazardRate = marketEntryHazardRate; // hazard rate for exit
	public static double gammaMarketExit = 3.0;  // Intensity parameter for exit  
	
	/*Strategy parameters*/
	
	public static double strategyParameterEntry = 0.0;
	public static double strategyParameterExit = 0.0;
	public static double strategyParameterSwitch = 0.0;
	public static double strategyParameter = -0.1;
	public static double strategyParameterInnoDefault= 0.0;
	public static double strategyParameterInnoEntryDefault= 0.0;
	public static double strategyParameterInnoExitDefault= 0.0;
	public static double strategyParameterInnoSwitchDefault= 0.0;
	
	public static double strategyParameterImiDefault= 0.0;
	public static double strategyParameterImiEntryDefault= 0.0;
	public static double strategyParameterImiExitDefault= 0.0;
	public static double strategyParameterImiSwitchDefault= 0.0;
	
	public static int strategyFirmID = 6;
	
	
	/*Initial values firm variables - cannot be changed by the GUI */
	public double initialQualityConcept = 1.0;
	public double initialInnoEff = 0.04;//  
	public double initialImiEff = 0.05;// 
	
	
	/*Misc*/
	public static int numProcrssors =3;
	public static int iteration =0;
	
	
	long randomSeed;
	
	public static java.util.Random threadsafeRandomGenerator = new java.util.Random();;

	/*Aggregated variables*/
	int numActiveFirms = 5;
	int numInnovators;
	int numImitators;
	int overallFirmCounter;
	int numExitsNegativeOutput; // Count the firms that exit the market due to negative output
	double totalOutput; // total output
	double priceIndex ; // output weighted average price
	double qualityIndex; //output weighted average quality
	double totalProfit; // total profits including entering costs and inno as well as imitation expenditures
	double sumEquilProfits;
	double averageFirmLocations; //Average numbe rof locations
	double welfare; // Welfare measure
	double consumerSurplus; // Consumer surplus
	double sdOutput;
	double sdQuality;
	double sdProfit;
	double sdAverageLocation;
	double averageInnoProbability;
	double averageImiProbability;
	double sdInnoProbability;
	double sdImiProbability;
	static double maxQuality;
	static double minQuality;
	static double maxPublicKnowledge;
	static double minPublicKnowledge;
	
 
	double clusteringCoefficient;
	double weightedClusteringCoefficient;
	
	
	static boolean negativeOutput = false;
	public static boolean marketEntryHappend = false;
	public static boolean marketExitHappend = false;

	 public static int[][] networkMatrix ; 
	 
	

	/*Container for objects of class Firm*/
	static ArrayList<Firm> firmList = new ArrayList<Firm>();
	
	/*Container for objects of class Firm*/
	static ArrayList<Firm> exitedFirms = new ArrayList<Firm>();
	
	/*Container for  objects of class Location*/
	static ArrayList<Location> locationList = new ArrayList<Location>();
	
	ArrayList<Node> nodeList = new ArrayList<Node>();
	
	
	/*Used to write data to a firm specific text file (in debug mode)*/
	//static FirmWriter firmWriter;
	
	
	
	/*+++++++++++++Aux functions+++++++++++++*/
	/**
	 * Returns a pseudo-random number between min and max, inclusive.
	 *This function can be called from other classes
	 */
	public static int randomInt(int min, int max) {

		Random.createUniform();
	    return Random.uniform.nextIntFromTo(min, max);
	    				
	}
	
	/**
	 * Returns a uniformly distributed pseudo-random number between 0 and 1, inclusive.
	 *This function can be called from other classes
	 */
	public static double random(){
		
		double random;
	
			Random.createUniform();
			random = Random.uniform.nextDouble();
		
		 return random;
		
	}

	/**
	 * Returns a normally distributed pseudo-random number with mean and sd.
	 *This function can be called from other classes
	 */
	public static double randomNormal(double mean, double sd){
	
			Random.createNormal(mean, sd);	
			return   Random.normal.nextDouble();
		
	}

	/**
	 * This function returns the index of a firm with id within the firmList 
	 * */
	public static int returnFirmIndex(int id){
		
		int index= -1;
		
		for(int i=0; i < firmList.size(); i++){

			if(firmList.get(i).firmID==id){
				index = i;
				break;
			}			
		}
		return index;	
	}

	/**
	 * Returns the index from the firm with the maximum of variable varName
	 */
	public static int returnFirmIndexMax(String varName){
		
		int id = 0;
		
		double max = -99999999.0;
		for(int i=0; i < firmList.size(); i++){
		    Field field;
			try {
				field = firmList.get(i).getClass().getField(varName);
				
				Class clazzType = field.getType();
			    
			    if (clazzType.toString().equals("double")){
			    	
			    	double temp = Double.parseDouble(field.get(firmList.get(i)).toString());  
			    	if(temp>max){
			    		max = temp;
			    		id = i;
			    	}	
			    }else if( clazzType.toString().equals("int")){
				
			    	int temp = Integer.parseInt( field.get(firmList.get(i)).toString());  
			    	if(temp>max){
			    		max = temp;
			    		id = i;
			    	}
			    }
			}
			 catch (NoSuchFieldException e) {
				
				e.printStackTrace();
			} catch (SecurityException e) {
				
				e.printStackTrace();
			} catch (NumberFormatException e) {
				
				e.printStackTrace();
			} catch (IllegalArgumentException e) {
				
				e.printStackTrace();
			} catch (IllegalAccessException e) {
				
				e.printStackTrace();
			}
		    
		}

		 return(id); 
		
	}
	
	
	/**
	 * This function returns the median of a list of doubles
	 * */
	public static double returnMedianDoubleList(ArrayList<Double> list){

		double median = 0.0;
		
		Collections.sort(list);
		int mid = list.size()/2;
		median = (Double)list.get(mid);
		if (list.size()%2 == 0) {
		median = (median + (Double)list.get(mid-1))/2; 
		}
		
		return median;
		
	}
	
	
	/**
	 * This function returns the mean of a list of doubles
	 * */
		public static double returnMeanDoubleList(ArrayList<Double> list){

		double mean = 0.0;
		
		for(int i=0; i < list.size();i++){
			mean+= list.get(i);
		}
		
		if(list.size()>0)
			mean = mean / list.size();
		
		return mean;
		
	}

	/**
	 * Returns the index from the firm with the minimum of variable varName
	 */
	public static int returnFirmIndexMin(String varName){
		
		int id = 0;
		
		double min = 99999999.0;

		for(int i=0; i < firmList.size(); i++){

		    Field field;
			try {
				field = firmList.get(i).getClass().getField(varName);
				
				Class clazzType = field.getType();
			    
			    if (clazzType.toString().equals("double")){
			    	
			    	double temp = Double.parseDouble(field.get(firmList.get(i)).toString());  
			    	if(temp<min && temp!=0.0){
			    		min = temp;
			    		id = i;
			    	}
			    	    	
			    }else if( clazzType.toString().equals("int")){
				
			    	int temp = Integer.parseInt( field.get(firmList.get(i)).toString());  
			    	if(temp<min && temp!=0){
			    		min = temp;
			    		id = i;
			    	}
			    }
			}
			 catch (NoSuchFieldException e) {
				
				e.printStackTrace();
			} catch (SecurityException e) {
				
				e.printStackTrace();
			} catch (NumberFormatException e) {
				
				e.printStackTrace();
			} catch (IllegalArgumentException e) {
				
				e.printStackTrace();
			} catch (IllegalAccessException e) {
				
				e.printStackTrace();
			}
		    
		}
				 return   (id); 
		
	}

	
	/**
	 * Computes the percentage standard deviation for variable with name variableName contained in tempFirmList
	 * 
	 * */
	
	public static double standardDeviationFirmVariable(ArrayList<Firm> tempFirmList, String variableName){
		
		
		double sum = 0;
		double stdev = 0;

		double mean =  meanFirmVariable(tempFirmList,  variableName);

		for(int i=0; i < tempFirmList.size(); i++){
     
			try {
				Field field = tempFirmList.get(i).getClass().getField(variableName);
				
				Class clazzType = field.getType();
	  
			    if (clazzType.toString().equals("double")){

			    	double temp = Double.parseDouble(field.get(tempFirmList.get(i)).toString());  
			    	
				    sum+= Math.pow(temp-mean, 2);
	
			    }else if( clazzType.toString().equals("int")){
				
			    	int temp = Integer.parseInt( field.get(tempFirmList.get(i)).toString());  
			    	sum+= Math.pow(temp-mean, 2);
			    
			    }
			}
			 catch (NoSuchFieldException e) {
				
				e.printStackTrace();
			} catch (SecurityException e) {
				
				e.printStackTrace();
			}catch (IllegalArgumentException e){
				
				e.printStackTrace();
			} catch (IllegalAccessException e) {
				
				e.printStackTrace();
			}
		    
		}
		
		if(tempFirmList.size()>1)
			stdev = Math.pow(1.0/(tempFirmList.size()-1)*sum, 0.5);
		else
			stdev = 0.0; 
		
		
		if(Math.abs(mean)>1e-8)
			return stdev/mean;
		else 
			return 0.0;
	}
	

/**
 * Computes mean of variable variableName contained in tempFirmList
 */
	
public static double meanFirmVariable(ArrayList<Firm> tempFirmList, String variableName){

		double mean = 0;

		for(int i=0; i < tempFirmList.size(); i++){

		    Field field;
			try {
				field = tempFirmList.get(i).getClass().getField(variableName);
				
				Class clazzType = field.getType();
			    
			    if (clazzType.toString().equals("double")){
			    	
			    	double temp = Double.parseDouble(field.get(tempFirmList.get(i)).toString());  
			    	mean += temp; 

			    }else if( clazzType.toString().equals("int")){
				
			    	int temp = Integer.parseInt( field.get(tempFirmList.get(i)).toString());  
			    	mean += temp; 
			    
			    }
			}
			 catch (NoSuchFieldException e) {
				
				e.printStackTrace();
			} catch (SecurityException e) {
				
				e.printStackTrace();
			} catch (NumberFormatException e) {
				
				e.printStackTrace();
			} catch (IllegalArgumentException e) {
				
				e.printStackTrace();
			} catch (IllegalAccessException e) {
				
				e.printStackTrace();
			}
		    
		}
		
		mean = mean / tempFirmList.size();
		
		return mean ;
		
	}

/**
 * 
 * Add an empty location if there are none. This function is only executed if fixedNumLocations ==0
 * 
 */
	void addEmptyLocation(){
		
		//Check if there is an empty location
		int numEmpty = 0;
	
		
		for(int i=0; i < locationList.size(); i++){
			
			if(locationList.get(i).firmList.size()==0){
				
				numEmpty++;	
				
			}
			
		}
		
		System.out.println("Number locations before: "+numLocations);
		System.out.println("Number numEmpty: "+numEmpty);
		
		if(numEmpty==0){
			
			boolean found = false;
			int id=0;
			
			/*Find lowest non-used id*/
			for(int i =1; i < locationList.size()+1;i++){
				
				found = false;
				
				for(int j =0; j < locationList.size();j++){
					
				
					if(i ==  locationList.get(j).locationID){
						
						found = true;
						break;
						
					}				
				}
				
				if(!found){
					
					id =i;
					break;
				}
				
			}
			
			/*Otherwise use numLocations+1 as id*/
			if(found){
				
				id = locationList.size()+1;
				
			}

			if((isMaxNumberLocation && locationList.size()<maxNumberLoations)){
				
			
				
				Location aLocation = new Location(id,acadamecicActivity);
				locationList.add(aLocation);
				
				for(int i=0; i < firmList.size(); i++){
		
					firmList.get(i).passiveLocationList.add(new FirmLocation(firmList.get(i).firmID,aLocation.locationID ,aLocation.academicActivity));
					
				}
			
			}
		
		}else if(numEmpty>1 && locationList.size() > 5){
			
			for(int i=0; i < locationList.size();i++){
				
				if(locationList.get(i).firmList.size()==0){
					
					
					for(int j=0; j< firmList.size(); j++){
						
						for(int k =0; k < firmList.get(j).passiveLocationList.size();k++){
						
							if( firmList.get(j).passiveLocationList.get(k).locationID== locationList.get(i).locationID){
								
								
								 firmList.get(j).passiveLocationList.remove(k);
								 break;
					
							}
			
						}
					}
					
					
					locationList.remove(i);
					break;
					
				}
				
				
			}
			
		}

		numLocations = locationList.size();
		
		System.out.println("Number locations: "+numLocations);

	}
		

/**
 * 
 * This function updates the network matrix. Only used for illutration
 * 
 */
void updateNetworkMatrix(){
	
	networkMatrix = new int[firmList.size()][firmList.size()]; 
	
	for(int i=0; i < firmList.size();i++){
		
		for(int j=0; j < firmList.size();j++){
			
			Firm aFirm = firmList.get(i);
			
			if(i!=j){
				
				int counter = 0;
				
				for(int k=0; k < aFirm.activeLocationList.size();k++){
					
					for(int l=0; l < aFirm.activeLocationList.get(k).competitorList.size();l++){
					
						if(aFirm.activeLocationList.get(k).competitorList.get(l).firmID==firmList.get(j).firmID){
							
							counter++;
							break;
							
						}
	
					}
		
				}
				
				networkMatrix[i][j]= counter;
				
			}else
				networkMatrix[i][j]=1;
	
		}

	}

}
	

	/*++++++++++++ Generic repast functions to init and run simulations++++++++++++++*/
	
	
	
	/**
	 * 
	 * This method updates the firm network (only used in GUI mode)
	 * 
	 */
	public void updateFirmnetwork(){
		
		
		updateNetworkMatrix();
		
		/*Initialize the network*/
		
		nodeList.clear();
		
		for(int i=0; i < firmList.size();i++){
			
			firmList.get(i).computeLinks();
			
			RectNetworkItem drawable = new RectNetworkItem (i*10, i*10);
			FirmNode node = new FirmNode(firmList.get(i),drawable);
			
			node.setNodeLabel(Integer.toString(firmList.get(i).firmID));
			node.setLabelColor(Color.black);
			nodeList.add(node);
			
			
		}
		
		
		
		for(int i=0; i < nodeList.size();i++){
			
			Node node = (Node)nodeList.get(i);
			Node prevNode = null;
			
			
			for(int j =0; j < nodeList.size();j++){
				
					prevNode = (Node)nodeList.get(j);
					
					
					if(networkMatrix[i][j]>0){
			
						//EdgeFactory.createDrawableEdge(prevNode, node, Integer.toString(networkMatrix[i][j]),networkMatrix[i][j]);
						
						FirmEdge edge = new FirmEdge(prevNode,node);
						
						edge.setStrength(networkMatrix[i][j]);
						
						prevNode.addOutEdge(edge);
						node.addInEdge(edge);
				
					}
					
			}
			
		
		}
	
	}
	
	
	/**
	 * 
	 * This method analyzes the firm network; compute clustering coefficient
	 * 
	 */
	
	public void networkAnalysis(){
		
		
		double maxTriangles =0;
		int triangles =0;
		
		double valueClosedTriplets = 0;
		double totalValueTriplets = 0;
		
		double value1, value2;
		
		value1 = 0;
		value2=0;
		
		for(int i=0; i < nodeList.size(); i++){
			
			Node aNode = nodeList.get(i);
	
			for(int j=0; j < nodeList.size(); j++ ){
				
				if(aNode.hasEdgeTo(nodeList.get(j)) ){
					
					
					for(int k=0; k < nodeList.size(); k++ ){
						
						if(k!=i && k!= j &&  aNode.hasEdgeTo(nodeList.get(k))){
							
							for(int l =0; l < aNode.getOutEdges().size();l++){
								
								DefaultEdge edge = (DefaultEdge) aNode.getOutEdges().get(l);
								if(edge.getTo().equals(nodeList.get(k))){
									
									value1 = edge.getStrength();
									
								}else if (edge.getTo().equals(nodeList.get(j))){
									
									value2 = edge.getStrength();
									
								}
								
							}
							
							totalValueTriplets+= Math.sqrt(value1*value2);
							maxTriangles++;
							
							
							if((nodeList.get(j).hasEdgeTo(nodeList.get(k)))){
								
								valueClosedTriplets +=Math.sqrt(value1*value2);
								triangles++;
								
							}
							
							
						}
						
						
					}
							
					
				}
				
				
			}
			
		}
	
		
		
		clusteringCoefficient = 0.0;
		
		if(maxTriangles>0)
			clusteringCoefficient =   triangles/maxTriangles;
		
		weightedClusteringCoefficient=0;
		
		if(totalValueTriplets>0)
		  weightedClusteringCoefficient =  valueClosedTriplets /totalValueTriplets;
		
	}
	
	
	
	/**
	 * This method resets the firm lists in locations and the location arrays within the firm memory.
	 * */
	public void updateLocations(){
		
		
	/*  1. Global location list  */
		
		/*Clear location lists*/
		for(int i=0; i < locationList.size();i++){
			
			locationList.get(i).firmList.clear();	
			locationList.get(i).numberFirmsInLocation =0;
			
		}
		
		
		/*Fill location lists*/
		for(int i =0; i < firmList.size(); i++){
			
			Firm aFirm = (Firm) firmList.get(i); 
			
			for(int j=0; j< aFirm.activeLocationList.size();j++){
				
				for(int k =0; k < locationList.size();k++ ){
					
					if(locationList.get(k).locationID == aFirm.activeLocationList.get(j).locationID){
						
						locationList.get(k).firmList.add(aFirm);
						locationList.get(k).numberFirmsInLocation++;
						
					}
					
				}
					
			}
		}
		
		
		

		for (int i=0; i < locationList.size();i++){
			
			locationList.get(i).orderedIDList.clear();
			
			for(int j=0; j < locationList.get(i).firmList.size();j++){
				
				
				locationList.get(i).orderedIDList.add(locationList.get(i).firmList.get(j).firmID);
				
				
			}
			
			 Collections.sort(locationList.get(i).orderedIDList);
			
			
		}
			
		
		
		
		/* 2. Firm specific location lists */
		
		/*Update competitor list in FirmLocation object*/
		
		for(int i=0; i < firmList.size();i++){
			
			Firm aFirm = (Firm) firmList.get(i);
			
			for(int j=0; j < aFirm.activeLocationList.size();j++){
				
				FirmLocation aLocation = aFirm.activeLocationList.get(j);
				
				aLocation.competitorList.clear();
				
				for(int k=0; k < locationList.size();k++){
					
					if(locationList.get(k).locationID==aLocation.locationID){
						
						aLocation.orderedIDList.clear();
						
						for(int l =0; l < locationList.get(k).orderedIDList.size();l++){
							
							aLocation.orderedIDList .add(locationList.get(k).orderedIDList.get(l));
							
						}
						
						aLocation.type = 0;
						if(!fixedNumLocations)
							aLocation.toBeConsidered = false;
						else
							aLocation.toBeConsidered = true;
						
						
						for(int l=0; l <locationList.get(k).firmList.size();l++){
							
							if(locationList.get(k).firmList.get(l).firmID!= aFirm.firmID){
							for(int m =0; m < locationList.get(k).firmList.get(l).activeLocationList.size();m++ ){
								
								if(locationList.get(k).firmList.get(l).activeLocationList.get(m).locationID==locationList.get(k).locationID){
									
									aLocation.competitorList.add(new FirmReference(locationList.get(k).firmList.get(l).firmID,null, locationList.get(k).firmList.get(l).qualityConcept, locationList.get(k).firmList.get(l).innovator, locationList.get(k).firmList.get(l).maxPendingQualityConcept));
									break;
								}
							}
						
						}
						}
						
					}
					
					
				}
				
				aFirm.activeLocationList.get(j).countFirms();
				
			}
			
	
			for(int j=0; j < aFirm.passiveLocationList.size();j++){
				
				FirmLocation aLocation = aFirm.passiveLocationList.get(j);
				
				aLocation.competitorList.clear();
				
				for(int k=0; k < locationList.size();k++){
					
					if(locationList.get(k).locationID==aLocation.locationID){
						
						aLocation.orderedIDList.clear();
						
						for(int l =0; l < locationList.get(k).orderedIDList.size();l++){
							
							aLocation.orderedIDList .add(locationList.get(k).orderedIDList.get(l));
							
						}
						aLocation.type = 0;	
						
						if(!fixedNumLocations)
							aLocation.toBeConsidered = false;
						else
							aLocation.toBeConsidered = true;
						
						
						
						for(int l=0; l <locationList.get(k).firmList.size();l++){
							
							
							for(int m =0; m < locationList.get(k).firmList.get(l).activeLocationList.size();m++ ){
								
								if(locationList.get(k).firmList.get(l).activeLocationList.get(m).locationID==locationList.get(k).locationID){
									
									aLocation.competitorList.add(new FirmReference(locationList.get(k).firmList.get(l).firmID,null,locationList.get(k).firmList.get(l).qualityConcept,locationList.get(k).firmList.get(l).innovator, locationList.get(k).firmList.get(l).maxPendingQualityConcept));
									break;
								}
							}
						
						
						}
						
					}
					
					
				}
				
				aFirm.passiveLocationList.get(j).countFirms();
				
			}
			
			
		}
		
	
		
	}
	
	
	
	
	
	/**
	 * This method goes through the firm array and executes the location decisions of firms
	 * */
	public void executeFirmLocationDecision(){
		

		for(int i=0; i< firmList.size(); i++){
					
					Firm aFirm = (Firm) firmList.get(i);
					aFirm.firmSetup();
					aFirm.locationDecision();
				
		}
		
		updateLocations();
	
	}

	/**
	 * This method determines the profits, prices, and quantities of firms; if there are firms with negative quantities those are removed from the simulation
	 * */
	void executeFirmProfitDetermination(){
		
		negativeOutput =false;
		numExitsNegativeOutput = 0;
		
		/*Determine outputs*/
		for(int i=0; i< firmList.size(); i++){
			
			Firm aFirm = (Firm) firmList.get(i);
			aFirm.determineOutput();  
			
		}
		while(negativeOutput){
			
			numExitsNegativeOutput++;
			
			negativeOutput = false;
			double tempqual = 99999.0;
			int tempID = 0;
			
			/*Find firm with lowest quality*/
				
			for(int i=0; i< firmList.size(); i++){

				if(firmList.get(i).active && firmList.get(i).equilQuantity<tempqual){
						
						tempqual = firmList.get(i).equilQuantity;
						tempID = firmList.get(i).firmID;
	
				}				
			}
			
			/*Remove this firm*/
			
			for(int i=0; i< firmList.size(); i++){
				
				if(firmList.get(i).firmID==tempID){
					
					firmList.remove(i);
					break;
					
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
		
		
		/*Determine prices and profits*/
		
		/*Adjust competitor arrays*/
		for(int i=0; i< firmList.size(); i++){
			
			Firm aFirm = (Firm) firmList.get(i);
			aFirm.determinePriceProfits();
			
		}
		
		//System.out.println("Negative output: "+negativeOutput);
			
		}
		
	
	
	
	/**
	 * This method executes firms' innovation and imitation
	 * */
	void executeFirmInnovativeActivities(){
		
		
		/*Innovative and imitative activities of firm*/
		for(int i=0; i< firmList.size(); i++){
			
			Firm aFirm = (Firm) firmList.get(i);
		
					aFirm.innovationEffort();
					aFirm.imitationEffort();
					
		}
		
		/*Updating qualities for next iteration*/
		for(int i=0; i< firmList.size(); i++){
			
			Firm aFirm = (Firm) firmList.get(i);
			aFirm.updateQuality();
		
		}
		
		
		
		/*Updates the qualities in competitor list*/
		for(int i =0; i < firmList.size(); i++){
			
			for(int j=0; j < firmList.get(i).competitorList.size();j++){
				
				for(int k=0; k < firmList.size();k++){
					
					if(firmList.get(i).competitorList.get(j).firmID==firmList.get(k).firmID){
						
						
						firmList.get(i).competitorList.get(j).quality = firmList.get(k).quality;
						firmList.get(i).competitorList.get(j).qualityConcept=  firmList.get(k).qualityConcept;
						firmList.get(i).competitorList.get(j).maxPendingQualityConcept=  firmList.get(k).maxPendingQualityConcept;
						
					}
					
					
				}
	
			}
			
					
		}
		
		
	}
	
	
	
	

	/**
	 * This method executes firm market exists
	 * */
void firmMarketExit(){
		
		
		marketExitHappend = false;
		
		exitedFirms.clear();

		for(int i =0; i < firmList.size();i++){
			
			boolean exit = false;
			
			/* We have a minimum number of firms*/
			if(firmList.size()<= minNumFirms)
				break;
				
				
			/*The overall hazard rate is given by marketExitHazardRate; the probability for a single firm is depending on 
			the profit, where firms with ölower profit have a higher probability to exit*/
		
			double logit = 0.0, prob =0.0, random;
			double denominator = 0;
			
			for(int j=0; j < firmList.size();j++){
				
				
				denominator += Math.exp((-1)*gammaMarketExit*firmList.get(j).equilProfit);	
				
			}
			
			logit = Math.exp((-1)*gammaMarketExit*firmList.get(i).equilProfit)/denominator;
			
			prob = marketExitHazardRate * logit;
					

			random = random();
			
			
			if(random < prob)
				exit = true;
			
			/*If exit happened...*/
			if(exit){
		
				marketExitHappend = true;
				exitedFirms.add(firmList.get(i));
				firmList.remove(i);
				i--;
				numFirms--;

			}
			
			
		}
		
		/*Adjust the competitor list of remaining firms*/
		if(marketExitHappend){
			
			for(int i=0; i < firmList.size();i++){
			
				for(int k=0; k < exitedFirms.size();k++){
					for(int j=0; j < firmList.get(i).competitorList.size();j++){
					if(firmList.get(i).competitorList.get(j).firmID==exitedFirms.get(k).firmID){
						
						firmList.get(i).competitorList.remove(j);
						j--;
						
							
						}

					
				}
				
				
			}
			
			
		}
		
		
		}
		
	
	}
	
	
	/*
	 * This function determines the market entry of firms 
	 */
	void firmMarketEntry(boolean strategy, boolean innovator){
		
		double rand1, rand2;
		
		marketEntryHappend = false;

		/*If this function is not called by the exit function in case of replacing a firm with strategy profile, then the entry is purely random*/
		if(!strategy){
			rand1 = random();
			rand2 = random();
		}else{
		 
			rand1 = -1;
			if(innovator)
				rand2 = -1;
			else
				rand2 = 1;
		}
		
		double prob ;
		
		/*Exit rate is endogenous*/
		double marketEntryGamma = marketExitHazardRate / (Math.pow(householdBudget/targetNumberOfFirms,marketEntryNu));
		
		prob = Math.min(1,marketEntryGamma * (Math.pow(householdBudget/numActiveFirms,marketEntryNu)) );

		if(rand1 < prob){	
			
			marketEntryHappend =true;
	
			Firm aFirm;
		
			double meanQualityConcept = 0;
			double minQualityConcept = 0.0;
			int counter = 0;
			
			for(int i=0; i < firmList.size(); i++){
			
				if(firmList.get(i).qualityConcept>0){
					meanQualityConcept += firmList.get(i).qualityConcept;
					counter++;
					
					if(minQualityConcept==0.0 || firmList.get(i).qualityConcept< minQualityConcept){
						minQualityConcept = firmList.get(i).qualityConcept;
					}
				}
	
			}
			
			if(counter>0)
				meanQualityConcept = meanQualityConcept/counter;
			
			
			double qualityConcept = meanQualityConcept  + random()*(firmList.get(returnFirmIndexMax("qualityConcept")).qualityConcept - meanQualityConcept);
					
			
			/*Choose type of firm: here by chance  Changed 26.09.2016 -> strategic firm is not replaced immediately*/
			
				if(rand2 < 0.5){
						
					aFirm = new Firm(overallFirmCounter+1,initialInnoEff,fractionEffectivityInnovators*initialImiEff, qualityConcept,enteringCosts, exitingCosts, locationCosts,sigma,true );
					
					aFirm.activateStrategy(strategyParameterInnoEntryDefault, strategyParameterInnoExitDefault, strategyParameterInnoSwitchDefault, false);
				

							
				}else{
					aFirm = new Firm(overallFirmCounter+1,fractionEffectivityImitators*initialInnoEff,initialImiEff, minQualityConcept,enteringCosts, exitingCosts, locationCosts,sigma, false );
					aFirm.activateStrategy(strategyParameterImiEntryDefault, strategyParameterImiExitDefault, strategyParameterImiSwitchDefault, false);
				}
				
	
	
			aFirm.currentImiProb = aFirm.imitationProbability(1);
			aFirm.currentInnoProb = aFirm.innovationProbability();
			
			
			for(int j=0; j < locationList.size(); j++){
				
				aFirm.passiveLocationList.add(new FirmLocation(aFirm.firmID,locationList.get(j).locationID,locationList.get(j).academicActivity));
				
				
				
			}
			
		
			aFirm.firmSetup();
			aFirm.updateQuality();
			
			
			/*Add new firm to the firm list, create a file for storing the data*/
			firmList.add(aFirm);
			overallFirmCounter++;
			File file = new File("Firm"+(overallFirmCounter)+".txt");
			//firmWriter.writerList.add(file);
			numFirms++;
			
			
			/*Update the competitor list of other firms*/
			for(int i=0; i < firmList.size(); i++){
				
				if(firmList.get(i).firmID!=aFirm.firmID){
				
					firmList.get(i).competitorList.add(new FirmReference(aFirm.firmID,aFirm.quality,aFirm.qualityConcept, aFirm.innovator, 0.0));
				
				}
						
						
			}
			
			/*Create competitor list of new firm*/
			for(int j=0; j < firmList.size(); j++){
				
				if(firmList.get(j).firmID!=aFirm.firmID){
					
					Firm bFirm = (Firm) firmList.get(j);
					
					aFirm.competitorList.add(new FirmReference(bFirm.firmID,bFirm.quality,bFirm.qualityConcept, bFirm.innovator, 0.0));
					
				}
				
			}
			
			updateLocations();
	
		}
	
	}
	
	
	
	/*
	 * This determines the entry of the strategic firm in case of a strategy analysis
	 */
	void entryStrategicFirm(){
		
	
		Firm aFirm;
		
		double meanQualityConcept = 0;
		double minQualityConcept = 0.0;
		
		int counter = 0;
		
		// Determine minimum and average quality
		for(int i=0; i < firmList.size(); i++){
		
			if(firmList.get(i).qualityConcept>0){
				meanQualityConcept += firmList.get(i).qualityConcept;
				counter++;
				
				if(minQualityConcept==0.0 || firmList.get(i).qualityConcept< minQualityConcept){
					minQualityConcept = firmList.get(i).qualityConcept;
				}
			}

		}
		
		if(counter>0)
			meanQualityConcept = meanQualityConcept/counter;
		
		// Quality of strategic firm (if innovator is randomly drawn from the interval between maximum and average). The innovator enters with the publicly available one (minimum)
		double qualityConcept = meanQualityConcept  + random()*(firmList.get(returnFirmIndexMax("qualityConcept")).qualityConcept - meanQualityConcept);
				
		
			if(strategyInnovator){
					
				aFirm = new Firm(overallFirmCounter+1,initialInnoEff,fractionEffectivityInnovators*initialImiEff, qualityConcept,enteringCosts, exitingCosts, locationCosts,sigma,true );
				
				
				aFirm.activateStrategy(strategyParameterEntry, strategyParameterExit, strategyParameterSwitch, true);

						
			}else{
				aFirm = new Firm(overallFirmCounter+1,fractionEffectivityImitators*initialInnoEff,initialImiEff, minQualityConcept,enteringCosts, exitingCosts, locationCosts,sigma, false );
				aFirm.activateStrategy(strategyParameterEntry, strategyParameterExit, strategyParameterSwitch, true);
			}

		aFirm.currentImiProb = aFirm.imitationProbability(1);
		aFirm.currentInnoProb = aFirm.innovationProbability();
		
		
		for(int j=0; j < locationList.size(); j++){
			
			aFirm.passiveLocationList.add(new FirmLocation(aFirm.firmID,locationList.get(j).locationID,locationList.get(j).academicActivity));
	
		}
		
	
		aFirm.firmSetup();
		aFirm.updateQuality();
	
		
		/*Add new firm to the firm list, create a file for storing the data*/
		firmList.add(aFirm);
		overallFirmCounter++;
		File file = new File("Firm"+(overallFirmCounter)+".txt");
		//firmWriter.writerList.add(file);
		numFirms++;
		
		
		/*Update the competitor list of other firms*/
		for(int i=0; i < firmList.size(); i++){
			
			if(firmList.get(i).firmID!=aFirm.firmID){
			
				firmList.get(i).competitorList.add(new FirmReference(aFirm.firmID,aFirm.quality,aFirm.qualityConcept, aFirm.innovator, 0.0));
			
			}
					
					
		}
		
		/*Create competitor list of new firm*/
		for(int j=0; j < firmList.size(); j++){
			
			if(firmList.get(j).firmID!=aFirm.firmID){
				
				Firm bFirm = (Firm) firmList.get(j);
				
				aFirm.competitorList.add(new FirmReference(bFirm.firmID,bFirm.quality,bFirm.qualityConcept, bFirm.innovator, 0.0));
				
			}
			
		}
		
		updateLocations();

	}
	
	
	/*
	 * This method manages the market entry and exit
	 */
	
	void marketEntryExit(){
	
	//If market entry and exit is switched on: only potential entry after passing the transition phase
		if(marketEntryExit){
			if(iteration > transitionPhase){
				firmMarketExit();
				firmMarketEntry(false, false);
			}
	
		}else{
			
			marketEntryHappend = false;
			
			//In case of strategy experiment: Only at t=transitionPhase single entry of strategic firm
			if(iteration==transitionPhase){
	
				entryStrategicFirm();
				
				marketEntryHappend = true;
				
			}
			
		}
	
	}
	

	/*+++++++++++ Functions do manage agents actions (called within step())+++++++++++++++++++++++++++++++++++++++++++++++*/
	
	/**
	 *This method sets up the simulation
	 */
	public void setup(){
		
		super.setup();
		
	}
	
	
	/**
	 * This method initializes the simulation
	 * */
	public void buildModel(){
		
		iteration = 0;
		
		/*Distinguish two cases: Strategy experiment with fixed number of firms or industry simulation with market entry and exit*/
		if(strategyExperiment ){
			
			marketEntryExit = false;
		
			/*Scenario INS*/
			if(industryScenario==1){
				
				numFirms = 5;
				numInnovators = 3;	
				strategyInnovator =true;
			
				/*Scenario IMS*/	
			}else if(industryScenario==2){
				
				numFirms = 5;
				numInnovators = 4;
				strategyInnovator =false;
				
				/*Scenario INW*/	
			}else if(industryScenario==3){
				
				numFirms = 5;
				numInnovators = 1;
				strategyInnovator =true;
				
				/*Scenario IMW*/
			}else if(industryScenario==4){
				
				numFirms = 5;
				numInnovators = 2;
				strategyInnovator =false;
				
			
			}else{
				
				
				numInnovators = (int) Math.floor(numFirms/2.0);
				
			}
		
		
			/*Values for non-strategic firms*/
			strategyParameterImiEntryDefault = strategyParameterImiDefault;
			strategyParameterImiExitDefault = (-1)*strategyParameterImiDefault;
			strategyParameterImiSwitchDefault = strategyParameterImiDefault;
			
			strategyParameterInnoEntryDefault = strategyParameterInnoDefault;
			strategyParameterInnoExitDefault = (-1)*strategyParameterInnoDefault;
			strategyParameterInnoSwitchDefault = strategyParameterInnoDefault;
		
			
			/*values for strategic firm*/
			strategyParameterEntry = strategyParameter ;
			strategyParameterExit= (-1)*strategyParameter;
			
			if(strategyParameter<0.0)
				strategyParameterSwitch = 0.0;
			else
				strategyParameterSwitch = strategyParameter;
	
		}else{
			
			/*In industry analysis, the number of innovators is random*/
			
			numInnovators = randomInt(1,numFirms-1);
		
		}
	
		/*Delete all data bases*/
		if(printDebugLocation){
			
			 File directory = new File(System.getProperty("user.dir"));
			
			for(File f: directory.listFiles())
			    if(f.getName().startsWith("montecarlo"))
			        f.delete();
		}

// If modelModeTimeToMarket = true; then we have a certain period until an innovation can be introduced (time to market), otherwise the innovation is introduced next period
// Time to imitate, in contrast, is the period until the innovation becomes part of the standard product
		if(modelModeTimeToMarket){
			
			timeToImitation  = 0;
			timeToMarket = parTimeToMarket;
		}else{
			
			timeToImitation = parTimeToMarket;
			timeToMarket =1;
		}
	
	//Set random seeds
		Random.generateNewSeed();	
		randomSeed =  getRngSeed();
		System.out.println("Random seed:"+getRngSeed());
		
		/*This clears the lists from the previous run*/
		firmList.clear();
		locationList.clear();
		
		/*Create firms*/
		for(int i=0; i < numFirms; i++){
		
			Firm aFirm;
				if(i < numInnovators ){
						
					aFirm = new Firm(i+1,initialInnoEff,fractionEffectivityInnovators*initialImiEff, initialQualityConcept,enteringCosts, exitingCosts, locationCosts,sigma, true);
					
				}else{
					aFirm = new Firm(i+1,fractionEffectivityImitators*initialInnoEff,initialImiEff, initialQualityConcept, enteringCosts, exitingCosts,locationCosts, sigma, false );
				}
			
			firmList.add(aFirm);
			
		}

			/*If we conduct a strategy analysis*/
			if(strategyExperiment ){
				
				/*Go through the firm list and find the firm that is the firm considered*/
				for(int i=0; i < firmList.size(); i++){
					
					if(firmList.get(i).innovator){
					
						firmList.get(i).activateStrategy(strategyParameterInnoEntryDefault, strategyParameterInnoExitDefault, strategyParameterInnoSwitchDefault, false);
						
					}else{
				
						firmList.get(i).activateStrategy(strategyParameterImiEntryDefault, strategyParameterImiExitDefault, strategyParameterImiSwitchDefault, false);
						
					}
				
				}
		
		}
				
		/*Create Locations*/
		for(int i=0; i < numLocations; i++){
			
			double factor = 1.0;
			
			if(starRegionScenario && i==0)
				factor = 1.2;
			
			Location aLocation = new Location(i+1,factor*acadamecicActivity);
			locationList.add(aLocation);
			
		}		
			

		/*Distributing the firms among locations
		 * we assume that imitators are located in 
		 * those locations where innovators are already located */
	
		ArrayList<Integer> locIndexInnovators = new ArrayList<Integer> ();
		int locIndex = 0;
	
		
		/*First distribute the innovators */
		for(int i=0; i < firmList.size(); i++){
			Firm aFirm = (Firm) firmList.get(i);
			
			if(aFirm.innovator){
		
				aFirm.qualityConcept = aFirm.qualityConcept;
			
				aFirm.activeLocationList.add(new FirmLocation(aFirm.firmID,locationList.get(locIndex%numLocations).locationID , locationList.get(locIndex%numLocations).academicActivity));
				locationList.get(locIndex%numLocations).firmList.add(aFirm);
			
				for(int j=0; j < locationList.size(); j++){
					
					if((locIndex%numLocations) !=j)
					aFirm.passiveLocationList.add(new FirmLocation(aFirm.firmID,locationList.get(j).locationID ,locationList.get(j).academicActivity));
					
				}
		
				/*Compute initial quality*/
				
				aFirm.firmSetup();
				aFirm.updateQuality();
				
				locIndexInnovators.add(locIndex%numLocations);
				locIndex++;
				
				}
	
			}
				
		/*Second: distribute imitators*/		
		for(int i=0; i < firmList.size(); i++){
			Firm aFirm = (Firm) firmList.get(i);
			
				if(!aFirm.innovator){
	
				aFirm.qualityConcept = aFirm.qualityConcept;
				
				if(locIndexInnovators.size()>0)
					locIndex = locIndexInnovators.get(Model.randomInt(0, locIndexInnovators.size()-1));
				else
					locIndex = Model.randomInt(0, numLocations-1);
			
				aFirm.activeLocationList.add(new FirmLocation(aFirm.firmID,locationList.get(locIndex).locationID , locationList.get(locIndex).academicActivity));
				locationList.get(locIndex).firmList.add(aFirm);
		
			for(int j=0; j < locationList.size(); j++){
				
				if(locIndex !=j)
				aFirm.passiveLocationList.add(new FirmLocation(aFirm.firmID,locationList.get(j).locationID,locationList.get(j).academicActivity));
				
			}
			
			
			/*Compute initial quality*/
			
			aFirm.firmSetup();
			aFirm.updateQuality();
			
			}
			
		
			
		}
			
		overallFirmCounter =  firmList.size();
	
		/*Initialize the competitor lists of each firm*/		
		for(int i=0; i < firmList.size(); i++){
			
			Firm aFirm = (Firm) firmList.get(i);
			
			for(int j=0; j < numFirms; j++){
				
				if(i!=j){
					
					Firm bFirm = (Firm) firmList.get(j);
					
					aFirm.competitorList.add(new FirmReference(bFirm.firmID,initialQualityConcept, 0.0, bFirm.innovator, 0.0));  
					
				}
				
			}
					
					
		}
	
	//Update location information in each firm
		updateLocations();

	totalProfit =0;
	priceIndex = 0;
	sumEquilProfits=0;

	/*Setup the firm network*/
	updateFirmnetwork();
		
}
	
	
	/**
	 * This method defines the time line; it calls the methods that are sequentially executed in each iteration
	 * */
	public void step(){
		
		iteration = (int) getTickCount();
		
		System.out.println("Iteration "+iteration);

		if(!fixedNumLocations ){
			addEmptyLocation();
		}
		executeFirmLocationDecision();
		executeFirmProfitDetermination();
		executeFirmInnovativeActivities();
		marketEntryExit();
		updateFirmnetwork();
	
		result();
		
		System.out.println();
	
	}
	

	/**
	 * 
	 * This method computes and reports results at the end of each iteration
	 * */
	public void result(){
		
		numActiveFirms = 0;
		numInnovators =0;
		numImitators = 0;
		totalOutput = 0;
		priceIndex = 0;
		qualityIndex = 0;
		averageFirmLocations = 0;
		totalProfit =0;
		sumEquilProfits = 0;
		maxQuality = 0;
		minQuality = 999999.0;
		
		
		/*Go through firms' memory and aggregate firm variables*/
		for(int i=0; i < firmList.size(); i++){
			
			Firm aFirm = (Firm) firmList.get(i);
			
			
			if(aFirm.active){
			
				numActiveFirms++;
				
				if(aFirm.innovator)
					numInnovators++;
				else
					numImitators++;
				
				totalOutput+= aFirm.equilQuantity;
				priceIndex += aFirm.equilPrice;
				qualityIndex += aFirm.quality;
				averageFirmLocations += aFirm.numLocationsActive;
				totalProfit+= aFirm.totalProfit;
				sumEquilProfits += aFirm.equilProfit;
				
				maxQuality = Math.max(maxQuality, aFirm.quality);
				minQuality = Math.min(minQuality, aFirm.quality);
				
				
			
			}
		}
		
		priceIndex = priceIndex/numActiveFirms;
		qualityIndex = qualityIndex/numActiveFirms;
		averageFirmLocations = averageFirmLocations/numActiveFirms;
		
	
		
		
		double w1=0;
		double w2=0;
		double w3=0;
		
		for(int i=0; i< firmList.size();i++){
			
			if(firmList.get(i).equilQuantity>0){
				
				w1 += reservationPrice * firmList.get(i).equilQuantity;
				w2 += 0.5* Math.pow(firmList.get(i).equilQuantity, 2)/Math.pow(firmList.get(i).quality, 2);
				for(int j=0; j< firmList.size();j++){
					if(i!=j){
						w3 += productDifferentiation*(firmList.get(i).equilQuantity*firmList.get(j).equilQuantity)/(firmList.get(i).quality*firmList.get(j).quality);
					}
				}
			}
			
		}
		
		consumerSurplus = w1 - w2 - w3;
		welfare = consumerSurplus + totalProfit;
		
		sdOutput = standardDeviationFirmVariable(firmList, "equilQuantity"); 
		sdQuality = standardDeviationFirmVariable(firmList, "equilQuantity");
		sdAverageLocation =  standardDeviationFirmVariable(firmList, "numLocationsActive");
		averageInnoProbability = meanFirmVariable(firmList, "currentInnoProb");
		averageImiProbability = meanFirmVariable(firmList, "currentImiProb");
		sdInnoProbability = standardDeviationFirmVariable(firmList, "currentInnoProb");
		sdImiProbability = standardDeviationFirmVariable(firmList, "currentImiProb");
		sdProfit = standardDeviationFirmVariable(firmList, "equilProfit");
		
		networkAnalysis();
	}
	
	
	
	
	
	/**
	 * *** This is the main method ***
	 * */
	public static void main(String[] args){
		
		SimInit init = new SimInit();
		Model m = new Model();
		init.loadModel(m, null, false);
		
	}

}