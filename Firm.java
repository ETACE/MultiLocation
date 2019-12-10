import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Collections;

import uchicago.src.sim.gui.Drawable;
import uchicago.src.sim.gui.SimGraphics;

/**
*This class represents the agent type firm, that is interacting with other
* firms on a central goods market but can carry out its R&D in different locations 
*/

public class Firm implements Cloneable, Drawable {
    

	
	
	/****************Cloning   to avoid problems in deep copying************************/
	
	public Firm clone()  {
        try {
			return (Firm) super.clone();
		} catch (CloneNotSupportedException e) {
			
			e.printStackTrace();
			return null;
		}
		
        
    }
	
	
	/*Declaration of internal memory*/
	
		/*Constant firm parameters*/
		
		public int firmID;
		public int x_coord;
		public int y_coord;
		boolean active;
		boolean innovator = false;
		boolean firmStratgy = false;
		public double efficiencyInnovation;
		public double efficiencyImitation;
	
		/*Firm variables*/
		public int numLocationsActive;
		public int age;
		public int directLinksInnovators;
		public int directLinksImitators;
		public double qualityConcept;
		public double tempQualityConcept;
		public double maxPendingQualityConcept;
		public double quality;
		public double equilQuantity;
		public double equilPrice;
		public double equilProfit;
		public double marketProfit;
		public double totalProfit;
		public double strategyParameterEntry;
		public double strategyParameterExit;
		public double strategyParameterSwitch;
		public double testedNPVBefore;
		public double testedNPVAfter;
		public double movingCosts;
		public double locationCosts;
		public double enteringCosts;
		public double exitingCosts;
		public double currentInnoProb;
		public double currentImiProb;
		public double cumProfit;
		public double sigma;
		public double lambda;
		

		/*This holds the pending innovation processes*/
		ArrayList<InnovationProcess> pendingInnovationsList = new ArrayList<InnovationProcess>();

		/*Container of objects location: active locations*/
		protected ArrayList<FirmLocation> activeLocationList = new ArrayList<FirmLocation>();
		
		/*Container of objects location: passive locations*/
		protected ArrayList<FirmLocation> passiveLocationList = new ArrayList<FirmLocation>();
		
		/*Container of objects competitors*/
		protected ArrayList<FirmReference> competitorList = new ArrayList<FirmReference>();
		
		/*Container of objects entry*/
		protected ArrayList<EntryCharacteristics> entryList= new ArrayList<EntryCharacteristics>();
		
		/*Container of objects exits*/
		protected ArrayList<ExitCharacteristics> exitList= new ArrayList<ExitCharacteristics>();
	
	
	
	
	/*Constructor of class Firm*/
	
	Firm(int id, double effInno, double effImi,  double initQualityConcept,  double entCosts, double exCosts, double locCosts, double sig){
		
		 firmID = id;
		 efficiencyInnovation = effInno;
		 efficiencyImitation = effImi;
		 qualityConcept = initQualityConcept;
		 tempQualityConcept = initQualityConcept;
		 active = true;		
		 currentInnoProb  = effInno;
		 currentImiProb = effImi;
		 age = 0;
		 strategyParameterEntry = 0.0;
		 strategyParameterExit = 0.0;
		 strategyParameterSwitch = 0.0;
		 sigma = sig;
		 locationCosts = locCosts;
		 enteringCosts = entCosts;
		 exitingCosts= exCosts; 
		 cumProfit = 0;
		
			
		}
	
	
	
	
	Firm(int id, double effInno, double effImi, double initQualityConcept,   double entCosts, double exCosts, double locCosts, double sig, boolean inno){
		
		firmID = id;
		efficiencyInnovation = effInno;
		efficiencyImitation = effImi;
		active = true;
		innovator = inno;
		currentInnoProb  = effInno;
		currentImiProb = effImi;
		sigma = sig;
		qualityConcept = initQualityConcept;
		tempQualityConcept = initQualityConcept;
		locationCosts = locCosts;
		enteringCosts = entCosts;
		exitingCosts= exCosts;
		strategyParameterEntry = 0;
		strategyParameterExit = 0;
		strategyParameterSwitch = 0;
		cumProfit = 0;
		age = 0;
	}
	
	
	/*++++++++++++++++++++++++Methods+++++++++++++++++++++++++++++++++++++++++++++*/
	
	
	
	/**
	 * This function is used to activate a particular firm for the strategy analysis
	 * */
	void activateStrategy( double stratEntry, double stratExit, double stratSwitch , boolean strategicFirm ){
		
		firmStratgy = strategicFirm;
		
		
		strategyParameterEntry = stratEntry;
		strategyParameterExit = stratExit;
		strategyParameterSwitch = stratSwitch;
		
		
	}
	
	
	/**
	 * This function is executed first in the time line of a firm. Here we determine the maximum pending quality concept
	 * */
	void firmSetup(){
		
		
		active = true;
		age++;

		movingCosts = 0;
		
		maxPendingQualityConcept = 0.0;
		
		for(int i=0; i < pendingInnovationsList.size(); i++)
		
			maxPendingQualityConcept = Math.max(maxPendingQualityConcept, pendingInnovationsList.get(i).qualityConcect ) ;
		
	}
	
	
	
	/**
	 * This method calls the entry and exit methods;
	 */
void locationDecision(){

		/*Location decision with rate Model.rateLocationDecision */
		
		/*If established firm: location decision (entry, exit switching with low frequency acc. to hazard rate)*/
		
		/*Clear entry and exit list*/
		entryList.clear();
		exitList.clear();

		//This if-else is used to distinguish between the relocation decision of incumbent firms and entrants
		//If incubent: decide about possible relocation
		if(activeLocationList.size()>0){
		
			/*In the transition phase, firms do the location decision more often!*/
			if((Model.random() <   Model.rateLocationDecision && Model.iteration>Model.transitionPhase) || ((Model.random() <   0.1 && Model.iteration<=Model.transitionPhase))){
				
		
				
				if(!Model.fixedNumLocations){
					
					// This procedure is used, in case of a flexible number of locations, to reduce the computational burden. If locations have the same profile, they are not explicitly reconsidered in the MC later on

					Collections.shuffle(activeLocationList);
					Collections.shuffle(passiveLocationList);
					
					/*Check which locations in the active location list have the same profile*/
					int type = 1;
					for(int i=0; i < activeLocationList.size(); i++ ){
						
						if(activeLocationList.get(i).type==0){
							
							activeLocationList.get(i).type = type;
							
							activeLocationList.get(i).toBeConsidered= true;
							
							for(int j=i+1; j < activeLocationList.size(); j++ ){
								
								if(activeLocationList.get(i).orderedIDList.equals(activeLocationList.get(j).orderedIDList) && activeLocationList.get(j).type==0){
									
									activeLocationList.get(j).type= type;			
									
								}
				
							}
							
							type++;
						}
	
					}
	
					/*Check which locations in the passive location list have the same profile*/
					type = 1;
					
					for(int i=0; i < passiveLocationList.size()-1; i++ ){
						
						if(passiveLocationList.get(i).type==0){
							
							passiveLocationList.get(i).type = type;
							passiveLocationList.get(i).toBeConsidered= true;
							
							
							for(int j=i+1; j < passiveLocationList.size(); j++ ){
								
								if(passiveLocationList.get(i).orderedIDList.equals(passiveLocationList.get(j).orderedIDList) && passiveLocationList.get(j).type==0){
									
									passiveLocationList.get(j).type= type;
							
								}

							}
			
							type++;
						}
						
					}
					
				}
	
				/***************Monte Carlo analysis of NPVs *****************************/
				
				
				ArrayList<Strategy> strategies = new ArrayList<Strategy>();
				double npvBaseline, tempNPV;
				
				// Determine entering and exiting costs
				double effectiveExitingCosts  = exitingCosts * (1-Math.pow(1-Model.discontFactor,Model.timeHorizon+1))/(Model.discontFactor);
				double effectiveEnteringCosts  = enteringCosts * (1-Math.pow(1-Model.discontFactor,Model.timeHorizon+1))/(Model.discontFactor);

				/*This creates an instance of the MonteCarlo simulator*/
				MonteCarloGeneric montecarloSimulator = new MonteCarloGeneric(this, competitorList);
				
				montecarloSimulator.setup();
			
				/*Baseline without entry or exit*/
				npvBaseline = montecarloSimulator.monteCarloBase();
				
				
				/*Strategy 1: Exit only*/

				/*Exit only if the firm is in at least in two locations*/
				if(activeLocationList.size()>1){

					for (int i=0; i< activeLocationList.size();i++){
						
						if(activeLocationList.get(i).toBeConsidered){
		
							// Compute NPV of exiting this particular region
							tempNPV = montecarloSimulator.monteCarloExit( activeLocationList.get(i).locationID);
					
							if(Math.abs(npvBaseline)==0){
								
								System.out.println("Stopp because npv is 0 in Firm.locationDecision  (firmID)"+firmID);
								//System.exit(0);
								
							}
							
							//If NPV is higher than threshold -> add to considered strategies
							if(((tempNPV - npvBaseline) - effectiveExitingCosts)/(npvBaseline) > strategyParameterExit ){
								
								strategies.add(new Strategy(activeLocationList.get(i).locationID, false,true, ((tempNPV - npvBaseline) - effectiveExitingCosts)/(npvBaseline) , ((tempNPV - npvBaseline) - effectiveExitingCosts)/(npvBaseline) - strategyParameterExit , npvBaseline,tempNPV));	
					
							}
							
							//Store information of this strategy
							exitList.add(new ExitCharacteristics( activeLocationList.get(i).locationID, activeLocationList.get(i).competitorList.size(),0, 0,((tempNPV - npvBaseline) - effectiveExitingCosts)/(npvBaseline) , ((tempNPV - npvBaseline) - effectiveExitingCosts)/(npvBaseline) - strategyParameterExit, npvBaseline, tempNPV));
						}	
						
					}
					
				}
				
				
				/*Strategy 2:  Entry and exit  (Switching)*/

				for (int i=0; i< activeLocationList.size();i++){
					
					if(activeLocationList.get(i).toBeConsidered){
					
						for(int j=0; j < passiveLocationList.size();j++){
							
							if(passiveLocationList.get(j).toBeConsidered){
							
						
								montecarloSimulator.setPendingEntry(true);
								montecarloSimulator.setEntryID( passiveLocationList.get(j).locationID);
							
								//Compute NPV
								tempNPV = montecarloSimulator.monteCarloEntryExit(activeLocationList.get(i).locationID, passiveLocationList.get(j).locationID);
			
								if(Math.abs(npvBaseline)==0){
									
									System.out.println("Stopp because npv is 0 in Firm.locationDecision  (firmID)"+firmID);
									//System.exit(0);
									
								}
								
								
								//If NPV is higher than threshold -> add to considered strategies	
								if(((tempNPV - npvBaseline) - effectiveEnteringCosts - effectiveExitingCosts)/npvBaseline > strategyParameterSwitch ){
						
									strategies.add(new Strategy(passiveLocationList.get(j).locationID,activeLocationList.get(i).locationID, true,true, ((tempNPV - npvBaseline) - effectiveEnteringCosts - effectiveExitingCosts)/npvBaseline , ((tempNPV - npvBaseline) - effectiveEnteringCosts - effectiveExitingCosts)/npvBaseline - strategyParameterSwitch, npvBaseline, tempNPV));	
				
								}
								
								//Store information of this strategy
								exitList.add(new ExitCharacteristics( activeLocationList.get(i).locationID,  activeLocationList.get(i).competitorList.size(),0,1, ((tempNPV - npvBaseline)- effectiveEnteringCosts - effectiveExitingCosts)/(npvBaseline) , ((tempNPV - npvBaseline)- effectiveEnteringCosts - effectiveExitingCosts)/(npvBaseline) - strategyParameterSwitch, npvBaseline, tempNPV));
								entryList.add(new EntryCharacteristics(passiveLocationList.get(j).locationID,passiveLocationList.get(j).competitorList.size(),0,0, 1, ((tempNPV - npvBaseline) - effectiveEnteringCosts - effectiveExitingCosts)/npvBaseline,((tempNPV - npvBaseline) - effectiveEnteringCosts - effectiveExitingCosts)/npvBaseline  - strategyParameterSwitch, npvBaseline, tempNPV));
							}
						}
					}
					
				}
					
				
				/*Strategy 3:  Entry */
				
				for (int i=0; i< passiveLocationList.size();i++){
					
					
					if(passiveLocationList.get(i).toBeConsidered){
					
						montecarloSimulator.setPendingEntry(true);
						montecarloSimulator.setEntryID( passiveLocationList.get(i).locationID);
					
						//Compute NPV
						tempNPV = montecarloSimulator.monteCarloEntry(passiveLocationList.get(i).locationID);
					
						
						if(Math.abs(npvBaseline)==0){
							
							System.out.println("Stopp because npv is 0 in Firm.locationDecision  (firmID)"+firmID);
							//System.exit(0);
							
						}

						//If NPV is higher than threshold -> add to considered strategies	
						if((tempNPV - npvBaseline - effectiveEnteringCosts)/npvBaseline > strategyParameterEntry ){
				
							strategies.add(new Strategy(passiveLocationList.get(i).locationID, true,false, (tempNPV - npvBaseline - effectiveEnteringCosts)/npvBaseline , (tempNPV - npvBaseline - effectiveEnteringCosts)/npvBaseline - strategyParameterEntry, npvBaseline, tempNPV));	
							
						}
						
						//Store information of this strategy
						entryList.add(new EntryCharacteristics(passiveLocationList.get(i).locationID, passiveLocationList.get(i).competitorList.size(),0,0, 0, ((tempNPV - npvBaseline) - effectiveEnteringCosts )/npvBaseline,((tempNPV - npvBaseline) - effectiveEnteringCosts )/npvBaseline  - strategyParameterEntry,npvBaseline, tempNPV));
					}
					
				}
				
				
	// Here we take the actual decision. If there are strategies whose NPV is larger than the threshold, we choose the one with largest NPV
				
				boolean locationChoice = false;
				
				if(strategies.size()>0){
					
					
					/*Shuffle list:*/
					
					Collections.shuffle(strategies);
				
					int index = 0;
					
					/*Choose best strategy*/
					for(int i=0; i < strategies.size(); i++){
						
						if(strategies.get(i).npvIncKappa>= strategies.get(index).npvIncKappa){
							
							index = i;
							locationChoice = true;
							
						}
					
					}
		
					/*If there has been a location choice...*/	
					
					// If exit is involved: change active and passive location list
					if(locationChoice && strategies.get(index).exit){
						
						for(int i=0; i < activeLocationList.size();i++){
							
							if(locationChoice && strategies.get(index).exit && activeLocationList.get(i).locationID==strategies.get(index).locationIDExit){
								
								
								for(int j=0; j < exitList.size(); j++){
									
									if(exitList.get(j).locationID == strategies.get(index).locationIDExit &&  Math.abs( strategies.get(index).npvIncKappa - exitList.get(j).npvInclKappa)<1e-8){
								
										exitList.get(j).exit = 1;
								
										break;
									}
								
								}
								
								passiveLocationList.add(activeLocationList.get(i));
								activeLocationList.remove(i);
								
								movingCosts += effectiveExitingCosts;
								
								break;
										
							}
						
						}
					
					}
							
					// If exit is involved: change active and passive location list
					if(locationChoice && strategies.get(index).entry){
			
							
							for(int i=0; i < passiveLocationList.size();i++){
								
								if( passiveLocationList.get(i).locationID==strategies.get(index).locationIDEntry){
						
									for(int j=0; j < entryList.size(); j++){
										
										if(entryList.get(j).locationID == strategies.get(index).locationIDEntry &&  Math.abs( strategies.get(index).npvIncKappa - entryList.get(j).npvInclKappa)<1e-8){
											
											
											entryList.get(j).entry = 1;
									
											break;
										}
								
									}
									
									activeLocationList.add(passiveLocationList.get(i));
									passiveLocationList.remove(i);
									movingCosts += effectiveEnteringCosts;
									break;		
									
							}
							
						}
					
					}
			
				}
					
				/*Update num firm locations*/
				montecarloSimulator.destructMontecarloSimulator();
				montecarloSimulator = null;

			}
		
		}else{
			/*Else:  If the firm is new to market; only decide in which location to enter*/
			
			ArrayList<Strategy> strategies = new ArrayList<Strategy>();
			
			MonteCarloGeneric montecarloSimulator = new MonteCarloGeneric(this, competitorList);
			double tempNPV;
			montecarloSimulator.setup();
			
			for (int i=0; i< passiveLocationList.size();i++){
	
				/*This creates an instance of the MonteCarlo simulator*/

				montecarloSimulator.setPendingEntry(true);
				montecarloSimulator.setEntryID( passiveLocationList.get(i).locationID);
				//Compute npv
				tempNPV = montecarloSimulator.monteCarloEntry(passiveLocationList.get(i).locationID);

				strategies.add(new Strategy(passiveLocationList.get(i).locationID, true,false, tempNPV, tempNPV, 0.0, tempNPV));	
		
			}
			
			/*Update num firm locations*/
		
			montecarloSimulator.destructMontecarloSimulator();
			montecarloSimulator = null;

			int index = 0;
			
			/*Choose best strategy*/
			for(int i=0; i < strategies.size(); i++){
				
				if(strategies.get(i).npv>= strategies.get(index).npv){
					
					index = i;
			
				}
			}
			
			for(int i=0; i < passiveLocationList.size();i++){
				
				if( passiveLocationList.get(i).locationID==strategies.get(index).locationIDEntry){
					
					entryList.add(new EntryCharacteristics(passiveLocationList.get(i).locationID,passiveLocationList.get(i).competitorList.size(),1,1,0, strategies.get(index).npv,strategies.get(index).npvIncKappa, strategies.get(index).npvBaseline, strategies.get(index).npvGross));
					
					activeLocationList.add(passiveLocationList.get(i));
					passiveLocationList.remove(i);
					break;		
					
				}
			
			}
			
				
			/*Write the no-entries into the entry list*/
			for(int i=0; i < passiveLocationList.size();i++){
							
				/*Do not add the location the firm has just left in case switching*/
				
					entryList.add(new EntryCharacteristics(passiveLocationList.get(i).locationID,passiveLocationList.get(i).competitorList.size(),0,1,0, strategies.get(index).npv,strategies.get(index).npvIncKappa, strategies.get(index).npvBaseline, strategies.get(index).npvGross));
				
			}
			
		}
				
}
	
	/**
	 * This method returns the current imitation probability
	 */
	double imitationProbability(int numCompetitorsLocatoion){
	
		double probability = 0.0;
		
		
		if(numCompetitorsLocatoion>0){
			
			probability = efficiencyImitation;
		}

		return probability;
		
	}
		
		
		
	/**
	 * This method returns the current innovation probability
	 */
	double innovationProbability(){

		
		double probability = 0.0;
		
		if(activeLocationList.size()>0){
		
			lambda = 0.0;
			
			for(int i=0; i < activeLocationList.size(); i++){
				
				lambda += Math.pow(activeLocationList.get(i).academicActivity,sigma);
				
				
			}
			
			lambda = Math.pow(lambda, 1.0/sigma)/(1.0*activeLocationList.size() ); 
			
			probability = efficiencyInnovation * lambda;
		}
		
		return probability;
		
	}
		
		
	/**
	 * This method computes the current output
	 */
	void determineOutput(){
		
		/*Auxiliary calculation: compute the sum of quality differentials*/
		
		double sum1 = 0;
		for(int i=0; i <competitorList.size(); i++ ){
	
			sum1 += (quality - competitorList.get(i).quality);
			
		}
		
		/*Compute output quantity*/
		 equilQuantity = Model.reservationPrice * quality * ((2 - Model.productDifferentiation)*quality + Model.productDifferentiation* sum1)/((2-Model.productDifferentiation)*(2+Model.productDifferentiation*(competitorList.size())));
		
		 if(equilQuantity<0){
				
			 Model.negativeOutput = true;
		
		 }
		
	}
		
		
	/**
	 * This method computes the current price and profit
	 */	
		void determinePriceProfits(){
			
			double lambda;
			double sum0 = 0;
			double sum1 = equilQuantity;
			double sum2 = Math.pow(equilQuantity/quality, 2);
			double sum3 = 0;
			
			for(int i=0; i < competitorList.size();i++){
					
				sum0 += (quality - competitorList.get(i).quality);
				sum1 += competitorList.get(i).output;
				sum2 +=  Math.pow(competitorList.get(i).output/competitorList.get(i).quality, 2);
				sum3 += (competitorList.get(i).output*equilQuantity)/(quality*competitorList.get(i).quality);
			}
			
			
			for(int i=0; i < competitorList.size();i++){
				
				sum3 += (competitorList.get(i).output*equilQuantity)/(quality*competitorList.get(i).quality);
				
				for(int j=0; j < competitorList.size();j++){
					
					if(i!=j){
						
						sum3 += (competitorList.get(i).output*competitorList.get(j).output)/(competitorList.get(i).quality*competitorList.get(j).quality);
						
						
					}
			
				}

			}
			
			lambda = Model.householdBudget/(Model.reservationPrice*sum1  - sum2 - Model.productDifferentiation*sum3); 
			
			//Price 
			equilPrice = (lambda* equilQuantity)/(Math.pow(quality, 2)); 

			/*Equ. profit*/
			
			marketProfit = (Math.pow(Model.reservationPrice, 2)*lambda)/(Math.pow(2-Model.productDifferentiation, 2)*Math.pow(2+Model.productDifferentiation*competitorList.size(), 2))*Math.pow(((2-Model.productDifferentiation)*quality+ Model.productDifferentiation*sum0),2) ;
		
			equilProfit = marketProfit - activeLocationList.size()*locationCosts;
			
			/*total profit = equi profit - moving costs - exp. for inno and imi*/
			totalProfit = equilProfit - movingCosts ;
			
			
			if(Model.iteration>200)
				cumProfit += totalProfit;
		
		}
		
		
		/**
		 * This method determines the innovation process of the firm
		 */
		void innovationEffort(){

			/*Local innovation effort*/

			currentInnoProb = innovationProbability();
	
			if(Model.random() < currentInnoProb){

				// Add innovation into the pipeline
				pendingInnovationsList.add(new InnovationProcess(1,tempQualityConcept,true));
				
			}

		}

	/**
	 * This method determines the imitation process of the firm, where imitation is done in each location
	 */
	void imitationEffort(){
		
		
		/*Local imitation effort*/

		for(int i=0; i<activeLocationList.size();i++){

			FirmLocation aLocation = (FirmLocation) activeLocationList.get(i);
			
			currentImiProb = imitationProbability(aLocation.competitorList.size()+1);
	
			for(int j=0; j< aLocation.competitorList.size();j++){
	
				/*Random draw to determine if imitation is successful*/
				
				if(Model.random()< currentImiProb ){
				
						/*Imitation max pending concept quality of local competitor*/
						pendingInnovationsList.add(new InnovationProcess(activeLocationList.get(i).locationID,aLocation.competitorList.get(j).maxPendingQualityConcept,false));
		
				}
			}

		}
			
	}
		
		
	/**
	 * This function updates the quality level of the firm and the local knowledge stocks
	 */
	
	void updateQuality(){

		/*Update the qualities better than the firm's current one, that might become publicly available*/
		for(int i=0; i< competitorList.size(); i++){

			if(competitorList.get(i).qualityConcept> qualityConcept && competitorList.get(i).qualityConcept> maxPendingQualityConcept ){				
				
				pendingInnovationsList.add(new InnovationProcess(-1,competitorList.get(i).qualityConcept, false));
			}
		}
		
		/*Decrement the time to market */
		
		for(int i=0; i < pendingInnovationsList.size(); i++){
			
			pendingInnovationsList.get(i).timeToMarket--;
			
			/*If counter equals zero, then implement the new concept quality*/
			if(pendingInnovationsList.get(i).timeToMarket==0){
				
				tempQualityConcept = Math.max(tempQualityConcept,pendingInnovationsList.get(i).qualityConcect);
				pendingInnovationsList.remove(i);
				i--;
				
			}
		}

		//Update quality
		qualityConcept = tempQualityConcept;
		quality = qualityConcept;

		/*Update the number of locations*/
		numLocationsActive = activeLocationList.size();

	}
		
		
	// Some get functions	
	double getQuality(){
		
		return quality;
	}
	
	double getOutput(){
		
		return equilQuantity;
		
	}
	
	double getProfits(){
		
		return equilProfit;
		
	}
	
	
	int getNumLocations(){
		
		return numLocationsActive;
		
	}
	
	
		/**
		 * This function is used to compute the number of direct links to innovators and imitators. Used in the analysis
		 * 
		 */
		void computeLinks(){

			directLinksInnovators = 0;
			directLinksImitators = 0 ;
			
			for(int i=0; i < activeLocationList.size();i++){
				
				for(int j=0; j < activeLocationList.get(i).competitorList.size();j++){
					
					
					if(activeLocationList.get(i).competitorList.get(j).innovator){
						
						directLinksInnovators++;
						
					}else{
						
						directLinksImitators++;
			
						
					}

				}
	
			}
			
		}

		/**
		 * These functions is used in the GUI mode to draw the figures ()
		 * 
		 */

		@Override
		public void draw(SimGraphics g) {
		
			
		
			
			
			double value = Math.min((quality- Model.minQuality)/(Model.maxQuality-Model.minQuality),1.0);
			
	
			
			if(Double.isNaN(value)){
				
				value=1.0;
				
			}else if(Double.isInfinite(value)){
				
				value = 1.0;
			}
			
			int red, green, blue;
			
			if(value<0.5){
				
				
				red  = 255;
				green=(int) ((( value)/(0.5))*255);
				blue =0;
				
			}else{
				
				red  = 255;
				green =255;
				blue=(int) ((( value-0.5)/(1.0-0.5))*255);
				
				
			}
		
			
			Color color = new Color(red, green, blue);
		
			g.setFont(new Font("Arial", 16, 12));
			 
			 String label =Integer.toString(firmID);
			 
			 Rectangle2D r = g.getStringBounds(label);
		       // int height = (int)(r.getHeight()/g.getYScale()); 
		        //int width = (int)(r.getWidth()/g.getXScale());  
		      //  g.setDrawingCoordinates( getX(), getY(),0);
		      //  g.setDrawingParameters(width, height, 0);
		      
		        
		       if(innovator) 
		    	   g.drawStringInRect(color, Color.BLACK,  label);	
		       else
		    	   g.drawStringInOval(color, Color.BLACK,  label);	
		}
		
		



		@Override
		public int getX() {
			
			return ((int)(ModelGUI.dimX*0.5)) + x_coord;
		}


		@Override
		public int getY() {
			
			return ((int)(ModelGUI.dimY*0.5)) + y_coord;
		}	
		
		
		void setCoordinates(){
			
		
			
		}
		
		

		// Internal classes to hold the entry and exit characteristics. Those can be uesed to compute statistics 
			
			class EntryCharacteristics {
				
			
				int locationID;
				int numCompetitors;
				int numInnovators;
				int numImitators;
				int numLocations;
				
				int entry;
				int isInnovator;
				int marketEntry;
				int strategy;
				
				
				double averageQualityConcept;
				double lambda;
				double profits;
				double qualityProduct;
				double npv;
				double npvInclKappa;
				int switching;
				double npvBaseline;
				double npvGross;
				
			
				
				EntryCharacteristics(int locID, int numComp,  int doEntering, int marketEn, int sw , double np, double npka, double npvBase, double npGr){
					
				
					locationID = locID;
					
					/*Is innovator?*/
					if(innovator)
						isInnovator =1;
					else
						isInnovator=0;
					
					
					
					if(firmStratgy){

						strategy=1;
					}else{
						
						strategy=0;

					}
					
					/*Actually entered?*/
					entry=doEntering;
					switching = sw;
					/*Firm specific variables*/
					lambda = qualityConcept;
					profits = equilProfit;
					
					qualityProduct = quality;
					//Number locations before entry
					numLocations =  numLocationsActive;
					
					marketEntry =  marketEn;
					
					/*Location specific variables*/
					
					numCompetitors = numComp;
				
					averageQualityConcept = 0;
					numInnovators = 0;
					numImitators = 0;
					
					npv = np;
					npvInclKappa = npka;
					npvBaseline = npvBase;	
					npvGross = npGr;
			
					int counter =0;
		
					for(int i=0; i < Model.locationList.size();i++){
						
						if(Model.locationList.get(i).locationID==locationID){
							
							for(int j=0; j< Model.locationList.get(i).firmList.size();j++){
								
								for(int k=0; k<Model.locationList.get(i).firmList.get(j).activeLocationList.size();k++){
									
									if(Model.locationList.get(i).firmList.get(j).activeLocationList.get(k).locationID==locationID){
										
										
										averageQualityConcept += Model.locationList.get(i).firmList.get(j).qualityConcept;
										counter ++;
										
										if(Model.locationList.get(i).firmList.get(j).innovator)
											numInnovators++;
										else
											numImitators++;
										
									}
									
								}
								
								
								
							}
							
							
						}
						
						
						
					}
					
					if(counter >0 ){
						
						averageQualityConcept = averageQualityConcept/counter;
					}
					
				}
				
			}	
			
			
			
			class ExitCharacteristics {
				
				int locationID;
				int numLocations;
				double npv;
				double npvInclKappa;
				double lambda;
				double qualityProduct;
				int numCompetitors;
				
				double averageQualityConcept;
				int numInnovators;
				int numImitators;
				int exit;
				double profits;
			
				int isInnovator;
				int strategy;
				int switching;
				double npvBaseline;
				double npvGross;
				
				ExitCharacteristics(int locID, int numComp, int doExiting, int sw, double np, double npka, double npvBase, double npGr){
					
				
					locationID = locID;
					numLocations =  numLocationsActive;
					
					averageQualityConcept = 0;
					numCompetitors = numComp;
				
					int counter =0;
					numInnovators = 0;
					numImitators = 0;
					exit = doExiting;
					lambda = qualityConcept;
					profits = equilProfit;
					npv = np;
					npvInclKappa = npka;
					qualityProduct = quality;
					switching = sw;
					npvBaseline = npvBase;
					npvGross = npGr;
							
					if(innovator)
						isInnovator =1;
					else
						isInnovator=0;
					
					if(firmStratgy){

						strategy=1;
					}else{
						
						strategy=0;

					}
					
					for(int i=0; i < Model.locationList.size();i++){
						
						if(Model.locationList.get(i).locationID==locationID){
							
							for(int j=0; j< Model.locationList.get(i).firmList.size();j++){
								
								for(int k=0; k<Model.locationList.get(i).firmList.get(j).activeLocationList.size();k++){
									
									if(Model.locationList.get(i).firmList.get(j).activeLocationList.get(k).locationID==locationID){
										
										
										averageQualityConcept += Model.locationList.get(i).firmList.get(j).qualityConcept;
										counter ++;
										
										if(Model.locationList.get(i).firmList.get(j).innovator)
											numInnovators++;
										else
											numImitators++;
									}
									
								}
								
								
								
							}
							
							
						}
						
						
						
					}
					
					if(counter >0 ){
						
						averageQualityConcept = averageQualityConcept/counter;
					}
					
				}
				
			}
			
			
		
		
		class Strategy{
			
			boolean entry;
			boolean exit;
			
			int locationIDEntry;
			int locationIDExit;
			
			double npv;
			double npvIncKappa;
			double npvTp1;
			double npvBaseline;
			double npvGross;
			
			Strategy(int id, boolean en, boolean ex, double pv, double pvKa, double npBa, double npGr){
				
				if(en)
					locationIDEntry = id;
				else
					locationIDExit =id;
				
				entry = en;
				exit = ex;
				npv = pv;
				npvIncKappa = pvKa;
				npvBaseline = npBa;
				npvGross= npGr;
				
			}
			
			
			Strategy(int idEn,int idEx, boolean en, boolean ex, double pv, double pvKa, double npBa, double npGr){
				
				
					locationIDEntry = idEn;
			
					locationIDExit =idEx;
				
				entry = en;
				exit = ex;
				npv = pv;
				npvIncKappa = pvKa;
				npvBaseline = npBa;
				npvGross = npGr;
			}
			
			
		}
		

}
