import java.util.ArrayList;
import uchicago.src.sim.analysis.OpenSequenceGraph;
import uchicago.src.sim.analysis.Sequence;
import uchicago.src.sim.engine.SimInit;
import uchicago.src.sim.gui.CircularGraphLayout;
import uchicago.src.sim.gui.DisplaySurface;
import uchicago.src.sim.gui.Network2DDisplay;
import uchicago.src.sim.gui.Object2DDisplay;
import uchicago.src.sim.space.Object2DGrid;

/**This class is used for the GUI mode for running an illustrative single run. It opens a GUI which can be used to edit model parameters
 * and, after launching the simulations, opens some informative plots (time series of specific variables as well as simple plots 
 * illustrating the spatial interaction of firms)*/


public class ModelGUI extends SetGetModel{
	
	public OpenSequenceGraph graphOutput, graphPrice, graphQuality, graphFirmLocations,graphNumFirms,
	graphEnteringCosts, graphAverageProbability, 
	graphTotalProfit, graphSingleOutput, graphSinglePrice, graphSingleQuality, graphSingleFirmLocations,
	 graphSingleProbability, graphSinglelProfit, graphSingleQualityConcept, graphClusteringCoefficient, graphCumProfit;  
	
	
	
	ArrayList<Object2DGrid>  gridList ;
	ArrayList<DisplaySurface>  dsurfList ;
	ArrayList<GridSpaceItem> gridSpaceList;
	ArrayList<Coordinate> freeCoordinates;
	
	DisplaySurface surface;
	Network2DDisplay display;
	CircularGraphLayout layout;
	
	
	 public static int dimX, dimY;
	 
	
	/*Sequence class for total output*/
	class  SeqTotalOutput implements Sequence {
			
			
			public double getSValue(){
				
				return (double) (totalOutput);
				
			}
			
	}
	
	
	
	/*Sequence class for numActiveFirms*/
	class  SeqNumActiveFirms implements Sequence {
			
			
			public double getSValue(){
				
				return (double) (numActiveFirms);
				
			}
			
	}
	
	
	/*Sequence class for numInnovators*/
	class  SeqNumInnovators implements Sequence {
			
			
			public double getSValue(){
				
				return (double) (numInnovators);
				
			}
			
	}
	
	/*Sequence class for numImitators*/
	class  SeqNumImitators implements Sequence {
			
			
			public double getSValue(){
				
				return (double) (numImitators);
				
			}
			
	}
	
	
	
	
	/*Sequence class for price index*/
	class  SeqPriceIndex implements Sequence {
			
			
			public double getSValue(){
				
				return (double) (priceIndex);
				
			}
			
	}
	
	
	
	/*Sequence class for quality index*/
	class  SeqQualityIndex implements Sequence {
			
			
			public double getSValue(){
				
				return (double) (qualityIndex);
				
			}
			
	}
	
	
	
	/*Sequence class for clustering coefficient*/
	class  SeqClusteringCoefficient implements Sequence {
			
			
			public double getSValue(){
				
				return (double) (clusteringCoefficient);
				
			}
			
	}
	
	
	
	/*Sequence class for weighted clustering coefficient*/
	class  SeqWeightedClusteringCoefficient implements Sequence {
			
			
			public double getSValue(){
				
				return (double) (weightedClusteringCoefficient);
				
			}
			
	}
	
	/*Sequence class for average number of firm locations*/
	class SetAverageFirmLocations  implements Sequence {
		
		
		public double getSValue(){
			
			return (double) (averageFirmLocations);
			
		}
		
}
	


	
	/*Sequence class for sdOutput*/
	class SetSdOutput  implements Sequence {
		
		
		public double getSValue(){
			
			return (double) (sdOutput);
			
		}
		
}
	

	
	/*Sequence class for sdQuality*/
	class SetSdQuality  implements Sequence {
		
		
		public double getSValue(){
			
			return (double) (sdQuality);
			
		}
		
}
	
	

	
	/*Sequence class for sdAverageLocation*/
	class SetSdAverageLocation  implements Sequence {
		
		
		public double getSValue(){
			
			return (double) (sdAverageLocation);
			
		}
		
}
	
	
	

	
	/*Sequence class for averageInnoProbability*/
	class SetAverageInnoProbability  implements Sequence {
		
		
		public double getSValue(){
			
			return (double) (averageInnoProbability);
			
		}
		
}
	
	
	

	/*Sequence class for averageImiProbability*/
	class SetAverageImiProbability  implements Sequence {
		
		
		public double getSValue(){
			
			return (double) (averageImiProbability);
			
		}
		
}
	
	
	/*Sequence class for total profit*/
	class SetTotalProfit  implements Sequence {
		
		
		public double getSValue(){
			
			return (double) (totalProfit);
			
		}
		
}
	
	
	
	
	/*Sequence class for total profit*/
	class SetSinglePrice  implements Sequence {
		
int index, id;
		
SetSinglePrice(int i){	
			id = i;	
		}
		
		
		public double getSValue(){

			index = returnFirmIndex(id);
			
			if(index!=-1)
				return (double) (firmList.get(index).equilPrice);
			else
				return Double.NaN;
		}
		
	
		
}
	
	
	/*Sequence class for single profit*/
	class SetSingleProfit implements Sequence {
		
		int index, id;
		
		SetSingleProfit(int i){	
				id = i;	
			}
		
		
		public double getSValue(){

			index = returnFirmIndex(id);
			
			if(index!=-1)
				return (double) (firmList.get(index).equilProfit);
			else
				return Double.NaN;
		}
		
		
		
}
	
	
	/*Sequence class for output*/
	class SetSingleOutput  implements Sequence {
		
		int index, id;
		
		SetSingleOutput(int i){
			
			id = i;
			
		}
		
		
		public double getSValue(){
			
			index = returnFirmIndex(id);
			
			if(index!=-1)
				return (double) (firmList.get(index).equilQuantity);
			else
				return Double.NaN;
		}
		
}
	
	
	
	/*Sequence class for cumulated profit*/
	class SetCumProfit  implements Sequence {
		
		int index, id;
		

		SetCumProfit(int i){
			
			id = i;
			
		}
		
		
		public double getSValue(){
			
			index = returnFirmIndex(id);
			
			if(index!=-1)
				return (double) (firmList.get(index).cumProfit);
			else
				return Double.NaN;
		}
		
}
	
	
	
	/*Sequence class for single qualities*/
	class SetSingleQuality  implements Sequence {
		
		int index, id;
		
		SetSingleQuality(int i){	
			id = i;	
		}
		
		
		public double getSValue(){

			index = returnFirmIndex(id);
			
			if(index!=-1)
				return (double) (firmList.get(index).quality);
			else
				return Double.NaN;
		}
		
		
	
		
}
	
	
	/*Sequence class for inno probability*/
	class SetSingleInnoProbability  implements Sequence {
		
int index, id;
		
SetSingleInnoProbability(int i){	
			id = i;	
		}
		
		
		public double getSValue(){

			index = returnFirmIndex(id);
			
			if(index!=-1)
				return (double) (firmList.get(index).currentInnoProb);
			else
				return Double.NaN;
		}
		
		
}
	
	
	/*Sequence class for imi probs*/
	class SetSingleImiProbability  implements Sequence {
		
int index, id;
		
SetSingleImiProbability(int i){	
			id = i;	
		}
		
		
		public double getSValue(){

			index = returnFirmIndex(id);
			
			if(index!=-1)
				return (double) (firmList.get(index).currentImiProb);
			else
				return Double.NaN;
		}
		
		
}
	
	
	/*Sequence class for single number of locations*/
	class SetSingleNumLocations  implements Sequence {
		
int index, id;
		
SetSingleNumLocations(int i){	
			id = i;	
		}
		
		
		public double getSValue(){

			index = returnFirmIndex(id);
			
			if(index!=-1)
				return (double) (firmList.get(index).numLocationsActive);
			else
				return Double.NaN;
		}
		
		
}
	
	
	/*Sequence class for single qualities*/
	class SetSingleQualityConcept  implements Sequence {
		
int index, id;
		
SetSingleQualityConcept(int i){	
			id = i;	
		}
		
		
		public double getSValue(){

			index = returnFirmIndex(id);
			
			if(index!=-1)
				return (double) (firmList.get(index).qualityConcept);
			else
				return Double.NaN;
		}
		
		
		
}
	
	

	
	/*++++++++++++++++++++++++++  Methods ++++++++++++++++++++++++++++++++++++++++++*/
	
	/**
	 * setup method inherited from Model class
	 * */
	public void setup(){
		

		super.setup();
		
		params = new String[] {"numLocations","numFirms","timeHorizon","enteringCosts",
				 "reservationPrice", "productDifferentiation",
				 "adjustEnteringCosts",
				 "sdNoiseImitationKnowledge",
				 "sdNoiseImitationConcept",
				 "sigma","discontFactor","initialInnoEff","fractionEffectivityImitators",
				 "initialImiEff","exitingCosts","fractionEffectivityInnovators",
				 "locationCosts","isolatedLocation", "maxQualityConceptProgress","heterogeneityMode","randomAllocation", "numMontecarloRuns",
				 "maxQualityConceptProgress","rateLocationDecision", "marketEntryCosts",
				 "marketEntryHazardRate", "parTimeToMarket","ModelModeTimeToMarket","scrapValueParameterInno","scrapValueParameterImi","strategyParameter","targetNumFirms"};

	

		
		if(graphOutput!=null)
			graphOutput.dispose();
		
		if(graphPrice!=null)
			graphPrice.dispose();
		
		if(graphQuality!=null)
			graphQuality.dispose();
		if(graphNumFirms!=null)
			graphNumFirms.dispose();
		
		
		
		if(graphFirmLocations!=null)
			graphFirmLocations.dispose();
		
		if(graphEnteringCosts!=null)
			graphEnteringCosts.dispose();
		
		
		if(graphAverageProbability!=null)
			graphAverageProbability.dispose();
		
		if(graphTotalProfit!=null)
			graphTotalProfit.dispose();
		
		if(graphClusteringCoefficient!=null)
			graphClusteringCoefficient.dispose();
		
		
		if(graphSingleOutput!=null)
		graphSingleOutput.dispose();
		
		if(graphSinglePrice!=null)
		graphSinglePrice.dispose(); 
		
		if(graphSingleQuality!=null)
		graphSingleQuality.dispose(); 
		
		if(graphSingleFirmLocations!=null)
		graphSingleFirmLocations.dispose();
		
		if(graphSingleProbability!=null)
		 graphSingleProbability.dispose(); 
		
		if(graphSinglelProfit!=null)
		graphSinglelProfit.dispose();
		
		
		if(graphSingleQualityConcept!=null)
			graphSingleQualityConcept.dispose();
		
		
		if(graphCumProfit!=null)
			graphCumProfit.dispose();
		
		
	    
	    if(dsurfList!=null){
		    for(int i=0; i<dsurfList.size();i++){
		    		
		    		dsurfList.get(i).dispose();
		    
		    }
	}
	    
	    gridList = null;
	    dsurfList = null;
	    
	     gridList = new  ArrayList<Object2DGrid>() ;
	     dsurfList	= new	ArrayList<DisplaySurface>()  ;
	     
	     
	     dimX = ((int) Math.sqrt(numFirms))+5;
	     dimY=dimX;
	  
	     
	     
	     if (surface != null)
	         surface.dispose ();
	       surface = null;
	     
	     surface = new DisplaySurface (this, "Firm network");
	     registerDisplaySurface ("Firm network", surface);
		
		
	}
	

	/**
	 * Method: Build displays 
	 * */
	public void buildDisplay() {
		
		
		graphOutput = new OpenSequenceGraph("Total output", this);
		graphOutput.addSequence("Total Output", new SeqTotalOutput());
		graphOutput.display();
		graphOutput.step();
	
		graphPrice = new OpenSequenceGraph("Price index", this);
		graphPrice.addSequence("Price index", new SeqPriceIndex());
		graphPrice.display();
		graphPrice.step();
		
		graphQuality = new OpenSequenceGraph("Quality index", this);
		graphQuality.addSequence("Quality index", new SeqQualityIndex());
		graphQuality.display();
		graphQuality.step();
		
		graphNumFirms = new OpenSequenceGraph("Num Firms", this);
		graphNumFirms.addSequence("Num Firms", new SeqNumActiveFirms());
		graphNumFirms.addSequence("Num Imitators", new SeqNumImitators());
		graphNumFirms.addSequence("Num Innovators", new SeqNumInnovators());
		graphNumFirms.display();
		graphNumFirms.step();
		
		graphFirmLocations = new OpenSequenceGraph("Av. number locations", this);
		graphFirmLocations.addSequence("Average number locations", new SetAverageFirmLocations());
		graphFirmLocations.display();
		graphFirmLocations.step();
		
		graphClusteringCoefficient = new OpenSequenceGraph("Clustering Coefficient", this);
		graphClusteringCoefficient.addSequence("Clustering Coefficient", new SeqClusteringCoefficient());
		graphClusteringCoefficient.addSequence("Weighted Clustering Coefficient", new SeqWeightedClusteringCoefficient());
		graphClusteringCoefficient.display();
		graphClusteringCoefficient.step();
		
		graphSingleOutput = new OpenSequenceGraph("Single output", this);
		
		graphSinglePrice = new OpenSequenceGraph("Single price", this);
		graphSingleQuality = new OpenSequenceGraph("Single quality", this);
		graphSingleFirmLocations = new OpenSequenceGraph("Single num locations", this);
		graphSingleProbability = new OpenSequenceGraph("Single probabilities", this);
		
		graphSinglelProfit = new OpenSequenceGraph("Single profit", this);
		graphSingleQualityConcept =new OpenSequenceGraph("Single quality concept", this);
		graphCumProfit =	new OpenSequenceGraph("Cum Profit", this);
			
			
			
		for(int i=0; i < firmList.size();i++){

			graphSingleOutput.addSequence("Firm"+(i+1), new SetSingleOutput(firmList.get(i).firmID));
			
			graphSinglePrice.addSequence("Firm"+(i+1), new SetSinglePrice(firmList.get(i).firmID));
			graphSingleQuality.addSequence("Firm"+(i+1), new SetSingleQuality(firmList.get(i).firmID));
			graphSingleFirmLocations.addSequence("Firm"+(i+1), new SetSingleNumLocations(firmList.get(i).firmID));
		    graphSingleProbability.addSequence("Inno prob Firm"+(i+1), new SetSingleInnoProbability(firmList.get(i).firmID));
		    graphSingleProbability.addSequence("Imi Prob Firm"+(i+1), new SetSingleImiProbability(firmList.get(i).firmID));
		    graphSinglelProfit.addSequence("Firm"+(i+1), new SetSingleProfit(firmList.get(i).firmID));
		    graphSingleQualityConcept.addSequence("Firm"+(i+1), new SetSingleQualityConcept(firmList.get(i).firmID));
		    
		    graphCumProfit.addSequence("Firm"+(i+1), new SetCumProfit(firmList.get(i).firmID));
		    
			
		}
			
		graphSingleOutput.display();
		graphSingleOutput.step();
		

		graphSinglePrice.display();
		graphSinglePrice.step();	
		

		graphSingleQuality.display();
		graphSingleQuality.step();	
		
	
		graphSingleFirmLocations.display();
		graphSingleFirmLocations.step();
		

		graphSingleProbability.display();
		graphSingleProbability.step();
		

		graphSinglelProfit.display();
		graphSinglelProfit.step();
		
		graphSingleQualityConcept.display();
		graphSingleQualityConcept.step();
		
		graphCumProfit.display();
		graphCumProfit.step();

		graphAverageProbability = new OpenSequenceGraph("Prob inno - prob imi", this);
		graphAverageProbability.addSequence("Prob inno", new SetAverageInnoProbability());
		graphAverageProbability.addSequence("Prob imi", new SetAverageImiProbability());
		graphAverageProbability.setYRange(-0.2, 1.0);
		graphAverageProbability.display();
		graphAverageProbability.step();
		
		graphTotalProfit= new OpenSequenceGraph("Profits", this);
		graphTotalProfit.addSequence("Profits", new SetTotalProfit());
		graphTotalProfit.display();
		graphTotalProfit.step();
	
		layout = new CircularGraphLayout(nodeList, 100,100);
	
	 	display = new Network2DDisplay (layout);

	 
		surface.addDisplayableProbeable (display, "Jiggle View");
		surface.addZoomable (display);
		surface.setBackground (java.awt.Color.white);
		addSimEventListener (surface);

	     for(int i=0; i<numLocations;i++){
	    	 
	    	 gridList.add(new Object2DGrid(dimX,dimY));
	    	 dsurfList.add( new DisplaySurface(gridList.get(i).getSize(),this, "Location"+(i+1)));
	     }
			
		gridSpaceList =  new ArrayList<GridSpaceItem>();
		
		for(int i=0; i < locationList.size();i++){
			
			gridSpaceList.add(new GridSpaceItem(locationList.get(i)));
	
			 Object2DDisplay agentDisplay = new Object2DDisplay(gridList.get(i));
			 agentDisplay.setObjectList(gridSpaceList.get(gridSpaceList.size()-1).agentList);
			 
			 
			 dsurfList.get(i).addDisplayableProbeable(agentDisplay, "Agents");
			 addSimEventListener(dsurfList.get(i));
		
		}
		
		
		int fi=0;
		
		for(int j =-1; j<=1;j++){
			
			for(int k =-1; k<=1;k++){
				
				
				if(!(k==0 && j ==0) ){
					firmList.get(fi).x_coord = k;
					firmList.get(fi).y_coord = j;
					
					fi++;
					
					for(int l=0; l < freeCoordinates.size();l++){
						
						if(freeCoordinates.get(l).xCo==k && freeCoordinates.get(l).yCo==j){
							
							freeCoordinates.remove(l);
							l--;
							break;
			
						}
					
					}
				
					if(fi==firmList.size()){
						
						break;
					}
				
				}
			
			}
			
			if(fi==firmList.size()){
				
				break;
			}
		}
   
	}
	
	
	
	/**
	 * This adds a new time series sequence to a graphic*/
	public void addSequencesToGraphs(){
		
		
		graphSingleOutput.addSequence("Firm"+firmList.get(firmList.size()-1).firmID, new SetSingleOutput(firmList.get(firmList.size()-1).firmID));
		
		graphSinglePrice.addSequence("Firm"+firmList.get(firmList.size()-1).firmID, new SetSinglePrice(firmList.get(firmList.size()-1).firmID));
		graphSingleQuality.addSequence("Firm"+firmList.get(firmList.size()-1).firmID, new SetSingleQuality(firmList.get(firmList.size()-1).firmID));
		graphSingleFirmLocations.addSequence("Firm"+firmList.get(firmList.size()-1).firmID, new SetSingleNumLocations(firmList.get(firmList.size()-1).firmID));
	    graphSingleProbability.addSequence("Inno prob Firm"+firmList.get(firmList.size()-1).firmID, new SetSingleInnoProbability(firmList.get(firmList.size()-1).firmID));
	    graphSingleProbability.addSequence("Imi Prob Firm"+firmList.get(firmList.size()-1).firmID, new SetSingleImiProbability(firmList.get(firmList.size()-1).firmID));
	    graphSinglelProfit.addSequence("Firm"+firmList.get(firmList.size()-1).firmID, new SetSingleProfit(firmList.get(firmList.size()-1).firmID));
	    graphSingleQualityConcept.addSequence("Firm"+firmList.get(firmList.size()-1).firmID, new SetSingleQualityConcept(firmList.get(firmList.size()-1).firmID));
	    graphCumProfit.addSequence("Firm"+firmList.get(firmList.size()-1).firmID, new SetCumProfit(firmList.get(firmList.size()-1).firmID));
	   
	    firmList.get(firmList.size()-1).x_coord = freeCoordinates.get(0).xCo;
	    firmList.get(firmList.size()-1).y_coord = freeCoordinates.get(0).yCo;
	    freeCoordinates.remove(0);
	   
		
	}
	
	/**
	 * Method executed at the beginning
	 * */
	 public void begin(){
		 
		
		 freeCoordinates = new ArrayList<Coordinate>();
		 

			for(int i = 1; i <25; i++){
				for(int j =-2; j<=2;j++){
					
					for(int k =-2; k<=2;k++){
						
						if(!(k==0 && j ==0) ){
							
							
						
							freeCoordinates.add(new Coordinate(k,j));
							
						
						}

						
					}
					
					
				}

				
			}
			
		    buildModel();
		    buildSchedule();
		    buildDisplay();
		    
		    for(int i=0; i < dsurfList.size();i++){
		    	dsurfList.get(i).display();
		    }
		    
		   layout.updateLayout();
		    surface.display ();
		    
		  }
	
	
	/**
	 * Build model method inherited from Model class
	 * */
	public void buildModel(){
		
		super.buildModel();
		
	}
	
	
	/**
	 * step method inherited from Model class
	 * */
	public void step(){
		
		
		super.step();
		
		
		
		if(marketEntryHappend){
			
			
			addSequencesToGraphs();
			
			
			
		}
		
		
		if(marketExitHappend){
			
			for(int i =0; i < exitedFirms.size();i++){
			
				freeCoordinates.add(new Coordinate(exitedFirms.get(i).x_coord,exitedFirms.get(i).y_coord)); 
			
			}
			
		}
	
		graphOutput.step();
		graphPrice.step();
		graphQuality.step();
		graphFirmLocations.step();
		//graphEnteringCosts.step();
		graphAverageProbability.step();
		graphTotalProfit.step();
		graphNumFirms.step();
		
		graphSingleOutput.step();
		graphSinglePrice.step(); 
		graphSingleQuality.step(); 
		graphSingleFirmLocations.step();
	    graphSingleProbability.step(); 
		graphSinglelProfit.step();
		graphSingleQualityConcept.step();
	
		graphClusteringCoefficient.step();
		
		graphCumProfit.step();
		
		for(int i=0; i < dsurfList.size();i++){
			
			
			gridSpaceList.get(i).update(locationList.get(i));;
			
		    dsurfList.get(i).updateDisplay();
		    }
		
			layout.setList(nodeList);
		
		 layout.updateLayout();
		 surface.updateDisplay();
		 
		
		
	}
	
/*
 * Get and set methods for model parameters
 * */
	
	
	public double getExitingCosts(){
		
		return exitingCosts;
		
	}
	
	public void setExitingCosts( double exco){
		
		exitingCosts = exco;
		
	}
	

	public boolean getModelModeTimeToMarket(){
		
		return modelModeTimeToMarket;
		
	}
	
	public void setModelModeTimeToMarket( boolean exco){
		
		modelModeTimeToMarket = exco;
		
	}
	
	
	public void setMaxQualityConceptProgress(double ma){
		
		
		maxQualityConceptProgress = ma;
		
	}
	
	
	public double getMaxQualityConceptProgress(){
		
		return maxQualityConceptProgress;
		
		
	}
	
	public double getReservationPrice(){
			
		 return reservationPrice;
			
			
	}
	 
	public void setReservationPrice(double resPr){
			
		 reservationPrice = resPr;
			
			
	}
	
	public double getProductDifferentiation(){


		return productDifferentiation;
		
	}
	
	public void setProductDifferentiation(double proDiff){
			
		productDifferentiation = proDiff;
			
	}
	
	public void setMarketEntryHazardRate(double hr){
		
		
		marketEntryHazardRate = hr;
		
	}
	
	public double getMarketEntryHazardRate(){
		
		
		return marketEntryHazardRate;
		
	}
	

public void setMarketEntryCosts(double hr){
		
		
	marketEntryCosts = hr;
		
	}
	
	public double getMarketEntryCosts(){
		
		
		return marketEntryCosts;
		
	}
	
	
 
	public double getSdNoiseImitationKnowledge(){
			
			return sdNoiseImitationKnowledge;
			
	}
	 
	 
	public void setSdNoiseImitationKnowledge(double sdIm){
			
			sdNoiseImitationKnowledge = sdIm;
			
	}
	 
	public double getSdNoiseImitationConcept(){
			
			return sdNoiseImitationConcept;
			
	}
	 
	 
	public void setSdNoiseImitationConcept(double sdIm){
			
			sdNoiseImitationConcept = sdIm;
			
	}
	 
	 
	public double getSigma(){
			
			return sigma;
			
	}
	 
	public void setSigma(double be){
			
		sigma = be;
			
	}
	

	public double getDiscontFactor(){
			
			return discontFactor;
			
	}
	 
	public void setDiscontFactor(double dr){
			
			 discontFactor = dr;
			
	}
	
	
public double getInitialImiEff(){
		
		return initialImiEff;
		
}
 
	public void setInitialImiEff(double dr){
		
		initialImiEff = dr;
		
}
	
	
	public double getInitialInnoEff(){
		
		return initialInnoEff;
		
}
 
	public void setInitialInnoEff(double dr){
		
		initialInnoEff = dr;
		
}

public double getExogenousQualityConceptProgress(){
		
		return maxQualityConceptProgress;
		
}
 
	public void setExogenousQualityConceptProgress(double dr){
		
		maxQualityConceptProgress = dr;
		
}
	

	/**
	 * *** This is the main method ***
	 * */
	public static void main (String[] args){
		
		
		SimInit init = new SimInit();
		ModelGUI m = new ModelGUI();
		init.loadModel(m, null, false);
		
	}
	
	
	
	class GridSpaceItem{
		
		ArrayList<Object> agentList ;
		Location location;
		
		GridSpaceItem(Location loc){
			
			location = loc;
			
			agentList = new ArrayList<Object>();
			
			agentList.add(location);
			agentList.addAll(location.firmList);
			
			
		}
		
		void update(Location loc){
			
			agentList.clear();
			agentList.add(location);
			agentList.addAll(location.firmList);
			
		};
		
		
		
	}
	
	
	class Coordinate{
		
		int xCo, yCo;
		
		Coordinate(int x, int y){
			
			xCo=x;
			yCo=y;
			
		}
		
	}

}
