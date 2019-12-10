/**
 * This class holds some set and get methods (used for manipulating model parameters), that are used at several
 * places of the code.
 * */
public class SetGetModel extends Model{

	/*++++++++Get and set methods for parameters++++++++*/
	public int getNumLocations(){
	
		return numLocations;
	
	}
	
	public void setNumLocations(int numloc){
		
		numLocations = numloc;
		
	}
	
	

	
public double getFractionEffectivityImitators(){
		
		return fractionEffectivityImitators;
		
	}
	
	public void setFractionEffectivityImitators(double frac){
			
		fractionEffectivityImitators = frac;
			
	}

public double getFractionEffectivityInnovators(){
		
		return fractionEffectivityInnovators;
		
	}
	
	public void setFractionEffectivityInnovators(double frac){
			
		fractionEffectivityInnovators = frac;
			
	}

	public int getNumFirms(){
	
		return numFirms;
	
	}
	
	public void setNumFirms(int nofirms){
		
		numFirms = nofirms;
		
	}
	
	public int getTimeHorizon(){
		
		return timeHorizon;
		
	}
	
	public void setTimeHorizon(int tiho){
		
	 timeHorizon = tiho;
		
	}
	
	
	public int getNumMontecarloRuns(){
		
		return numMontecarloRuns;
		
	}
	
	public void setNumMontecarloRuns(int tiho){
		
		numMontecarloRuns = tiho;
		
	}
	
public double getRateLocationDecision(){
		
		
		return rateLocationDecision;
	}
	
	
	public void setRateLocationDecision(double f){
		
		
		rateLocationDecision = f;
	}
	
	

	public double getSigma(){
			
			return sigma;
			
	}
	 
	public void setSigma(double be){
			
		sigma = be;
			
	}
	
	
	
	public double getEnteringCosts(){
		
		return enteringCosts;
		
	}
	
	public void setEnteringCosts( double enco){
		
		 enteringCosts = enco;
		
	}

public void setMarketEntryCosts(double hr){
		
		
	marketEntryCosts = hr;
		
	}
	
	public double getMarketEntryCosts(){
		
		
		return marketEntryCosts;
		
	}

public double getLocationCosts(){
		
		return locationCosts;
		
	}
	
	public void setLocationCosts( double loco){
		
		locationCosts = loco;
		
	}

public int getParTimeToMarket(){
		
		return parTimeToMarket;
		
	}
	
	public void setParTimeToMarket( int loco){
		
		parTimeToMarket = loco;
		
	}
	
	
public double getStrategyParameter(){
		
		return strategyParameter;
		
	}
	
	public void setStrategyParameter( double loco){
		
		strategyParameter = loco;
		
	}

}
